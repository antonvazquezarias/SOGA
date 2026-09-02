import pandas as pd
from collections import defaultdict
import re
import subprocess
from pathlib import Path
import tempfile
import shutil
from Bio import SeqIO
from BCBio import GFF

# Functions
def fasta_contig_ids(fasta_file, directory: str = ""):
    """
    Returns the sequence IDs of a FASTA file, in order of appearance.

    The ID is the first whitespace-delimited token of each header line, which
    is how Aragorn refers to each contig. MFannot instead renames them C_0,
    C_1, ... in the same order, so this list maps one onto the other.
    """
    ids = []
    with open(Path(directory) / fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids

def reorder_fasta(fasta_file, strand, cut_pos, new_seq_name=None, old_seq_name=None, directory: str = ""):
    """
    Reorders a circular sequence starting from a specific cut position
    and optionally adjusts the strand.

    Parameters:
    -----------
    fasta_file : str
        Name of the input FASTA file.
    strand : str
        Target strand ('+' or '-'). If '-', the sequence is reverse complemented.
    cut_pos : int
        The 0-based index at which to cut the sequence. The new sequence will start
        from this position.
    new_seq_name : str, optional
        Name for the output file. If None, it overwrites the input path (unless changed elsewhere).
    old_seq_name : str, optional
        Name to save the original sequence as backup.
    directory : str or Path
        Directory containing the files.
    """
    fasta_path = Path(directory, fasta_file)
    
    with open(fasta_path) as f:
        header = f.readline().strip()
        seq = "".join(line.strip() for line in f)

    seq = seq[cut_pos:] + seq[:cut_pos]

    if strand == "-":
        complement = str.maketrans("ACGTacgt", "TGCAtgca")
        seq = seq.translate(complement)[::-1]
    
    if old_seq_name:
        old_seq_path = Path(directory, old_seq_name)
        fasta_path.rename(old_seq_path)

    if new_seq_name:
        fasta_path = Path(directory, new_seq_name)

    with open(fasta_path, "w") as f:
        f.write(header + "\n")
        f.write(seq + "\n")
        
# HMM models shipped inside the MFannot image, used to locate rns cheaply
RNS_MODELS = [
    "/MFannot_data/models/HMM_models/RNA/5-rns-other.hmm",
    "/MFannot_data/models/HMM_models/RNA/5-rns-fungi.hmm",
]

def find_rns_start(fasta_file, directory: str = "", evalue="1e-3"):
    """
    Locates the 5' end of rns with nhmmer, without running a full annotation.

    Unlike MFannot, this also finds an rns that spans the origin of the linear
    representation: MFannot discards a feature whose parts are not in
    increasing coordinate order, whereas nhmmer scores each match on its own.

    Parameters:
    -----------
    fasta_file : str
        Single-sequence FASTA file to search.
    directory : str or Path
        Directory containing the file.
    evalue : str
        Reporting threshold passed to nhmmer.

    Returns:
    --------
    tuple or None
        (strand, position) where position is the 1-based coordinate of the
        first base of rns on that strand, or None if no match was found.
    """
    path = Path(directory) / fasta_file

    length = 0
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                length += len(line.strip())

    best = None
    for model in RNS_MODELS:
        if not Path(model).exists():
            continue
        proc = subprocess.run(
            ["nhmmer", "--dna", "--tblout", "/dev/stdout", "-o", "/dev/null",
             "-E", evalue, model, str(path)],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            universal_newlines=True)
        if proc.returncode != 0:
            continue
        for line in (proc.stdout or "").splitlines():
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            if len(fields) < 13:
                continue
            hmmfrom, alifrom, strand = int(fields[4]), int(fields[6]), fields[11]
            score = float(fields[12])
            if best is None or score < best[0]:
                best = (score, hmmfrom, alifrom, strand)

    if best is None:
        return None

    _, hmmfrom, alifrom, strand = best
    # The match may start inside the model, so walk back to model position 1.
    offset = hmmfrom - 1
    start = alifrom - offset if strand == "+" else alifrom + offset
    # Walking back can run off either end of a circular sequence.
    start = ((start - 1) % length) + 1
    return strand, start

def annotate_mfannot(fasta_file, directory, organelle):
    """
    Runs MFannot using a Docker container to annotate the sequence.

    Parameters:
    -----------
    fasta_file : str
        Input FASTA file.
    directory : str or Path
        Directory containing the file (will be mounted to Docker).
    organelle : str
        Type of organelle ('chloroplast' or 'mitochondrion') to determine genetic code.
    """
    folder = Path(directory).resolve()
    file_name = Path(fasta_file).stem

    if organelle.lower() == "chloroplast":
        genetic_code = 11
    elif organelle.lower() in ["mitochondrion", "mitochondria"]:
        genetic_code = 4
    else:
        genetic_code = None
    
    # Run MFannot
    proc = subprocess.run([
        "sh", "-c",
        f"cd {folder} && /mfannot/mfannot -g {genetic_code} --tbl {file_name}.fasta"
    ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)

    candidates = sorted(folder.glob(f"{file_name}.fasta.new*.tbl"))
    if proc.returncode != 0 or not candidates:
        log = (proc.stdout or "").strip()
        tail = "\n".join(log.splitlines()[-20:]) or "(MFannot produced no output)"
        raise RuntimeError(
            f"MFannot failed for {file_name} in {folder}\n"
            f"  exit code: {proc.returncode}\n"
            f"  .tbl files produced: {len(candidates)}\n"
            f"  last lines of MFannot output:\n{tail}"
        )
    tbl_file = candidates[-1]
    shutil.move(tbl_file, folder / f"{file_name}.tbl")
    return proc.stdout or ""

def annotate_aragorn(fasta_file, directory, organelle, circular=True):
    """
    Runs Aragorn using a Docker container to identify tRNAs and tmRNAs.

    Parameters:
    -----------
    fasta_file : str
        Input FASTA file.
    directory : str or Path
        Directory containing the file.
    organelle : str
        Organelle type for genetic code selection.
    circular : bool
        Topology of the sequence. True for circular, False for linear.
    """
    folder = Path(directory).resolve()
    file_name = Path(fasta_file).stem

    if organelle.lower() == "chloroplast":
        genetic_code = 11
    elif organelle.lower() in ["mitochondrion", "mitochondria"]:
        genetic_code = 4
    else:
        genetic_code = None
    
    if circular:
        topology = "c"
    else:
        topology = "l"  

    # Run Aragorn in batch mode
    proc = subprocess.run([
        "aragorn",
        "-t",       # tRNA
        "-m",       # tmRNA
        "-i",       # introns
        f"-{topology}",
        f"-g{genetic_code}",
        "-w",       # batch mode
        "-o", f"{folder}/{file_name}.txt",
        f"{folder}/{file_name}.fasta"
    ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)

    if proc.returncode != 0:
        log = (proc.stdout or "").strip()
        tail = "\n".join(log.splitlines()[-20:]) or "(Aragorn produced no output)"
        raise RuntimeError(
            f"Aragorn failed for {file_name} in {folder}\n"
            f"  exit code: {proc.returncode}\n"
            f"  last lines of Aragorn output:\n{tail}"
        )
    return proc.stdout or ""

def mfannot_to_gff3(tbl_file, seq_name="Default_name", export=True, directory: str = "", organelle="", contig_ids=None):
    """
    Parses MFannot .tbl output and converts it to GFF3 format.

    Parameters:
    -----------
    tbl_file : str
        Name of the MFannot .tbl file.
    seq_name : str
        Fallback sequence identifier, used when contig_ids is not given.
    export : bool
        If True, writes the result to a .gff3 file.
    directory : str or Path
        Working directory.
    organelle : str
        Used to set translation table attributes.
    contig_ids : list of str, optional
        Real sequence IDs in FASTA order. MFannot renames contigs to C_0, C_1,
        ... in that order, so the Nth ">Feature" section is given
        contig_ids[N]. If None, seq_name is used for every feature.

    Returns:
    --------
    pd.DataFrame
        DataFrame containing the GFF3 data.
    """

    gff_columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    # A table holds one section per contig, each introduced by a
    # ">Feature C_<n> ..." line in FASTA order.
    with open(Path(directory) / tbl_file) as f:
        raw_lines = f.read().splitlines()

    lines = [l for l in raw_lines if l.strip() and not l.startswith(">")]

    if len(lines) == 0:
        empty_df = pd.DataFrame(columns=gff_columns)
        
        if export:
            file_name = Path(directory) / f"{seq_name}_MF.gff3"
            with open(file_name, "w") as f:
                f.write("##gff-version 3\n")
            # This ensures the file exists and is valid GFF3 even if empty
            empty_df.to_csv(file_name, sep="\t", header=False, index=False, mode="a")
            
        return empty_df


    # Parse features and their qualifiers, tracking which contig each belongs to
    features = []
    current_feature = None
    contig_index = -1

    for line in raw_lines:
        if not line.strip():
            continue
        if line.startswith(">"):
            # ">Feature C_<n> Table1": trust the number MFannot writes rather
            # than counting sections, so a missing section cannot shift the map.
            match = re.match(r">Feature\s+C_(\d+)", line)
            contig_index = int(match.group(1)) if match else contig_index + 1
            current_feature = None
            continue
        elements = line.split("\t")
        if elements[0]:
            current_feature = {"feature": elements, "qualifiers": [],
                               "contig_index": contig_index}
            features.append(current_feature)
        elif current_feature is not None:
            current_feature["qualifiers"].append(elements)

    rows = []
    ID_count = defaultdict(int)
    gene_ID = feature_ID = None
    cumulative_CDS_length = 0

    for feat in features:
        f = feat["feature"]
        qual = feat["qualifiers"]

        # Map this feature onto its real sequence ID via the contig section order
        ci = feat["contig_index"]
        seqid = contig_ids[ci] if (contig_ids and 0 <= ci < len(contig_ids)) else seq_name

        # Safe extraction of row_type (inherit from previous if missing)
        if len(f) > 2 and f[2]:
            row_type = f[2]
        elif qual and len(qual[0]) >= 5 and "ycf" in qual[1][4]:
            row_type = "CDS"
            qual[0][4] = "hypothetical protein"
        elif qual and len(qual[0]) >= 5 and "rnz" in qual[1][4]:
            row_type = "CDS"
            qual[0][4] = "Ribonuclease Z"
        elif qual and len(qual[0]) >= 5 and "odp" in qual[1][4].lower():
            row_type = "CDS"
            if "odpa" in qual[1][4].lower():
                qual[0][4] = "Pyruvate dehydrogenase E1 component subunit alpha"
            elif "odpb" in qual[1][4].lower():
                qual[0][4] = "Pyruvate dehydrogenase E1 component subunit beta"
        else:
            row_type = rows[-1][2] if rows else ""
        
        # Determine start/end and strand based on coordinate order
        start, end = (int(f[0]), int(f[1])) if len(f) > 1 and int(f[0]) < int(f[1]) else (int(f[1]), int(f[0]))
        strand = "+" if int(f[0]) < int(f[1]) else "-"

        qualifiers_dict = {}

        # Construct IDs and Names based on feature type
        if row_type == "gene":
            gene_ID = qual[0][4] if qual else "gene_unknown"
            ID = gene_ID
            qualifiers_dict["name"] = f"{gene_ID} gene"
            cumulative_CDS_length = 0
        elif row_type not in ["gene", "intron", "exon"]:
            feature_ID = f"{gene_ID}_{row_type}"
            ID = feature_ID
            qualifiers_dict["name"] = f"{gene_ID} {row_type}"
            qualifiers_dict["parent"] = gene_ID
        else:
            ID = f"{gene_ID}_{row_type}"
            qualifiers_dict["name"] = f"{gene_ID} {row_type}"
            qualifiers_dict["parent"] = feature_ID
        
        if row_type == "RNA":
            row_type = "ncRNA"

        # Handle duplicate feature IDs
        ID_count[ID] += 1
        if ID_count[ID] > 1:
            qualifiers_dict["ID"] = f"{ID}{ID_count[ID]}"
            if row_type == "gene":
                gene_ID = f"{ID}{ID_count[ID]}"
            elif row_type not in ["gene", "intron", "exon"]:
                feature_ID = f"{ID}{ID_count[ID]}"
        else:
            qualifiers_dict["ID"] = ID

        # Calculate Phase for CDS
        phase = "."
        if row_type == "CDS":
            phase = (3 - (cumulative_CDS_length % 3)) % 3
            cumulative_CDS_length += end - start + 1
            
            if organelle.lower() == "chloroplast":
                qualifiers_dict["transl_table"] = 11
            elif organelle.lower() in ["mitochondrion", "mitochondria"]:
                qualifiers_dict["transl_table"] = 4

        # Parse original MFannot qualifiers
        for q in qual:
            if len(q) < 5:
                continue

            key = q[3]
            val = q[4]

            # Skip gene names containing "orf" or Clean protein IDs
            if key == "gene" and "orf" in val:
                continue
            if key == "protein_id":
                qualifiers_dict[key] = val.replace("lcl| ","")
                continue

            qualifiers_dict[key] = val

        row = [
            seqid,
            "MFannot",
            row_type,
            start,
            end,
            ".",
            strand,
            phase,
            ";".join(f"{k}={v}" for k, v in qualifiers_dict.items())
        ]
        rows.append(row)

    df = pd.DataFrame(rows)
    df.columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    df = sort_gff3(df)

    if export:
        file_name = f"{directory}{seq_name}_MF.gff3"
        with open(file_name, "w") as f:
            f.write("##gff-version 3\n")
        df.to_csv(file_name, sep="\t", header=False, index=False, mode="a")

    return df

def aragorn_to_gff3(txt_file, seq_name="Default_name", export=True, directory: str = "", organelle=""):
    """
    Parses Aragorn output text file and converts it to GFF3 format.

    Parameters:
    -----------
    txt_file : str
        Name of the Aragorn output file.
    seq_name : str
        Fallback sequence identifier, used only if a data row precedes any
        contig header, which does not happen in normal Aragorn output.
    export : bool
        If True, writes the result to a .gff3 file.
    directory : str or Path
        Working directory.
    organelle : str
        Used to set translation table attributes for tmRNA.

    Returns:
    --------
    pd.DataFrame
        DataFrame containing the GFF3 data.
    """
    gff_columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

    # Batch (-w) output holds one section per contig:
    #   >contig_id len=...
    #   N genes found
    #   1  tRNA-Xxx  [a,b]  pos  (codon)
    # and closes with a ">end ..." line. Keep the contig id for every data row,
    # since a contig with no genes must not discard the others.
    with open(Path(directory) / txt_file) as f:
        raw_lines = f.read().splitlines()

    seqids = []
    data_rows = []
    current_seqid = seq_name
    for line in raw_lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if not stripped.startswith(">end"):
                current_seqid = stripped[1:].split()[0]
            continue
        if stripped.endswith("genes found"):
            continue
        seqids.append(current_seqid)
        data_rows.append(re.sub(r"[ \t]+", " ", stripped).split(" "))

    if not data_rows:
        empty_df = pd.DataFrame(columns=gff_columns)

        if export:
            file_name = Path(directory) / f"{seq_name}_AR.gff3"
            with open(file_name, "w") as f:
                f.write("##gff-version 3\n")
            # This ensures the file exists and is valid GFF3 even if empty
            empty_df.to_csv(file_name, sep="\t", header=False, index=False, mode="a")

        return empty_df

    df = pd.DataFrame(data_rows)
    df['seqid'] = seqids

    # Cleaning and Pre-processing columns. Column 0 is Aragorn's per-section
    # row number, which restarts on every contig.
    df = df.drop(columns=0).reset_index(drop=True)
    df['strand'] = df[2].str.startswith('c').map({True: '-', False: '+'})
    
    start_end = df[2].str.extract(r'c?\[(\d+),(\d+)\]')
    df['start'] = pd.to_numeric(start_end[0], errors='coerce')
    df['end'] = pd.to_numeric(start_end[1], errors='coerce')
    
    intron_match = df[4].str.extract(r'i\((\d+),(\d+)\)')
    df['intron_distance'] = pd.to_numeric(intron_match[0], errors='coerce')
    df['intron_length'] = pd.to_numeric(intron_match[1], errors='coerce')
    
    df[4] = df[4].str.replace(r'i\(\d+,\d+\)', '', regex=True).str.strip()
    
    type_aa = df[1].str.split('-', n=1, expand=True)
    df['type'] = type_aa[0]
    df['aminoacid'] = type_aa[1]
    
    df = df.rename(columns={4: 'codon', 3: 'position'})
    
    # Calculate relative start/end for tmRNA and tRNA features
    df['feature_start'] = pd.NA
    df['feature_end'] = pd.NA

    for idx, row in df.iterrows():
        raw = row['position']
        if row['type'] == 'tmRNA':
            rel_start, rel_end = map(int, str(raw).split(','))
            if row['strand'] == '+':
                df.at[idx, 'feature_start'] = row['start'] + rel_start - 1
                df.at[idx, 'feature_end']   = row['start'] + rel_end   - 1
            else:
                df.at[idx, 'feature_start'] = row['end'] - rel_end   + 1
                df.at[idx, 'feature_end']   = row['end'] - rel_start + 1
        elif row['type'] == 'tRNA':
            raw = int(raw)
            if row['strand'] == '+':
                df.at[idx, 'feature_start'] = row['start'] + raw - 1
                df.at[idx, 'feature_end']   = row['start'] + raw + 1
            else:
                df.at[idx, 'feature_start']   = row['end'] - raw - 1
                df.at[idx, 'feature_end'] = row['end'] - raw + 1
            
    # Compute absolute intron coordinates
    df['intron_start'] = pd.NA
    df['intron_end'] = pd.NA

    for idx, row in df.iterrows():
        if pd.notna(row['intron_distance']):
            if row['strand'] == '+':
                df.at[idx, 'intron_start'] = row['start'] + row['intron_distance'] - 1
                df.at[idx, 'intron_end'] = row['start'] + row['intron_distance'] + row['intron_length'] - 2
            else:
                df.at[idx, 'intron_end'] = row['end'] - row['intron_distance'] + 1
                df.at[idx, 'intron_start'] = row['end'] - row['intron_distance'] - row['intron_length'] + 2
    
    numeric_cols = ['start', 'end', 'position', 'intron_distance', 'intron_length']
    df[numeric_cols] = df[numeric_cols].apply(pd.to_numeric, errors='coerce')

    df = df[['seqid', 'type', 'start', 'end', 'strand', 'codon', 'aminoacid', 'feature_start', 'feature_end', 'intron_start', 'intron_end']]

    # Construct GFF3 rows
    gff_rows = []
    name_counter = defaultdict(int)

    aminoacids = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Sec": "U", "Pyl": "O"
    }

    for idx, row in df.iterrows():
        # Gene logic
        if row['type'] == "tRNA":
            gene_base = f"trn{aminoacids.get(str(row['aminoacid']), '')}{row['codon']}"
        elif row['type'] == "tmRNA":
            gene_base = "ssrA"
        else:
            gene_base = row['type']

        name_counter[gene_base] += 1
        gene_id = gene_base if name_counter[gene_base] == 1 else f"{gene_base}{name_counter[gene_base]}"
        tmRNA_id = f"{gene_id}_tmRNA"

        gff_rows.append({
            "seqid": row['seqid'],
            "source": "ARAGORN",
            "type": "gene",
            "start": row['start'],
            "end": row['end'],
            "score": ".",
            "strand": row['strand'],
            "phase": ".",
            "attributes": f"ID={gene_id};name={gene_base};gene={gene_base}"
        })

        # tRNA Logic (handling introns)
        if row['type'] == "tRNA":
            position = f"complement({row['feature_start']}..{row['feature_end']})" if row['strand'] == "-" else f"{row['feature_start']}..{row['feature_end']}"
            anticodon_attr = f"anticodon=(pos:{position},aa:{row['aminoacid']},seq:{str(row['codon']).replace('(', '').replace(')', '')})"
            
            if pd.notna(row['intron_start']):
                # Split tRNA with intron
                anticodon_attr = "" # Simplified for split features
                if row['strand'] == '-':
                    anticodon_attr = ""

                # Part 1 (Exon 1)
                gff_rows.append({
                    "seqid": row['seqid'],
                    "source": "ARAGORN",
                    "type": "tRNA",
                    "start": row['start'],
                    "end": row['intron_start'] - 1,
                    "score": ".",
                    "strand": row['strand'],
                    "phase": ".",
                    "attributes": ";".join([x for x in [
                        f"ID={gene_id}_tRNA",
                        f"name={gene_id} tRNA",
                        f"parent={gene_id}",
                        f"product=tRNA-{row['aminoacid']}",
                        anticodon_attr
                    ] if x])
                })
                gff_rows.append({
                    "seqid": row['seqid'],
                    "source": "ARAGORN",
                    "type": "exon",
                    "start": row['start'],
                    "end": row['intron_start'] - 1,
                    "score": ".",
                    "strand": row['strand'],
                    "phase": ".",
                    "attributes": f"ID={gene_id}_exon1;name={gene_id} exon1;parent={gene_id}_tRNA"
                })
                # Intron
                gff_rows.append({
                    "seqid": row['seqid'],
                    "source": "ARAGORN",
                    "type": "intron",
                    "start": row['intron_start'],
                    "end": row['intron_end'],
                    "score": ".",
                    "strand": row['strand'],
                    "phase": ".",
                    "attributes": f"ID={gene_id}_intron;name={gene_id} intron;parent={gene_id}"
                })
                # Part 2 (Exon 2)
                gff_rows.append({
                    "seqid": row['seqid'],
                    "source": "ARAGORN",
                    "type": "tRNA",
                    "start": row['intron_end'] + 1,
                    "end": row['end'],
                    "score": ".",
                    "strand": row['strand'],
                    "phase": ".",
                    "attributes": ";".join([x for x in [
                        f"ID={gene_id}_tRNA",
                        f"name={gene_id} tRNA",
                        f"parent={gene_id}",
                        f"product=tRNA-{row['aminoacid']}",
                        anticodon_attr
                    ] if x])
                })
                gff_rows.append({
                    "seqid": row['seqid'],
                    "source": "ARAGORN",
                    "type": "exon",
                    "start": row['intron_end'] + 1,
                    "end": row['end'],
                    "score": ".",
                    "strand": row['strand'],
                    "phase": ".",
                    "attributes": f"ID={gene_id}_exon2;name={gene_id} exon2;parent={gene_id}_tRNA"
                })
            else:
                # Standard tRNA
                gff_rows.append({
                    "seqid": row['seqid'],
                    "source": "ARAGORN",
                    "type": "tRNA",
                    "start": row['start'],
                    "end": row['end'],
                    "score": ".",
                    "strand": row['strand'],
                    "phase": ".",
                    "attributes": ";".join([
                        f"ID={gene_id}_tRNA",
                        f"name={gene_id} tRNA",
                        f"parent={gene_id}",
                        f"product=tRNA-{row['aminoacid']}",
                        anticodon_attr
                    ])
                })

        # tmRNA Logic
        elif row['type'] == "tmRNA":
            gff_rows.append({
                "seqid": row['seqid'],
                "source": "ARAGORN",
                "type": "tmRNA",
                "start": row['start'],
                "end": row['end'],
                "score": ".",
                "strand": row['strand'],
                "phase": ".",
                "attributes": f"ID={tmRNA_id};name=ssrA tmRNA;parent={gene_id}"
            })
            if organelle.lower() == "chloroplast":
                transl_table = 11
            elif organelle.lower() in ["mitochondrion", "mitochondria"]:
                transl_table = 4
            else:
                transl_table = None
            
            gff_rows.append({
                "seqid": row['seqid'],
                "source": "ARAGORN",
                "type": "CDS",
                "start": row['feature_start'],
                "end": row['feature_end'],
                "score": ".",
                "strand": row['strand'],
                "phase": "0",
                "attributes": f"ID={tmRNA_id}_CDS;name={tmRNA_id} CDS;parent={tmRNA_id};transl_table={transl_table}"
            })

    gff_df = pd.DataFrame(gff_rows, columns=["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])

    gff_df['start'] = gff_df['start'].astype(int)
    gff_df['end'] = gff_df['end'].astype(int)

    if export:
        file_name = f"{directory}{seq_name}_AR.gff3"
        with open(file_name, "w") as f:
            f.write("##gff-version 3\n")
        gff_df.to_csv(file_name, sep="\t", header=False, index=False, mode="a")

    return gff_df

def sort_gff3(gff3_object):
    """
    Sorts a GFF3 DataFrame based on a custom biological order and coordinates.

    Parameters:
    -----------
    gff3_object : pd.DataFrame
        The GFF3 data.

    Returns:
    --------
    pd.DataFrame
        Sorted DataFrame.
    """
    custom_order = ["gene", "ncRNA", "rRNA", "tmRNA", "CDS", "tRNA", "intron", "exon"]
    category_type = pd.CategoricalDtype(categories=custom_order, ordered=True)
    gff3_object = gff3_object.copy()
    gff3_object['type'] = gff3_object['type'].astype(category_type)
    # Keep contigs grouped, in their order of first appearance (FASTA order)
    seqid_order = list(dict.fromkeys(gff3_object['seqid']))
    gff3_object['seqid'] = pd.Categorical(gff3_object['seqid'],
                                          categories=seqid_order, ordered=True)
    gff3_object_sorted = gff3_object.sort_values(
        by = ['seqid','start','type'],
        ascending=[True,True,True]
    )
    return gff3_object_sorted

def merge_annotations(mfannot_gff3, aragorn_gff3, seq_name="Default_name", export=True, directory: str = ""):
    """
    Merges MFannot and Aragorn GFF3 DataFrames.
    Aragorn tRNA predictions supersede MFannot tRNA predictions.

    Parameters:
    -----------
    mfannot_gff3 : pd.DataFrame
        MFannot GFF3 data.
    aragorn_gff3 : pd.DataFrame
        Aragorn GFF3 data.
    seq_name : str
        Sequence identifier.
    export : bool
        If True, writes the merged result to a file.
    directory : str or Path
        Output directory.

    Returns:
    --------
    pd.DataFrame
        Merged and sorted GFF3 DataFrame.
    """
    # Remove tRNAs from the MFannot object (Aragorn is preferred for tRNAs)
    # 1. Filter MFannot tRNAs safely
    if not mfannot_gff3.empty:
        filtered_mf = mfannot_gff3.loc[~mfannot_gff3['attributes'].str.contains("ID=trn", na=False)]
    else:
        filtered_mf = mfannot_gff3

    # 2. Avoid a pandas FutureWarning by choosing the non-empty frame
    if filtered_mf.empty and aragorn_gff3.empty:
        merged_df = filtered_mf.copy() # Keeps the column structure
    elif filtered_mf.empty:
        merged_df = aragorn_gff3.copy()
    elif aragorn_gff3.empty:
        merged_df = filtered_mf.copy()
    else:
        merged_df = pd.concat([filtered_mf, aragorn_gff3], ignore_index=True)

    # 3. Only process and sort if there's data. Each feature keeps its own
    #    seqid, so multi-contig genomes stay separated by contig.
    if not merged_df.empty:
        sorted_df = sort_gff3(merged_df)
    else:
        sorted_df = merged_df

    # 4. Standard GFF3 export
    if export:
        file_name = Path(directory) / f"{seq_name}.gff3"
        with open(file_name, "w") as f:
            f.write("##gff-version 3\n")
        # Only append data if sorted_df is not empty
        if not sorted_df.empty:
            sorted_df.to_csv(file_name, sep="\t", header=False, index=False, mode="a")

    return sorted_df

def annotate_fasta(fasta_file, format, directory: str = "", organelle="", circular=False, keep_old=True, keep_intermediates=False):
    """
    Main pipeline: Annotates a FASTA file using MFannot and Aragorn.
    Handles sequence rotation if the genome is circular.

    Parameters:
    -----------
    fasta_file : str
        Input FASTA filename.
    directory : str or Path
        Working directory containing the file.
    organelle : str
        Organelle type ('chloroplast' or 'mitochondrion').
    circular : bool
        If True, attempts to rotate the genome to start at rns.
    keep_old : bool
        If True, keeps the original sequence file before rotation.
    keep_intermediates : bool
        If True, also writes the MFannot table and the Aragorn output next to
        the results instead of discarding them with the temporary directory.
    """
    import time
    start_time = time.time()
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_path = Path(tmpdir)
        
        # Copy fasta to temp directory to work safely
        input_fasta = Path(directory) / fasta_file
        working_fasta = tmp_path / input_fasta.name
        shutil.copy(input_fasta, working_fasta)
        
        fasta_stem = input_fasta.stem
        tbl_file = tmp_path / f"{fasta_stem}.tbl"

        log_path = Path(directory) / f"{fasta_stem}.log"

        def log(text, section=None):
            """Append to the run log, printing anything that is not tool output."""
            with open(log_path, "a") as fh:
                if section:
                    fh.write(f"\n----- {section} -----\n{text.rstrip()}\n")
                else:
                    fh.write(f"{text}\n")
            if not section:
                print(text)

        log_path.write_text(
            f"SOGA annotation of {fasta_file}\n"
            f"  started:    {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            f"  organelle:  {organelle}\n"
            f"  format:     {format}\n"
            f"  circular:   {circular}\n"
            f"  keep_old:   {keep_old}\n"
        )

        def run_mfannot(path, section):
            try:
                out = annotate_mfannot(str(path), tmp_path, organelle)
            except RuntimeError as e:
                log(str(e), section)
                raise
            log(out, section)

        # Real sequence IDs in FASTA order, used to give each feature the right
        # seqid. reorder_fasta treats the file as a single sequence, so rotation
        # is only attempted for single-contig genomes; annotation itself is
        # contig-aware either way.
        contig_ids = fasta_contig_ids(working_fasta)
        # A fragmented assembly is a set of linear pieces, so --circular only
        # takes effect for a single-sequence file.
        treat_circular = circular and len(contig_ids) == 1
        reordered = False

        if circular and not treat_circular:
            log(f"{fasta_stem}: {len(contig_ids)} contigs, annotating each "
                f"without rotation to rns")
        elif circular:
            # Locate rns with nhmmer rather than a full MFannot pass. MFannot
            # loses any gene that spans the origin, so rns is invisible to it
            # in exactly the case where rotation matters most.
            hit = find_rns_start(working_fasta)
            if hit is None:
                log(f"No rns found in {fasta_stem}, proceeding without reordering")
            else:
                strand, start = hit
                if strand == "+" and start == 1:
                    log(f"rns already at position 1 in {fasta_stem}, no rotation needed")
                else:
                    cut_pos = start - 1 if strand == "+" else start
                    reorder_fasta(str(working_fasta), strand=strand, cut_pos=cut_pos,
                                  new_seq_name=str(working_fasta))
                    reordered = True
                    log(f"Rotated {fasta_stem} to start at rns "
                        f"(strand {strand}, position {start})")

        # Generate final annotations
        run_mfannot(working_fasta, "MFannot")
        log(annotate_aragorn(str(working_fasta), tmp_path, organelle, circular), "Aragorn")
        
        mf_df = mfannot_to_gff3(str(tbl_file), fasta_stem, False, tmp_path, organelle,
                                contig_ids=contig_ids)
        ar_df = aragorn_to_gff3(f"{fasta_stem}.txt", fasta_stem, False, tmp_path, organelle)
        if mf_df.empty and ar_df.empty:
            log(f"File {fasta_stem} has no features")
            return
        else:
            merge_annotations(mf_df, ar_df, fasta_stem, True, tmp_path)
        
        #Export in adequate format
        export_gff3 = tmp_path / f"{fasta_stem}.gff3"        
        
        if format == "gb":
            # One GenBank record per contig, keyed by the seqid used in the GFF3
            records = {}
            for rec in SeqIO.parse(str(working_fasta), "fasta"):
                rec.name = rec.id[:16]  # GenBank LOCUS name has a length cap
                rec.annotations["molecule_type"] = "DNA"
                # Without this the LOCUS line has no topology and readers
                # default to linear.
                rec.annotations["topology"] = "circular" if treat_circular else "linear"
                records[rec.id] = rec

            with open(export_gff3) as gff_handle:
                annotated = list(GFF.parse(gff_handle, base_dict=records))

            for rec in annotated:
                flat = []
                stack = list(rec.features)
                while stack:
                    f = stack.pop(0)
                    flat.append(f)
                    if hasattr(f, "sub_features") and f.sub_features:
                        stack = list(f.sub_features) + stack
                records[rec.id].features = flat

            SeqIO.write([records[i] for i in contig_ids if i in records],
                        str(Path(directory) / f"{fasta_stem}.gb"), "gb")

        elif format == "gff3":
            gff3_dst = Path(directory) / f"{fasta_stem}.gff3"
            shutil.move(export_gff3, gff3_dst)

        # If the genome was rotated, the annotation coordinates refer to the
        # rotated sequence, so it has to replace the input file. Whether the
        # pre-rotation sequence is kept as a backup is a separate choice.
        if reordered:
            if keep_old:
                input_fasta.rename(Path(directory) / f"{fasta_stem}_old.fasta")
            shutil.move(working_fasta, Path(directory) / fasta_file)

        if keep_intermediates:
            for src in (tmp_path / f"{fasta_stem}.tbl", tmp_path / f"{fasta_stem}.txt"):
                if src.exists():
                    shutil.copy(src, Path(directory) / src.name)

        elapsed = time.time() - start_time
        m, s = divmod(elapsed, 60)
        log(f"Annotation completed in {int(m)} min {s:.1f} s")
