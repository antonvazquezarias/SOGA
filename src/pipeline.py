import pandas as pd
from collections import defaultdict
import re
import subprocess
from pathlib import Path
import tempfile
import traceback
import time
import shutil
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
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

def reorder_sequence(seq, strand, cut_pos):
    """
    Reorders a circular sequence starting from a specific cut position
    and optionally adjusts the strand.

    Parameters:
    -----------
    seq : str
        The input nucleotide sequence string.
    strand : str
        Target strand ('+' or '-'). If '-', the sequence is reverse complemented.
    cut_pos : int
        The 0-based index at which to cut the sequence. The new sequence will start
        from this position.

    Returns:
    --------
    str
        The rotated (and optionally reverse complemented) sequence string.
    """
    seq = seq[cut_pos:] + seq[:cut_pos]
    if strand == "-":
        complement = str.maketrans("ACGTacgt", "TGCAtgca")
        seq = seq.translate(complement)[::-1]
    return seq

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
    proc = subprocess.run(
        ["mfannot", "-g", str(genetic_code), "--tbl", f"{file_name}.fasta"],
        cwd=folder,                      # <-- Replaces "cd {folder} &&"
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True                        # Modern replacement for universal_newlines=True
    )

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
    ], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

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

    df = sort_gff3(df, contig_ids=contig_ids)

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

        gene_id = gene_base
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

def sort_gff3(gff3_object, contig_ids=None):
    """
    Sorts a GFF3 DataFrame based on a custom biological order and coordinates.

    Parameters:
    -----------
    gff3_object : pd.DataFrame
        The GFF3 data.
    contig_ids : list of str, optional
        Order of sequence IDs to respect when sorting.

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
    if contig_ids is not None:
        seqid_order = [c for c in contig_ids if c in gff3_object['seqid'].values]
        for s in dict.fromkeys(gff3_object['seqid']):
            if s not in seqid_order:
                seqid_order.append(s)
    else:
        seqid_order = list(dict.fromkeys(gff3_object['seqid']))

    gff3_object['seqid'] = pd.Categorical(gff3_object['seqid'],
                                          categories=seqid_order, ordered=True)
    gff3_object_sorted = gff3_object.sort_values(
        by = ['seqid','start','type'],
        ascending=[True,True,True]
    )
    return gff3_object_sorted

def merge_annotations(mfannot_gff3, aragorn_gff3, seq_name="Default_name", export=True, directory: str = "", contig_ids=None):
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
    contig_ids : list of str, optional
        Preserved contig order for sorting.

    Returns:
    --------
    pd.DataFrame
        Merged and sorted GFF3 DataFrame.
    """
    # Remove tRNAs from the MFannot object (Aragorn is preferred for tRNAs)
    # 1. Filter MFannot tRNAs safely
    if not mfannot_gff3.empty:
        filtered_mf = mfannot_gff3.loc[~mfannot_gff3['attributes'].str.contains(r"(?:ID|parent)=trn", na=False)]
    else:
        filtered_mf = mfannot_gff3

    # 2. Avoid a pandas FutureWarning by choosing the non-empty frame
    if filtered_mf.empty and aragorn_gff3.empty:
        merged_df = filtered_mf.copy() # Keeps the column structure
    elif filtered_mf.empty:
        merged_df = aragorn_gff3.copy().reset_index(drop=True)
    elif aragorn_gff3.empty:
        merged_df = filtered_mf.copy().reset_index(drop=True)
    else:
        merged_df = pd.concat([filtered_mf, aragorn_gff3], ignore_index=True)

    # 3. Disambiguate duplicate IDs across sequences and within sequences
    if not merged_df.empty:
        gene_indices = merged_df.index[merged_df['type'] == 'gene'].tolist()
        gene_occurrences = defaultdict(lambda: defaultdict(list))

        # Pass 1: Map where each gene appears across all sequences
        for idx in gene_indices:
            attr = merged_df.at[idx, 'attributes']
            match = re.search(r"(?:^|;)ID=([^;]+)", attr)
            base_id = match.group(1) if match else "unknown_gene"
            seqid = merged_df.at[idx, 'seqid']
            gene_occurrences[base_id][seqid].append(idx)

        # Pass 2: Compute disambiguated IDs and cascade strictly to ID and parent
        for base_id, seq_dict in gene_occurrences.items():
            multi_seq = len(seq_dict) > 1

            for seqid, idx_list in seq_dict.items():
                multi_in_seq = len(idx_list) > 1

                for copy_num, g_idx in enumerate(idx_list, start=1):
                    # Add {seqid}_ prefix if present in multiple sequences; add _1, _2 suffix if multiple on this sequence
                    prefix = f"{seqid}_" if multi_seq else ""
                    suffix = f"_{copy_num}" if multi_in_seq else ""
                    new_gene_id = f"{prefix}{base_id}{suffix}"

                    if new_gene_id == base_id:
                        continue  # Globally unique, leave unmodified

                    g_start = merged_df.at[g_idx, 'start']
                    g_end = merged_df.at[g_idx, 'end']

                    # Target only ID= and parent=/Parent= fields, leaving name= and gene= untouched
                    pattern = re.compile(rf"((?:^|;)(?:ID|parent|Parent)=){re.escape(base_id)}(?=[;_]|$)")

                    # Update the parent gene row ID
                    old_attr = merged_df.at[g_idx, 'attributes']
                    merged_df.at[g_idx, 'attributes'] = pattern.sub(rf"\g<1>{new_gene_id}", old_attr)

                    # Update child feature ID and parent references within the gene span
                    candidate_mask = (
                        (merged_df['seqid'] == seqid) &
                        (merged_df['start'] >= g_start) &
                        (merged_df['end'] <= g_end) &
                        (merged_df.index != g_idx) &
                        (merged_df['attributes'].str.contains(rf"(?:ID|parent|Parent)={re.escape(base_id)}", regex=True, na=False))
                    )

                    for child_idx in merged_df.index[candidate_mask]:
                        c_attr = merged_df.at[child_idx, 'attributes']
                        merged_df.at[child_idx, 'attributes'] = pattern.sub(rf"\g<1>{new_gene_id}", c_attr)

    # 4. Only process and sort if there's data. Each feature keeps its own
    #    seqid, so multi-contig genomes stay separated by contig.
    if not merged_df.empty:
        sorted_df = sort_gff3(merged_df, contig_ids=contig_ids)
    else:
        sorted_df = merged_df

    # 5. Standard GFF3 export
    if export:
        file_name = Path(directory) / f"{seq_name}.gff3"
        with open(file_name, "w") as f:
            f.write("##gff-version 3\n")
        # Only append data if sorted_df is not empty
        if not sorted_df.empty:
            sorted_df.to_csv(file_name, sep="\t", header=False, index=False, mode="a")

    return sorted_df

def annotate_sequence(
    seq_file,
    format,
    directory: str = "",
    organelle="",
    circular=None,
    keep_old=True,
    keep_intermediates=False,
    output_dir: str = None,
    save_log: bool = False,
    return_log: bool = False,
):
    """
    Main pipeline: Annotates a FASTA or GenBank sequence file using MFannot and Aragorn.
    Handles sequence rotation if circular sequences are present.
    """
    start_time = time.time()

    # ---------------------------------------------------------
    # 1. PATH RESOLUTION & VALIDATION
    # ---------------------------------------------------------
    in_dir = Path(directory)
    out_dir = Path(output_dir) if output_dir else in_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    is_in_place = in_dir.resolve() == out_dir.resolve()

    input_file_path = in_dir / seq_file
    input_suffix = input_file_path.suffix.lower()

    valid_fasta = {".fasta", ".fa", ".fna"}
    valid_gb = {".gb", ".gbk", ".genbank"}
    if input_suffix not in (valid_fasta | valid_gb):
        raise ValueError(
            f"Unsupported file format '{input_suffix}' for {seq_file}. "
            f"Supported formats are FASTA ({', '.join(valid_fasta)}) and GenBank ({', '.join(valid_gb)})."
        )
    is_gb_input = input_suffix in valid_gb

    seq_stem = input_file_path.stem.replace(" ", "_")
    log_path = out_dir / f"{seq_stem}.log"

    # In-memory log buffer
    log_lines = []

    def log(text, section=None):
        """Buffer log messages and print non-section messages to terminal immediately."""
        formatted = f"\n----- {section} -----\n{text.rstrip()}\n" if section else f"{text}"
        log_lines.append(formatted)
        if not section:
            print(text, flush=True)

    # Clean header: Detailed in single mode, compact in batch mode
    if not return_log:
        log(
            f"SOGA annotation of {seq_file}\n"
            f"  started:    {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            f"  organelle:  {organelle}\n"
            f"  format:     {format}\n"
            f"  circular:   {circular}\n"
            f"  keep_old:   {keep_old}\n"
            f"  output_dir: {out_dir}"
        )
    else:
        log(f"Started: {time.strftime('%H:%M:%S')} | Organelle: {organelle} | Circular: {circular}")

    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            working_fasta = tmp_path / f"{seq_stem}.fasta"

            # ---------------------------------------------------------
            # 2. INGEST INPUT & INITIALIZE RECORDS
            # ---------------------------------------------------------
            gb_records = {}
            if is_gb_input:
                for rec in SeqIO.parse(str(input_file_path), "genbank"):
                    rec.features = []  # Clear previous features, keep metadata
                    gb_records[rec.id] = rec
                SeqIO.write(list(gb_records.values()), str(working_fasta), "fasta")
            else:
                shutil.copy(input_file_path, working_fasta)

            records_dict = {rec.id: rec for rec in SeqIO.parse(str(working_fasta), "fasta")}
            contig_ids = list(records_dict.keys())
            num_seqs = len(contig_ids)
            total_bp = sum(len(r.seq) for r in records_dict.values())

            # Report contigs and total genome size
            log(f"Sequences: {num_seqs} contig(s) ({total_bp:,} bp total)")

            # ---------------------------------------------------------
            # 3. TOPOLOGY DETERMINATION
            # ---------------------------------------------------------
            is_circular = {}
            for cid in contig_ids:
                if is_gb_input:
                    rec_topo = gb_records[cid].annotations.get("topology", "linear").lower()
                    rec_is_circular = (rec_topo == "circular")
                    if circular is not None:
                        if rec_is_circular != circular:
                            msg = (
                                f"Contig '{cid}' has contrasting topology '{rec_topo}' in {seq_file}, "
                                f"using circular={circular} as specified."
                            )
                            warnings.warn(msg)
                            log(f"WARNING: {msg}")
                        is_circular[cid] = circular
                    else:
                        is_circular[cid] = rec_is_circular
                else:
                    is_circular[cid] = bool(circular) if circular is not None else False

            # ---------------------------------------------------------
            # 4. SEQUENCE ROTATION (TO RNS START)
            # ---------------------------------------------------------
            reordered = False
            for cid in contig_ids:
                if not is_circular[cid]:
                    continue

                single_tmp = tmp_path / "single_rns_check.fasta"
                SeqIO.write(records_dict[cid], str(single_tmp), "fasta")
                hit = find_rns_start(single_tmp.name, directory=tmp_path)

                if hit is None:
                    log(f"No rns found in contig {cid}, proceeding without reordering")
                else:
                    strand, start = hit
                    if strand == "+" and start == 1:
                        log(f"rns already at position 1 in contig {cid}, no rotation needed")
                    else:
                        cut_pos = start - 1 if strand == "+" else start
                        new_seq_str = reorder_sequence(str(records_dict[cid].seq), strand=strand, cut_pos=cut_pos)

                        records_dict[cid].seq = Seq(new_seq_str)
                        if is_gb_input and cid in gb_records:
                            gb_records[cid].seq = Seq(new_seq_str)

                        reordered = True
                        log(f"Rotated contig {cid} to start at rns (strand {strand}, position {start})")

            if reordered:
                SeqIO.write(list(records_dict.values()), str(working_fasta), "fasta")

            # ---------------------------------------------------------
            # 5. FEATURE ANNOTATION (MFannot & Aragorn)
            # ---------------------------------------------------------
            try:
                mfannot_out = annotate_mfannot(str(working_fasta), tmp_path, organelle)
                log(mfannot_out, "MFannot")
            except RuntimeError as e:
                log(str(e), "MFannot")
                raise

            has_circular = any(is_circular.values())
            has_linear = any(not c for c in is_circular.values())

            if has_circular and has_linear:
                circ_records = [records_dict[cid] for cid in contig_ids if is_circular[cid]]
                lin_records = [records_dict[cid] for cid in contig_ids if not is_circular[cid]]

                circ_fasta = tmp_path / f"{seq_stem}_circ.fasta"
                SeqIO.write(circ_records, str(circ_fasta), "fasta")
                log(annotate_aragorn(str(circ_fasta), tmp_path, organelle, circular=True), "Aragorn (circular)")
                ar_circ_df = aragorn_to_gff3(f"{seq_stem}_circ.txt", seq_stem, False, tmp_path, organelle)

                lin_fasta = tmp_path / f"{seq_stem}_lin.fasta"
                SeqIO.write(lin_records, str(lin_fasta), "fasta")
                log(annotate_aragorn(str(lin_fasta), tmp_path, organelle, circular=False), "Aragorn (linear)")
                ar_lin_df = aragorn_to_gff3(f"{seq_stem}_lin.txt", seq_stem, False, tmp_path, organelle)

                ar_df = pd.concat([ar_circ_df, ar_lin_df], ignore_index=True)
                ar_df = sort_gff3(ar_df, contig_ids=contig_ids)
            else:
                log(annotate_aragorn(str(working_fasta), tmp_path, organelle, circular=has_circular), "Aragorn")
                ar_df = aragorn_to_gff3(f"{seq_stem}.txt", seq_stem, False, tmp_path, organelle)

            tbl_file = tmp_path / f"{seq_stem}.tbl"
            mf_df = mfannot_to_gff3(str(tbl_file), seq_stem, False, tmp_path, organelle, contig_ids=contig_ids)

            if mf_df.empty and ar_df.empty:
                log(f"File {seq_stem} has no features")
                return "\n".join(log_lines) if return_log else None

            merge_annotations(mf_df, ar_df, seq_stem, True, tmp_path, contig_ids=contig_ids)
            export_gff3 = tmp_path / f"{seq_stem}.gff3"

            # ---------------------------------------------------------
            # 6. BACKUP & EXPORT
            # ---------------------------------------------------------
            if reordered and keep_old and is_in_place and input_file_path.exists():
                input_file_path.rename(in_dir / f"{seq_stem}_old{input_suffix}")

            if reordered and not (format == "gb" and is_gb_input):
                if is_gb_input:
                    SeqIO.write([gb_records[i] for i in contig_ids if i in gb_records], str(out_dir / seq_file), "gb")
                else:
                    shutil.copy(working_fasta, out_dir / seq_file)

            if format == "gb":
                if is_gb_input:
                    records = gb_records
                    for rec_id, rec in records.items():
                        rec.annotations["molecule_type"] = "DNA"
                        rec.annotations["topology"] = "circular" if is_circular.get(rec_id, False) else "linear"
                else:
                    records = {}
                    for rec in SeqIO.parse(str(working_fasta), "fasta"):
                        rec.name = rec.id[:16]
                        rec.annotations["molecule_type"] = "DNA"
                        rec.annotations["topology"] = "circular" if is_circular.get(rec_id, False) else "linear"
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

                SeqIO.write([records[i] for i in contig_ids if i in records], str(out_dir / f"{seq_stem}.gb"), "gb")

            elif format == "gff3":
                shutil.copy(export_gff3, out_dir / f"{seq_stem}.gff3")

            # ---------------------------------------------------------
            # 7. CLEANUP & FINISH
            # ---------------------------------------------------------
            if keep_intermediates:
                for filename in [f"{seq_stem}.tbl", f"{seq_stem}.txt", f"{seq_stem}_circ.txt", f"{seq_stem}_lin.txt"]:
                    src = tmp_path / filename
                    if src.exists():
                        shutil.copy(src, out_dir / src.name)

            elapsed = time.time() - start_time
            m, s = divmod(elapsed, 60)
            log(f"Annotation completed in {int(m)} min {s:.1f} s")

        # Save single-file log if explicitly requested
        if save_log and not return_log:
            log_path.write_text("\n".join(log_lines) + "\n")

        if return_log:
            return "\n".join(log_lines)

    except Exception as e:
        log(f"\n[CRITICAL ERROR] Failed during annotation: {e}\n{traceback.format_exc()}")
        # In single file mode, ALWAYS dump the log on error
        if not return_log:
            log_path.write_text("\n".join(log_lines) + "\n")
            print(f"\n[ERROR] An error occurred. Details saved to log: {log_path}", flush=True)
        raise