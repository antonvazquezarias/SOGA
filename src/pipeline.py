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
        
def reorder_fasta_by_rns(fasta_file, tbl_file, directory: str = "", keep_old=True):
    """
    Rotates a circular genome so that it starts with the small subunit ribosomal RNA (rns).

    Parameters:
    -----------
    fasta_file : str
        Input FASTA file.
    tbl_file : str
        MFannot output table file containing feature coordinates.
    directory : str or Path
        Working directory.
    keep_old : bool
        If True, keeps the original non-rotated fasta as a backup.

    Raises:
    -------
    ValueError
        If 'rns' is not found or is already at position 1.
    """
    seq_name = Path(fasta_file).stem
    fasta_path = Path(directory, fasta_file)

    # Load annotation to locate features
    annotation = mfannot_to_gff3(tbl_file, seq_name, export=False, directory=directory, organelle="")

    mask = annotation['attributes'].str.contains("ID=rns", na=False)
    if not mask.any():
        raise ValueError("No rns found")

    if (annotation.loc[mask, "start"] == 1).any():
        raise ValueError("rns already at position 1")

    rns_row = annotation.loc[mask].iloc[0]

    with open(fasta_path) as f:
        seq = "".join(line.strip() for line in f)
        genome_length = len(seq)
    
    rns_length = rns_row['end'] - rns_row['start']
    strand = rns_row['strand']

    # Check if rns crosses the origin (will lead to large length)
    if rns_length/genome_length > 0.5:
        rns_split = True
    else:
        rns_split = False

    # Determine cut position based on strand and whether gene is split
    if not rns_split:
        cut_pos = rns_row['start'] - 1 if strand == "+" else rns_row['end']
    else:
        strand = "-" if strand == "+" else "+"
        cut_pos = rns_row['end'] - 1 if strand == "+" else rns_row['start']

    old_seq_name = f"{seq_name}_old.fasta" if keep_old else None

    reorder_fasta(fasta_file, strand, cut_pos, new_seq_name=fasta_file, old_seq_name=old_seq_name, directory=directory)

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
    
    # Run MFannot via Docker, moving the .tbl to the mounted folder
    subprocess.run([
        "sh", "-c",
        f"cd {folder} && /mfannot/mfannot -g {genetic_code} --tbl {file_name}.fasta > /dev/null 2>&1"
    ], check=True)

    candidates = sorted(folder.glob(f"{file_name}.fasta.new*.tbl"))
    if not candidates:
        raise FileNotFoundError(
            f"MFannot did not produce a .tbl file for {file_name} in {folder}"
        )
    tbl_file = candidates[-1]
    shutil.move(tbl_file, folder / f"{file_name}.tbl")

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

    # Run Aragorn in Docker in batch mode
    subprocess.run([
        "aragorn",
        "-t",       # tRNA
        "-m",       # tmRNA
        "-i",       # introns
        f"-{topology}",
        f"-g{genetic_code}",
        "-w",       # batch mode
        "-o", f"{folder}/{file_name}.txt",
        f"{folder}/{file_name}.fasta"
    ], check=True)

def mfannot_to_gff3(tbl_file, seq_name="Default_name", export=True, directory: str = "", organelle=""):
    """
    Parses MFannot .tbl output and converts it to GFF3 format.

    Parameters:
    -----------
    tbl_file : str
        Name of the MFannot .tbl file.
    seq_name : str
        Sequence identifier to use in the GFF3.
    export : bool
        If True, writes the result to a .gff3 file.
    directory : str or Path
        Working directory.
    organelle : str
        Used to set translation table attributes.

    Returns:
    --------
    pd.DataFrame
        DataFrame containing the GFF3 data.
    """

    gff_columns = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    # Read file, skipping header and footer lines specific to MFannot output
    with open(Path(directory) / tbl_file) as f:
        lines = f.read().splitlines()[1:-1]
    
    if len(lines) == 0:
        empty_df = pd.DataFrame(columns=gff_columns)
        
        if export:
            file_name = Path(directory) / f"{seq_name}_MF.gff3"
            with open(file_name, "w") as f:
                f.write("##gff-version 3\n")
            # This ensures the file exists and is valid GFF3 even if empty
            empty_df.to_csv(file_name, sep="\t", header=False, index=False, mode="a")
            
        return empty_df


    # Parse features and their qualifiers
    features = []
    current_feature = None

    for line in lines:
        elements = line.split("\t")
        if elements[0]:
            current_feature = {"feature": elements, "qualifiers": []}
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
            seq_name,
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
        Sequence identifier for GFF3.
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
    
    file_path = Path(directory) / txt_file
    with open(file_path) as f:
        full_content = f.read()

    # Early exit if no genes are found
    if "0 genes found" in full_content or not full_content.strip():
        empty_df = pd.DataFrame(columns=gff_columns)
        
        if export:
            file_name = Path(directory) / f"{seq_name}_AR.gff3"
            with open(file_name, "w") as f:
                f.write("##gff-version 3\n")
            # This writes the column names even if there is no data
            empty_df.to_csv(file_name, sep="\t", header=False, index=False, mode="a")
            
        return empty_df

    # Process lines normally
    lines = full_content.splitlines()[2:]
    lines = [re.sub(r"[ \t]+", " ", line).strip() for line in lines if line.strip()]
    df = pd.DataFrame([line.split(" ") for line in lines])
    
    # Cleaning and Pre-processing columns
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

    df = df[['type', 'start', 'end', 'strand', 'codon', 'aminoacid', 'feature_start', 'feature_end', 'intron_start', 'intron_end']]

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
            "seqid": seq_name,
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
                    "seqid": seq_name,
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
                    "seqid": seq_name,
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
                    "seqid": seq_name,
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
                    "seqid": seq_name,
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
                    "seqid": seq_name,
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
                    "seqid": seq_name,
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
                "seqid": seq_name,
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
                "seqid": seq_name,
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
    gff3_object['type'] = gff3_object['type'].astype(category_type)
    gff3_object_sorted = gff3_object.sort_values(
        by = ['start','type'],
        ascending=[True,True]
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

    # 2. Avoid the FutureWarning by choosing the non-empty DF
    if filtered_mf.empty and aragorn_gff3.empty:
        merged_df = filtered_mf.copy() # Keeps the column structure
    elif filtered_mf.empty:
        merged_df = aragorn_gff3.copy()
    elif aragorn_gff3.empty:
        merged_df = filtered_mf.copy()
    else:
        merged_df = pd.concat([filtered_mf, aragorn_gff3], ignore_index=True)

    # 3. Only process and sort if there's data[cite: 1]
    if not merged_df.empty:
        merged_df['seqid'] = seq_name
        sorted_df = sort_gff3(merged_df)
    else:
        sorted_df = merged_df

    # 4. Standard GFF3 Export[cite: 1]
    if export:
        file_name = Path(directory) / f"{seq_name}.gff3"
        with open(file_name, "w") as f:
            f.write("##gff-version 3\n")
        # Only append data if sorted_df is not empty[cite: 1]
        if not sorted_df.empty:
            sorted_df.to_csv(file_name, sep="\t", header=False, index=False, mode="a")

    return sorted_df

def annotate_fasta(fasta_file, format, directory: str = "", organelle="", circular=False, keep_old=True):
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

        reordered = False
        tbl_ready = False

        if circular:
            reordered = True
            # Annotate first attempt to find landmarks
            annotate_mfannot(str(working_fasta), tmp_path, organelle)

            try:
                reorder_fasta_by_rns(str(working_fasta), str(tbl_file), tmp_path, keep_old=False)

            except ValueError as e:
                msg = str(e)

                # Case 1: no rns → rotate arbitrarily to avoid split features and try again
                if msg == "No rns found":
                    rotated_fasta = tmp_path / f"{fasta_stem}_rotated.fasta"
                    rotated_tbl = tmp_path / f"{fasta_stem}_rotated.tbl"

                    reorder_fasta(str(working_fasta), strand="+", cut_pos=3000, new_seq_name=str(rotated_fasta))
                    annotate_mfannot(str(rotated_fasta), tmp_path, organelle)

                    try:
                        reorder_fasta_by_rns(str(rotated_fasta), str(rotated_tbl), tmp_path, keep_old=False)
                        working_fasta.unlink()
                        rotated_fasta.rename(working_fasta)
                        tbl_file.unlink()
                        rotated_tbl.rename(tbl_file)
                    except ValueError:
                        print(f"No rns found in {fasta_stem}, proceeding without reordering")
                        reordered = False

                # Case 2: rns already at start → use existing annotation
                elif msg == "rns already at position 1":
                    print(f"rns already at position 1 in {fasta_stem}, using existing annotation")
                    reordered = False
                    tbl_ready = True

                else:
                    raise

        # Generate final annotations
        if not tbl_ready:
            annotate_mfannot(str(working_fasta), tmp_path, organelle)
        annotate_aragorn(str(working_fasta), tmp_path, organelle, circular)
        
        mf_df = mfannot_to_gff3(str(tbl_file), fasta_stem, False, tmp_path, organelle)
        ar_df = aragorn_to_gff3(f"{fasta_stem}.txt", fasta_stem, False, tmp_path, organelle)
        if mf_df.empty and ar_df.empty:
            print(f"File {fasta_stem} has no features")
            return
        else:
            merge_annotations(mf_df, ar_df, fasta_stem, True, tmp_path)
        
        #Export in adequate format
        export_gff3 = tmp_path / f"{fasta_stem}.gff3"        
        
        if format == "gb":
            bio_record = SeqIO.read(working_fasta, "fasta")
            bio_record.id = fasta_stem
            bio_record.name = fasta_stem[:16]  # GenBank LOCUS name has a length cap
            with open(export_gff3) as gff_handle:
                result = next(GFF.parse(gff_handle, base_dict={bio_record.id: bio_record}))
            flat = []
            stack = list(result.features)
            while stack:
                f = stack.pop(0)
                flat.append(f)
                if hasattr(f, "sub_features") and f.sub_features:
                    stack = list(f.sub_features) + stack

            bio_record.features = flat
            bio_record.annotations["molecule_type"] = "DNA"
            SeqIO.write(bio_record, Path(directory) / f"{fasta_stem}.gb", "gb")

        elif format == "gff3":
            gff3_dst = Path(directory) / f"{fasta_stem}.gff3"
            shutil.move(export_gff3, gff3_dst)
            if keep_old and reordered:
                old_path = Path(directory) / f"{fasta_stem}_old.fasta"
                input_fasta.rename(old_path)
                fasta_dst = Path(directory) / fasta_file
                shutil.move(working_fasta, fasta_dst)

    
    elapsed = time.time() - start_time
    m, s = divmod(elapsed, 60)
    print(f"Annotation completed in {int(m)} min {s:.1f} s")
