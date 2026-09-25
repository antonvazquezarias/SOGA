from collections import defaultdict, Counter
from dataclasses import dataclass, field
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
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

# Functions
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


def find_rns_start(genes, seqid):
    """
    Finds the 5' end of rns on a contig, from the genes found by nhmmer.

    A fragmented rns starts at its first module (rns_a). If rns is present
    more than once (e.g. in an inverted repeat), the highest-scoring copy
    is used.

    Parameters:
    -----------
    genes : list of Gene
        Output of parse_hmmer.
    seqid : str
        Contig to search.

    Returns:
    --------
    tuple or None
        (strand, position) where position is the 1-based coordinate of the
        first base of rns on that strand, or None if the contig has no rns.
    """
    candidates = [g for g in genes if g.seqid == seqid and g.name in ("rns", "rns_a")]
    if not candidates:
        return None
    rns = max(candidates, key=lambda g: g.children[0].score)
    start, end = rns.parts[0]
    return rns.strand, start if rns.strand == "+" else end


def annotate_hmmer(fasta_file, directory, genetic_code):
    """
    Runs nhmmer with the custom rnl/rns models to locate the rRNA genes.

    Writes <name>_nhmmer.tbl (hit table) and <name>_nhmmer.out (alignments,
    one line each), both read by parse_hmmer.

    Parameters:
    -----------
    fasta_file : str
        Input FASTA file.
    directory : str or Path
        Directory containing the file.
    genetic_code : int
        NCBI translation table. 11 uses the plastid models, anything else
        the mitochondrial ones.
    """
    folder = Path(directory).resolve()
    file_name = Path(fasta_file).stem
    organelle = "chloro" if genetic_code == 11 else "mito"

    proc = subprocess.run([
        "nhmmer",
        "--tblout", f"{file_name}_nhmmer.tbl",
        "-o", f"{file_name}_nhmmer.out",
        "--notextw",    # one line per alignment
        f"/SOGA_data/hmms/rrna_{organelle}.hmm",
        f"{file_name}.fasta"
    ], cwd=folder, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    if proc.returncode != 0:
        log = (proc.stdout or "").strip()
        tail = "\n".join(log.splitlines()[-20:]) or "(nhmmer produced no output)"
        raise RuntimeError(
            f"nhmmer failed for {file_name} in {folder}\n"
            f"  exit code: {proc.returncode}\n"
            f"  last lines of nhmmer output:\n{tail}"
        )
    return proc.stdout or ""


def annotate_mfannot(fasta_file, directory, genetic_code):
    """
    Runs MFannot using a Docker container to annotate the sequence.

    Parameters:
    -----------
    fasta_file : str
        Input FASTA file.
    directory : str or Path
        Directory containing the file (will be mounted to Docker).
    genetic_code : int
        NCBI translation table (e.g. 11 for plastids, 4 for mitochondria).
    """
    folder = Path(directory).resolve()
    file_name = Path(fasta_file).stem

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

def annotate_aragorn(fasta_file, directory, genetic_code, circular=True):
    """
    Runs Aragorn using a Docker container to identify tRNAs and tmRNAs.

    Parameters:
    -----------
    fasta_file : str
        Input FASTA file.
    directory : str or Path
        Directory containing the file.
    genetic_code : int
        NCBI translation table (e.g. 11 for plastids, 4 for mitochondria).
    circular : bool
        Topology of the sequence. True for circular, False for linear.
    """
    folder = Path(directory).resolve()
    file_name = Path(fasta_file).stem

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

# ============================================================
# Annotation model
# ============================================================
# The parsers turn MFannot, Aragorn and nhmmer output into a list of Gene objects,
# merging works on that list, and GFF3 or GenBank output is written from it
# at the very end. A child is linked to its parent by being stored inside it,
# not by an ID string, so two genes with the same name can never be confused.
# GFF3 IDs are only created when writing the file.

@dataclass(eq=False)
class Feature:
    """
    Something annotated inside a gene: CDS, tRNA, rRNA, tmRNA, exon, intron...

    parts holds one (start, end) pair per segment, 1-based and inclusive,
    with start <= end. Segments are listed 5' to 3', so in descending order
    on the minus strand. A tRNA interrupted by an intron is ONE Feature with
    two parts, not two features.
    """
    type: str
    parts: list
    strand: str
    attrs: dict = field(default_factory=dict)
    children: list = field(default_factory=list)
    score: float | None = None


@dataclass(eq=False)
class Gene:
    """A gene and everything annotated inside it (in children)."""
    seqid: str
    source: str
    name: str
    parts: list
    strand: str
    attrs: dict = field(default_factory=dict)
    children: list = field(default_factory=list)


def get_gene_features(gene):
    """
    Returns every feature inside a gene (its children, their children, and
    so on) as a flat list, each parent before its children.
    """
    features = []
    pending = list(gene.children)
    while pending:
        feat = pending.pop(0)
        features.append(feat)
        pending = feat.children + pending
    return features


# Order of feature types that start at the same position, used when sorting output
TYPE_ORDER = {"gene": 0, "ncRNA": 1, "rRNA": 2, "tmRNA": 3, "CDS": 4, "tRNA": 5, "intron": 6, "exon": 7}


# ============================================================
# Parsers
# ============================================================

def parse_mfannot(tbl_file, genetic_code=None, contig_ids=None, seq_name="Default_name"):
    """
    Parses MFannot .tbl output into a list of Gene objects.

    Parameters:
    -----------
    tbl_file : str or Path
        Path to the MFannot .tbl file.
    genetic_code : int, optional
        Written as transl_table on CDS features.
    contig_ids : list of str, optional
        Real sequence IDs in FASTA order. MFannot renames contigs to C_0, C_1,
        ... in that order, so the Nth ">Feature" section is given
        contig_ids[N]. If None, seq_name is used for every feature.
    seq_name : str
        Fallback sequence identifier.

    Returns:
    --------
    list of Gene
    """
    def qualifiers_to_attrs(qual):
        """Converts MFannot qualifier lines into an attribute dict."""
        # MFannot writes transl_except as "(pos:21144..21146, aa : Q)"; the
        # INSDC format is "(pos:21144..21146,aa:Gln)"
        three_letter = {
            "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
            "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
            "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
            "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
            "U": "Sec", "O": "Pyl", "*": "TERM"
        }
        attrs = {}
        for q in qual:
            if len(q) < 5:
                continue
            key, val = q[3], q[4]
            # Skip gene names containing "orf"
            if key == "gene" and "orf" in val:
                continue
            if key == "protein_id":
                val = val.replace("lcl| ", "")
            if key == "transl_except":
                match = re.match(r"\(pos:\s*(\S+?)\s*,\s*aa\s*:\s*(\S+?)\s*\)", val)
                if match:
                    aa = three_letter.get(match.group(2), match.group(2))
                    val = f"(pos:{match.group(1)},aa:{aa})"
            attrs[key] = val
        return attrs

    with open(tbl_file) as f:
        raw_lines = f.read().splitlines()

    # Pass 1: group each feature line with the qualifier lines below it, and
    # record which contig section (">Feature C_<n> ...") it belongs to.
    entries = []
    current = None
    contig_index = -1
    for line in raw_lines:
        if not line.strip():
            continue
        if line.startswith(">"):
            # Trust the number MFannot writes rather than counting sections,
            # so a missing section cannot shift the map.
            match = re.match(r">Feature\s+C_(\d+)", line)
            contig_index = int(match.group(1)) if match else contig_index + 1
            current = None
            continue
        elements = line.split("\t")
        if elements[0]:
            current = {"feature": elements, "qualifiers": [], "contig_index": contig_index}
            entries.append(current)
        elif current is not None:
            current["qualifiers"].append(elements)

    # Pass 2: build genes. In the .tbl format, a feature made of several
    # segments is written as one line with a type followed by lines without
    # one, and its qualifiers come after the last segment.
    genes = []
    gene = None           # gene being filled
    feature = None        # latest CDS/rRNA/tRNA... of that gene; introns and exons go inside it
    last = None           # latest gene or feature created, extended by continuation lines
    previous_type = ""

    for entry in entries:
        f, qual = entry["feature"], entry["qualifiers"]
        ci = entry["contig_index"]
        seqid = contig_ids[ci] if (contig_ids and 0 <= ci < len(contig_ids)) else seq_name

        # A new contig section never continues the previous gene
        if gene is not None and gene.seqid != seqid:
            gene = feature = last = None

        # Feature type. Some ORFs come without a type and are recognized by
        # their protein_id (which MFannot builds from the gene name, e.g.
        # "lcl| G-ycf3"); any other line without a type continues the
        # previous feature.
        has_type = len(f) > 2 and bool(f[2])
        protein_id = next((q[4] for q in qual if len(q) >= 5 and q[3] == "protein_id"), "")
        if has_type:
            row_type = f[2]
        elif "ycf" in protein_id:
            row_type = "CDS"
            qual[0][4] = "hypothetical protein"
        elif "rnz" in protein_id:
            row_type = "CDS"
            qual[0][4] = "Ribonuclease Z"
        elif "odp" in protein_id.lower():
            row_type = "CDS"
            if "odpa" in protein_id.lower():
                qual[0][4] = "Pyruvate dehydrogenase E1 component subunit alpha"
            elif "odpb" in protein_id.lower():
                qual[0][4] = "Pyruvate dehydrogenase E1 component subunit beta"
        else:
            row_type = previous_type
        if row_type == "RNA":
            row_type = "ncRNA"

        is_continuation = not has_type and last is not None and row_type == previous_type
        previous_type = row_type

        # Coordinates: MFannot writes minus-strand features as end..start
        a, b = int(f[0]), int(f[1])
        start, end = min(a, b), max(a, b)
        strand = "+" if a < b else "-"

        attrs = qualifiers_to_attrs(qual)

        if is_continuation:
            last.parts.append((start, end))
            last.attrs.update(attrs)
            if last is gene and gene.name == "gene_unknown" and qual:
                gene.name = qual[0][4]
        elif row_type == "gene":
            name = qual[0][4] if qual else "gene_unknown"
            gene = Gene(seqid, "MFannot", name, [(start, end)], strand, attrs)
            genes.append(gene)
            feature = None
            last = gene
        elif gene is None:
            continue  # MFannot always writes the gene line first
        elif row_type in ("intron", "exon"):
            last = Feature(row_type, [(start, end)], strand, attrs)
            parent = feature if feature is not None else gene
            parent.children.append(last)
        else:
            if row_type == "CDS" and genetic_code:
                attrs = {"transl_table": genetic_code, **attrs}
            feature = last = Feature(row_type, [(start, end)], strand, attrs)
            gene.children.append(feature)

    return genes


def parse_aragorn(txt_file, genetic_code=None, seq_name="Default_name", contig_lengths=None):
    """
    Parses Aragorn batch (-w) output into a list of Gene objects.

    On circular contigs, Aragorn writes a gene crossing the origin as
    [start,end] with start > end. Such genes are built with coordinates
    continuing past the contig end, then split at the origin.

    Parameters:
    -----------
    txt_file : str or Path
        Path to the Aragorn output file.
    genetic_code : int, optional
        Written as transl_table on the tmRNA coding region.
    seq_name : str
        Fallback sequence identifier, used only if a data row precedes any
        contig header, which does not happen in normal Aragorn output.
    contig_lengths : dict, optional
        seqid -> length, needed to split genes that cross the origin.

    Returns:
    --------
    list of Gene
    """
    aminoacids = {
        "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
        "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
        "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
        "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
        "Sec": "U", "Pyl": "O"
    }

    with open(txt_file) as f:
        raw_lines = f.read().splitlines()

    # Batch (-w) output holds one section per contig:
    #   >contig_id description
    #   N genes found
    #   1  tRNA-Xxx  [a,b]  pos  (codon)
    # and closes with a ">end ..." line. Keep the contig id for every data row,
    # since a contig with no genes must not discard the others.
    genes = []
    seqid = seq_name
    for line in raw_lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            if not stripped.startswith(">end"):
                seqid = stripped[1:].split()[0]
            continue
        if stripped.endswith("genes found"):
            continue

        # Columns: row number, type-aminoacid, [start,end] (c[...] on the minus
        # strand), anticodon position (tRNA) or coding region (tmRNA), and
        # codon (tRNA) or peptide tag (tmRNA)
        fields = re.sub(r"[ \t]+", " ", stripped).split(" ")
        coords = re.match(r"c?\[(\d+),(\d+)\]", fields[2]) if len(fields) > 2 else None
        if not coords:
            continue  # not a gene row
        start, end = int(coords.group(1)), int(coords.group(2))
        strand = "-" if fields[2].startswith("c") else "+"

        # Gene crossing the origin: continue the coordinates past the contig end
        contig_length = (contig_lengths or {}).get(seqid)
        crosses_origin = start > end and contig_length is not None
        if crosses_origin:
            end += contig_length
        mol_type, _, aminoacid = fields[1].partition("-")
        position = fields[3] if len(fields) > 3 else ""
        codon_field = fields[4] if len(fields) > 4 else ""

        # A tRNA intron is written after the codon as i(distance,length)
        intron = re.search(r"i\((\d+),(\d+)\)", codon_field)
        codon = re.sub(r"i\(\d+,\d+\)", "", codon_field).strip()

        if mol_type == "tRNA":
            name = f"trn{aminoacids.get(aminoacid, '')}{codon}"
        elif mol_type == "tmRNA":
            name = "ssrA"
        else:
            name = mol_type

        gene = Gene(seqid, "ARAGORN", name, [(start, end)], strand, {"gene": name})
        genes.append(gene)

        if mol_type == "tRNA":
            attrs = {"product": f"tRNA-{aminoacid}"}

            if intron:
                # tRNA split by an intron: one tRNA feature with two parts,
                # one exon per part, and the intron directly under the gene.
                # (No anticodon attribute for split tRNAs.)
                distance, length = int(intron.group(1)), int(intron.group(2))
                if strand == "+":
                    intron_start = start + distance - 1
                    intron_end = start + distance + length - 2
                else:
                    intron_end = end - distance + 1
                    intron_start = end - distance - length + 2
                exons = [(start, intron_start - 1), (intron_end + 1, end)]
                if strand == "-":
                    exons.reverse()  # parts are listed 5' to 3'
                trna = Feature("tRNA", exons, strand, attrs,
                               children=[Feature("exon", [part], strand) for part in exons])
                gene.children = [trna, Feature("intron", [(intron_start, intron_end)], strand)]
            else:
                # Anticodon: position is counted from the tRNA's 5' end
                pos = int(position)
                if strand == "+":
                    ac_start, ac_end = start + pos - 1, start + pos + 1
                else:
                    ac_start, ac_end = end - pos - 1, end - pos + 1
                if crosses_origin:
                    ac_start, ac_end = (ac_start - 1) % contig_length + 1, (ac_end - 1) % contig_length + 1
                location = f"complement({ac_start}..{ac_end})" if strand == "-" else f"{ac_start}..{ac_end}"
                attrs["anticodon"] = f"(pos:{location},aa:{aminoacid},seq:{codon.replace('(', '').replace(')', '')})"
                gene.children = [Feature("tRNA", [(start, end)], strand, attrs)]

        elif mol_type == "tmRNA":
            # Coding region, given relative to the tmRNA's 5' end
            rel_start, rel_end = map(int, position.split(","))
            if strand == "+":
                cds_start, cds_end = start + rel_start - 1, start + rel_end - 1
            else:
                cds_start, cds_end = end - rel_end + 1, end - rel_start + 1
            cds_attrs = {"transl_table": genetic_code} if genetic_code else {}
            cds = Feature("CDS", [(cds_start, cds_end)], strand, cds_attrs)
            gene.children = [Feature("tmRNA", [(start, end)], strand, children=[cds])]

        # Back to real coordinates, splitting whatever crosses the origin.
        # Parts stay 5' to 3': on the minus strand, the piece at the origin comes first.
        if crosses_origin:
            for item in [gene] + get_gene_features(gene):
                wrapped = []
                for s, e in item.parts:
                    if s > contig_length:
                        wrapped.append((s - contig_length, e - contig_length))
                    elif e > contig_length:
                        pieces = [(s, contig_length), (1, e - contig_length)]
                        wrapped += pieces if strand == "+" else pieces[::-1]
                    else:
                        wrapped.append((s, e))
                item.parts = wrapped

    return genes


def parse_hmmer(tbl_file, out_file, seed_evalue=1e-5, max_hole=50, max_overlap=150, max_gap=10_000):
    """
    Parses nhmmer output of the rnl/rns HMMs into a list of Gene objects.

    For each contig and gene, only the lineage model with the highest summed
    score is used. Its hits are chained into segments (the exons of one
    continuous stretch of the gene), and segments are grouped into gene copies:
    segments sharing model positions are separate copies (e.g. inverted
    repeats), segments that don't are modules of one fragmented gene, named
    rnl_a, rnl_b... in model order.

    Exon-intron junctions use alignment coordinates, which match curated
    boundaries closely; the outer ends of each gene or module use envelope
    coordinates, since alignments fade out at the variable rRNA termini.

    Parameters:
    -----------
    tbl_file : str or Path
        nhmmer --tblout file. May hold hits of several models and contigs.
    out_file : str or Path
        nhmmer main output (-o) of the same search, written with alignments
        (no --noali) and --notextw. Used to place the cut when consecutive
        exons overlap in the model.
    seed_evalue : float
        A segment must contain at least one hit with this E-value or lower.
        Weaker hits can only extend a segment, never start one.
    max_hole : int
        Max model positions skipped between consecutive exons.
    max_overlap : int
        Max model positions shared by consecutive exons, and by modules of
        the same fragmented gene.
    max_gap : int
        Max bases between consecutive exons.

    Returns:
    --------
    list of Gene
    """
    products = {
        "rnl": "large subunit ribosomal RNA",
        "rns": "small subunit ribosomal RNA",
    }

    def distance(a, b):
        """Bases between hit a and hit b, b being downstream of a on their strand."""
        if a["strand"] == "+":
            return b["start"] - a["end"] - 1
        return a["start"] - b["end"] - 1

    def follows(a, b):
        """True if hit b can be the next exon after hit a."""
        hole = b["hmm_from"] - a["hmm_to"] - 1  # negative = overlap in the model
        return (a["strand"] == b["strand"]
                and b["hmm_from"] > a["hmm_from"] and b["hmm_to"] > a["hmm_to"]
                and -max_overlap <= hole <= max_hole
                and 1 <= distance(a, b) <= max_gap)

    # ------------------------------------------------------------
    # 1. Read hits, grouped by contig and gene
    # ------------------------------------------------------------
    # Columns: target, accession, query, accession, hmmfrom, hmm to, alifrom,
    # ali to, envfrom, env to, sq len, strand, E-value, score, bias, description.
    # Minus-strand hits are written end..start.
    groups = defaultdict(list)
    with open(tbl_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            a, b, env_a, env_b = int(fields[6]), int(fields[7]), int(fields[8]), int(fields[9])
            gene_type = "rnl" if "rnl" in fields[2] else "rns"
            groups[(fields[0], gene_type)].append({
                "key": (fields[0], fields[2], min(a, b), max(a, b)),  # to find its alignment
                "model": fields[2],
                "hmm_from": int(fields[4]), "hmm_to": int(fields[5]),
                "start": min(a, b), "end": max(a, b),
                "env_start": min(env_a, env_b), "env_end": max(env_a, env_b),
                "strand": fields[11],
                "evalue": float(fields[12]), "score": float(fields[13]),
            })

    # ------------------------------------------------------------
    # 2. Read alignments, as one (model position, genome position, PP)
    #    tuple per column. Insertions have no model position, deletions no
    #    genome position. PP digits are 0-9 and "*" counts as 10.
    # ------------------------------------------------------------
    # Each alignment is four lines: model, match, target and PP:
    #   Rhodophyta_rnl_mito  2220 actcat...ggtc..cctA... 2434
    #                             actcat...ggtc   c A...
    #           NC_026905.1 25958 ACTCAT...GGTCgtGCGA... 26153
    #                             79****...***944589... PP
    alignments = {}
    with open(out_file) as f:
        lines = f.read().splitlines()
    query = None
    for i, line in enumerate(lines):
        fields = line.split()
        if line.startswith("Query:"):
            query = fields[1]
        elif len(fields) == 4 and fields[0] == query and fields[1].isdigit():
            model_line, target_line, pp_line = fields, lines[i + 2].split(), lines[i + 3].split()[0]
            a, b = int(target_line[1]), int(target_line[3])
            step = 1 if b >= a else -1
            model_pos, genome_pos = int(model_line[1]), a
            columns = []
            for m, t, p in zip(model_line[2], target_line[2], pp_line):
                pp = 10 if p == "*" else int(p) if p.isdigit() else 0
                columns.append((model_pos if m != "." else None,
                                genome_pos if t != "-" else None,
                                pp))
                if m != ".":
                    model_pos += 1
                if t != "-":
                    genome_pos += step
            alignments[(target_line[0], query, min(a, b), max(a, b))] = columns

    genes = []
    for (seqid, gene_type), hits in groups.items():

        # ------------------------------------------------------------
        # 3. Keep the lineage model with the highest summed score
        # ------------------------------------------------------------
        totals = defaultdict(float)
        for h in hits:
            if h["evalue"] <= seed_evalue:
                totals[h["model"]] += h["score"]
        if not totals:
            continue
        best_model = max(totals, key=totals.get)
        hits = sorted((h for h in hits if h["model"] == best_model),
                      key=lambda h: h["score"], reverse=True)

        # ------------------------------------------------------------
        # 4. Chain hits into segments, seeding from the strongest and
        #    always taking the nearest compatible hit
        # ------------------------------------------------------------
        unused = list(hits)
        segments = []
        for seed in hits:
            if seed not in unused or seed["evalue"] > seed_evalue:
                continue
            unused.remove(seed)
            segment = [seed]

            while True:  # towards 3'
                options = [h for h in unused if follows(segment[-1], h)]
                if not options:
                    break
                nxt = min(options, key=lambda h: distance(segment[-1], h))
                segment.append(nxt)
                unused.remove(nxt)

            while True:  # towards 5'
                options = [h for h in unused if follows(h, segment[0])]
                if not options:
                    break
                prv = min(options, key=lambda h: distance(h, segment[0]))
                segment.insert(0, prv)
                unused.remove(prv)

            segments.append(segment)

        # ------------------------------------------------------------
        # 5. Group segments: shared model positions = separate copies,
        #    no shared positions = modules of one fragmented gene
        # ------------------------------------------------------------
        copies = []
        for seg in sorted(segments, key=lambda s: sum(h["score"] for h in s), reverse=True):
            for copy in copies:
                overlaps = [min(seg[-1]["hmm_to"], other[-1]["hmm_to"])
                            - max(seg[0]["hmm_from"], other[0]["hmm_from"]) + 1
                            for other in copy]
                if max(overlaps) <= max_overlap:
                    copy.append(seg)
                    break
            else:
                copies.append([seg])

        for copy in copies:
            copy.sort(key=lambda s: s[0]["hmm_from"])
            for i, segment in enumerate(copy):
                name = gene_type if len(copy) == 1 else f"{gene_type}_{chr(ord('a') + i)}"
                strand = segment[0]["strand"]

                # ------------------------------------------------------------
                # 6. Resolve model overlaps between consecutive exons: prev
                #    keeps the shared positions up to the cut, hit keeps the
                #    rest. Take the cut with the highest summed PP of both
                #    alignments over the shared positions.
                # ------------------------------------------------------------
                for prev, hit in zip(segment, segment[1:]):
                    if prev["hmm_to"] < hit["hmm_from"]:
                        continue

                    prev_columns, hit_columns = alignments[prev["key"]], alignments[hit["key"]]
                    prev_pp = {m: p for m, g, p in prev_columns if m is not None}
                    hit_pp = {m: p for m, g, p in hit_columns if m is not None}
                    shared = range(hit["hmm_from"], prev["hmm_to"] + 1)
                    cut = max(range(shared.start - 1, shared.stop),
                              key=lambda c: sum(prev_pp[m] for m in shared if m <= c)
                                          + sum(hit_pp[m] for m in shared if m > c))

                    # Last genome base of prev at or before the cut, first of hit after it
                    prev_3 = [g for m, g, p in prev_columns if m is not None and m <= cut and g is not None][-1]
                    hit_5 = [g for m, g, p in hit_columns if m is not None and m > cut and g is not None][0]
                    prev["hmm_to"], hit["hmm_from"] = cut, cut + 1
                    if strand == "+":
                        prev["end"], hit["start"] = prev_3, hit_5
                    else:
                        prev["start"], hit["end"] = prev_3, hit_5

                # ------------------------------------------------------------
                # 7. Exons 5' to 3': envelope coords at the outer ends,
                #    alignment coords at the junctions
                # ------------------------------------------------------------
                exons = [[h["start"], h["end"]] for h in segment]
                if strand == "+":
                    exons[0][0] = segment[0]["env_start"]
                    exons[-1][1] = segment[-1]["env_end"]
                else:
                    exons[0][1] = segment[0]["env_end"]
                    exons[-1][0] = segment[-1]["env_start"]
                exons = [tuple(e) for e in exons]

                # ------------------------------------------------------------
                # 8. Gene > rRNA > exons and introns, as in MFannot
                # ------------------------------------------------------------
                children = []
                if len(exons) > 1:
                    for exon, next_exon in zip(exons, exons[1:]):
                        if strand == "+":
                            intron = (exon[1] + 1, next_exon[0] - 1)
                        else:
                            intron = (next_exon[1] + 1, exon[0] - 1)
                        children += [Feature("exon", [exon], strand), Feature("intron", [intron], strand)]
                    children.append(Feature("exon", [exons[-1]], strand))

                rrna = Feature("rRNA", exons, strand, {"product": products[gene_type]},
                               children, score=round(sum(h["score"] for h in segment), 1))
                span = (min(s for s, _ in exons), max(e for _, e in exons))
                genes.append(Gene(seqid, "nhmmer", name, [span], strand, {"gene": name}, [rrna]))

    return genes


def merge_annotations(mfannot_genes, aragorn_genes, hmmer_genes):
    """
    Combines MFannot, Aragorn and nhmmer genes. Aragorn tRNAs and nhmmer
    rnl/rns replace MFannot's, so those MFannot genes are dropped. Only the
    rRNA genes themselves are matched (rnl, rns, rnl_a...), so ORFs inside
    their introns, which MFannot writes as separate genes, are kept.
    """
    kept = [gene for gene in mfannot_genes
            if not gene.name.startswith("trn") and not re.fullmatch(r"rn[ls](_[a-z])?", gene.name)]
    return kept + aragorn_genes + hmmer_genes


def write_gff3(genes, path, contig_ids=None):
    """
    Writes genes to a GFF3 file, sorted by contig (in contig_ids order),
    then position, then feature type.

    Every gene gets a unique ID based on its name. Every feature gets an ID
    made of its gene's ID and its type, numbered when a gene has several of
    one type (cox1_intron1, cox1_intron2...). A feature with several parts
    is written as one row per part, all sharing the same ID.
    """

    def escape(value):
        """Percent-encodes the characters that have a meaning in GFF3 column 9."""
        value = str(value)
        for char, code in [("%", "%25"), (";", "%3B"), ("=", "%3D"), ("&", "%26"), (",", "%2C")]:
            value = value.replace(char, code)
        return value

    rows = []

    def add_rows(gene, feat_type, parts, strand, attrs, score=None):
        """Adds one row per part of a gene or feature."""
        # CDS phase depends on how much coding sequence came before (5' to 3')
        if feat_type == "CDS":
            phases, coding_length = [], 0
            for s, e in parts:
                phases.append((3 - coding_length % 3) % 3)
                coding_length += e - s + 1
        else:
            phases = ["."] * len(parts)
        column9 = ";".join(f"{key}={escape(val)}" for key, val in attrs.items())
        for (s, e), phase in zip(parts, phases):
            rows.append([gene.seqid, gene.source, feat_type, s, e,
                         "." if score is None else f"{score:g}", strand, phase, column9])

    # 1. Unique gene IDs. The contig is prepended when a name occurs on more
    #    than one contig, and _1, _2... is appended when it occurs more than
    #    once on the same contig.
    by_name = defaultdict(lambda: defaultdict(list))
    for gene in genes:
        by_name[gene.name][gene.seqid].append(gene)

    gene_ids = {}
    for name, per_contig in by_name.items():
        for seqid, copies in per_contig.items():
            copies.sort(key=lambda g: min(s for s, _ in g.parts))
            for n, gene in enumerate(copies, start=1):
                prefix = f"{seqid}_" if len(per_contig) > 1 else ""
                suffix = f"_{n}" if len(copies) > 1 else ""
                gene_ids[gene] = f"{prefix}{name}{suffix}"

    # 2. Rows of each gene and of the features inside it
    for gene in genes:
        gene_id = gene_ids[gene]
        add_rows(gene, "gene", gene.parts, gene.strand,
                 {"ID": gene_id, "Name": gene.name, **gene.attrs})

        type_counts = Counter(feat.type for feat in get_gene_features(gene))
        seen = Counter()

        # Go down the tree, parents before children, remembering each
        # feature's parent ID
        pending = [(feat, gene_id) for feat in gene.children]
        while pending:
            feat, parent_id = pending.pop(0)
            seen[feat.type] += 1
            number = seen[feat.type] if type_counts[feat.type] > 1 else ""
            feat_id = f"{gene_id}_{feat.type}{number}"
            attrs = {"ID": feat_id, "Name": f"{gene.name} {feat.type}{number}",
                     "Parent": parent_id, **feat.attrs}
            add_rows(gene, feat.type, feat.parts, feat.strand, attrs, feat.score)
            pending = [(child, feat_id) for child in feat.children] + pending

    # 3. Sort by contig, position and type
    contig_order = list(contig_ids or [])
    for gene in genes:
        if gene.seqid not in contig_order:
            contig_order.append(gene.seqid)
    rows.sort(key=lambda r: (contig_order.index(r[0]), r[3], TYPE_ORDER.get(r[2], len(TYPE_ORDER))))

    # 4. Write
    with open(path, "w") as f:
        f.write("##gff-version 3\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")


def genbank_features(genes):
    """
    Converts genes into Biopython SeqFeatures, grouped by contig and sorted
    by position. Every feature inherits its gene's /name and /gene qualifiers.

    Returns:
    --------
    dict
        seqid -> list of SeqFeature.
    """

    def make_seqfeature(feat_type, parts, strand, attrs):
        """Builds one SeqFeature; several parts become a join()."""
        sign = 1 if strand == "+" else -1
        locations = [FeatureLocation(s - 1, e, strand=sign) for s, e in parts]
        location = locations[0] if len(locations) == 1 else CompoundLocation(locations)
        qualifiers = {key: [str(val)] for key, val in attrs.items()}
        return SeqFeature(location, type=feat_type, qualifiers=qualifiers)

    per_contig = defaultdict(list)
    for gene in genes:
        # Every feature carries its gene's name, so ORFs (which get no /gene)
        # are still labelled. /name is not an INSDC qualifier: an NCBI
        # submission output will need /locus_tag or similar instead.
        gene_qualifiers = {"name": gene.name}
        if "gene" in gene.attrs:
            gene_qualifiers["gene"] = gene.attrs["gene"]
        per_contig[gene.seqid].append(
            make_seqfeature("gene", gene.parts, gene.strand, {**gene_qualifiers, **gene.attrs}))
        for feat in get_gene_features(gene):
            per_contig[gene.seqid].append(
                make_seqfeature(feat.type, feat.parts, feat.strand, {**gene_qualifiers, **feat.attrs}))

    for features in per_contig.values():
        features.sort(key=lambda sf: (int(sf.location.start), TYPE_ORDER.get(sf.type, len(TYPE_ORDER))))
    return per_contig


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
    Main pipeline: Annotates a FASTA or GenBank sequence file using nhmmer
    (rnl, rns), MFannot and Aragorn. Circular sequences are rotated to start
    at rns.

    organelle is an organelle name ('chloroplast', 'plastid', 'mitochondrion',
    'mitochondria') or a genetic code number (e.g. 4 or "4").
    format is 'gb' or 'gff3'.
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
    if format not in ("gb", "gff3"):
        raise ValueError(f"Unsupported output format '{format}'. Use 'gb' or 'gff3'.")
    is_gb_input = input_suffix in valid_gb

    # Genetic code: a number is used as given, organelle names are translated
    organelle_value = str(organelle).strip().lower()
    if organelle_value.isdigit():
        genetic_code = int(organelle_value)
    elif organelle_value in ("chloroplast", "plastid"):
        genetic_code = 11
    elif organelle_value in ("mitochondrion", "mitochondria"):
        genetic_code = 4
    else:
        raise ValueError(
            f"Unknown organelle '{organelle}' for {seq_file}. Use 'chloroplast', 'plastid', "
            f"'mitochondrion', 'mitochondria' or a genetic code number."
        )

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
            f"  started:      {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            f"  organelle:    {organelle}\n"
            f"  genetic code: {genetic_code}\n"
            f"  format:       {format}\n"
            f"  circular:     {circular}\n"
            f"  keep_old:     {keep_old}\n"
            f"  output_dir:   {out_dir}"
        )
    else:
        log(f"Started: {time.strftime('%H:%M:%S')} | Organelle: {organelle} "
            f"| Genetic code: {genetic_code} | Circular: {circular}")

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
            contig_lengths = {cid: len(rec.seq) for cid, rec in records_dict.items()}
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
            # 4. rRNA GENES (nhmmer) & ROTATION TO RNS START
            # ---------------------------------------------------------
            # Circular contigs are rotated so that the 5' end of rns (rns_a if
            # fragmented) becomes position 1 on the plus strand. nhmmer then
            # runs again on the rotated sequences, so its genes use the same
            # coordinates as MFannot's and Aragorn's.
            hmm_tbl = tmp_path / f"{seq_stem}_nhmmer.tbl"
            hmm_out = tmp_path / f"{seq_stem}_nhmmer.out"
            log(annotate_hmmer(str(working_fasta), tmp_path, genetic_code), "nhmmer")
            hmm_genes = parse_hmmer(hmm_tbl, hmm_out)

            reordered = False
            for cid in contig_ids:
                if not is_circular[cid]:
                    continue

                hit = find_rns_start(hmm_genes, cid)
                if hit is None:
                    log(f"No rns found in contig {cid}, proceeding without reordering")
                    continue

                strand, start = hit
                if strand == "+" and start == 1:
                    log(f"rns already at position 1 in contig {cid}, no rotation needed")
                    continue

                cut_pos = start - 1 if strand == "+" else start
                new_seq_str = reorder_sequence(str(records_dict[cid].seq), strand=strand, cut_pos=cut_pos)

                records_dict[cid].seq = Seq(new_seq_str)
                if is_gb_input and cid in gb_records:
                    gb_records[cid].seq = Seq(new_seq_str)

                reordered = True
                log(f"Rotated contig {cid} to start at rns (strand {strand}, position {start})")

            if reordered:
                SeqIO.write(list(records_dict.values()), str(working_fasta), "fasta")
                log(annotate_hmmer(str(working_fasta), tmp_path, genetic_code), "nhmmer (rotated)")
                hmm_genes = parse_hmmer(hmm_tbl, hmm_out)

            # ---------------------------------------------------------
            # 5. FEATURE ANNOTATION (MFannot & Aragorn)
            # ---------------------------------------------------------
            try:
                mfannot_out = annotate_mfannot(str(working_fasta), tmp_path, genetic_code)
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
                log(annotate_aragorn(str(circ_fasta), tmp_path, genetic_code, circular=True), "Aragorn (circular)")
                ar_genes = parse_aragorn(tmp_path / f"{seq_stem}_circ.txt", genetic_code, seq_stem, contig_lengths)

                lin_fasta = tmp_path / f"{seq_stem}_lin.fasta"
                SeqIO.write(lin_records, str(lin_fasta), "fasta")
                log(annotate_aragorn(str(lin_fasta), tmp_path, genetic_code, circular=False), "Aragorn (linear)")
                ar_genes += parse_aragorn(tmp_path / f"{seq_stem}_lin.txt", genetic_code, seq_stem, contig_lengths)
            else:
                log(annotate_aragorn(str(working_fasta), tmp_path, genetic_code, circular=has_circular), "Aragorn")
                ar_genes = parse_aragorn(tmp_path / f"{seq_stem}.txt", genetic_code, seq_stem, contig_lengths)

            tbl_file = tmp_path / f"{seq_stem}.tbl"
            mf_genes = parse_mfannot(tbl_file, genetic_code, contig_ids, seq_stem)

            if not mf_genes and not ar_genes and not hmm_genes:
                log(f"File {seq_stem} has no features")
                return "\n".join(log_lines) if return_log else None

            genes = merge_annotations(mf_genes, ar_genes, hmm_genes)

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
                        rec.annotations["topology"] = "circular" if is_circular.get(rec.id, False) else "linear"
                        records[rec.id] = rec

                features = genbank_features(genes)
                for rec_id, rec in records.items():
                    rec.features = features.get(rec_id, [])

                SeqIO.write([records[i] for i in contig_ids if i in records], str(out_dir / f"{seq_stem}.gb"), "gb")

            elif format == "gff3":
                write_gff3(genes, out_dir / f"{seq_stem}.gff3", contig_ids)

            # ---------------------------------------------------------
            # 7. CLEANUP & FINISH
            # ---------------------------------------------------------
            if keep_intermediates:
                for filename in [f"{seq_stem}.tbl", f"{seq_stem}.txt", f"{seq_stem}_circ.txt", f"{seq_stem}_lin.txt",
                                 f"{seq_stem}_nhmmer.tbl", f"{seq_stem}_nhmmer.out"]:
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