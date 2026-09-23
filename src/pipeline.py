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
# The parsers turn MFannot and Aragorn output into a list of Gene objects,
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


def parse_aragorn(txt_file, genetic_code=None, seq_name="Default_name"):
    """
    Parses Aragorn batch (-w) output into a list of Gene objects.

    Parameters:
    -----------
    txt_file : str or Path
        Path to the Aragorn output file.
    genetic_code : int, optional
        Written as transl_table on the tmRNA coding region.
    seq_name : str
        Fallback sequence identifier, used only if a data row precedes any
        contig header, which does not happen in normal Aragorn output.

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
    #   >contig_id len=...
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

    return genes


# ============================================================
# Merging
# ============================================================

def merge_annotations(mfannot_genes, aragorn_genes):
    """
    Combines MFannot and Aragorn genes. Aragorn tRNA predictions replace
    MFannot's, so MFannot tRNA genes are dropped.
    """
    kept = [gene for gene in mfannot_genes if not gene.name.startswith("trn")]
    return kept + aragorn_genes


# ============================================================
# GFF3 output
# ============================================================

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


# ============================================================
# GenBank output
# ============================================================

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


# ============================================================
# Main pipeline
# ============================================================

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
                ar_genes = parse_aragorn(tmp_path / f"{seq_stem}_circ.txt", genetic_code, seq_stem)

                lin_fasta = tmp_path / f"{seq_stem}_lin.fasta"
                SeqIO.write(lin_records, str(lin_fasta), "fasta")
                log(annotate_aragorn(str(lin_fasta), tmp_path, genetic_code, circular=False), "Aragorn (linear)")
                ar_genes += parse_aragorn(tmp_path / f"{seq_stem}_lin.txt", genetic_code, seq_stem)
            else:
                log(annotate_aragorn(str(working_fasta), tmp_path, genetic_code, circular=has_circular), "Aragorn")
                ar_genes = parse_aragorn(tmp_path / f"{seq_stem}.txt", genetic_code, seq_stem)

            tbl_file = tmp_path / f"{seq_stem}.tbl"
            mf_genes = parse_mfannot(tbl_file, genetic_code, contig_ids, seq_stem)

            if not mf_genes and not ar_genes:
                log(f"File {seq_stem} has no features")
                return "\n".join(log_lines) if return_log else None

            genes = merge_annotations(mf_genes, ar_genes)

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