import argparse
from pipeline import annotate_fasta

def main():
    parser = argparse.ArgumentParser(
        description="Annotate seaweed organelle genomes using MFannot and Aragorn."
    )
    parser.add_argument("fasta_file", help="Input FASTA filename (must be inside --directory)")
    parser.add_argument(
        "--organelle", required=True, choices=["chloroplast", "mitochondrion"],
        help="Organelle type, determines genetic code"
    )
    parser.add_argument(
        "--format", default="gb", choices=["gb", "gff3"],
        help="Output format (default: gb)"
    )
    parser.add_argument(
        "--directory", default="/data",
        help="Working directory containing for the input and output file (default: /data)"
    )
    parser.add_argument(
        "--circular", action="store_true",
        help="Treat the genome as circular and attempt rotation to start at rns"
    )
    parser.add_argument(
        "--no-keep-old", dest="keep_old", action="store_false",
        help="Discard the original (pre-rotation) fasta instead of keeping it. Only when exporting in gff3"
    )
    parser.add_argument(
        "--keep-intermediates", action="store_true",
        help="Also write the MFannot table (.tbl) and the Aragorn output (.txt)"
    )

    args = parser.parse_args()

    annotate_fasta(
        args.fasta_file,
        format=args.format,
        directory=args.directory,
        organelle=args.organelle,
        circular=args.circular,
        keep_old=args.keep_old,
        keep_intermediates=args.keep_intermediates,
    )

if __name__ == "__main__":
    main()