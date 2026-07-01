"""Command-line interface for ARGWAS."""

import argparse
import sys

from .clump import prepare_plink_clump_file
from .pipeline import run_argwas


def add_run_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments for the main ARGWAS scan."""

    parser.add_argument("-o", "--output", required=True, help="output directory for results")
    parser.add_argument("-p", "--ped", required=True, help="input PLINK PED file")
    parser.add_argument(
        "-m",
        "--map",
        dest="map_file",
        default=None,
        help="input PLINK MAP file matching the PED; defaults to the PED prefix with .map suffix",
    )
    parser.add_argument("-t", "--trees", required=True, help="directory containing .trees files")
    parser.add_argument("-a", "--ances", default=None, help="ancestral alleles file; required with --ancestral-mode file")
    parser.add_argument(
        "--ancestral-mode",
        choices=["file", "arg"],
        default="file",
        help="how to obtain ancestral alleles for observed SNPs",
    )
    parser.add_argument(
        "--ambiguous-ancestral",
        choices=["skip", "major", "reference"],
        default="skip",
        help="fallback policy for SNPs whose ancestral allele cannot be inferred cleanly from the ARG",
    )
    parser.add_argument(
        "--ancestral-coordinate-system",
        choices=["one-based", "zero-based"],
        default="one-based",
        help="coordinate system used by MAP physical positions when indexing the ancestral sequence",
    )
    parser.add_argument("-w", "--workers", type=int, default=4, help="number of worker processes")


def build_run_parser() -> argparse.ArgumentParser:
    """Build the legacy/default parser for running ARGWAS directly."""

    parser = argparse.ArgumentParser(
        prog="argwas",
        description="ARGWAS: scan ARG-derived putative QTLs for phenotype association.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_run_arguments(parser)
    return parser


def build_parser() -> argparse.ArgumentParser:
    """Build the full ARGWAS parser with subcommands."""

    parser = argparse.ArgumentParser(
        prog="argwas",
        description="ARGWAS command-line tools.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser(
        "run",
        help="run the ARGWAS scan",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    add_run_arguments(run_parser)

    clump_parser = subparsers.add_parser(
        "prepare-clump",
        help="prepare results.csv for PLINK --clump",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    clump_parser.add_argument("results_csv", help="ARGWAS results.csv file")
    clump_parser.add_argument("-o", "--output", required=True, help="output association file")
    clump_parser.add_argument(
        "--marker-column",
        default="marker_id",
        help="column containing marker IDs matching PLINK MAP/BIM SNP IDs",
    )
    clump_parser.add_argument("--r-squared-column", default="r_squared", help="ARGWAS R² column")

    return parser


def main(argv=None):
    """Parse command-line arguments and run the selected ARGWAS command."""

    argv = list(sys.argv[1:] if argv is None else argv)

    # Backward compatibility: keep allowing the old form:
    #   argwas -o results -p input.ped -t trees -a ancestral.txt
    # The subcommand form is preferred for documentation:
    #   argwas run -o results -p input.ped -t trees -a ancestral.txt
    if argv and argv[0] not in {"run", "prepare-clump", "-h", "--help"}:
        args = build_run_parser().parse_args(argv)
        run_argwas(
            output=args.output,
            ped=args.ped,
            trees=args.trees,
            ances=args.ances,
            workers=args.workers,
            map_file=args.map_file,
            ancestral_coordinate_system=args.ancestral_coordinate_system,
            ancestral_mode=args.ancestral_mode,
            ambiguous_ancestral=args.ambiguous_ancestral,
        )
        return

    args = build_parser().parse_args(argv)

    if args.command == "run":
        run_argwas(
            output=args.output,
            ped=args.ped,
            trees=args.trees,
            ances=args.ances,
            workers=args.workers,
            map_file=args.map_file,
            ancestral_coordinate_system=args.ancestral_coordinate_system,
            ancestral_mode=args.ancestral_mode,
            ambiguous_ancestral=args.ambiguous_ancestral,
        )
    elif args.command == "prepare-clump":
        output = prepare_plink_clump_file(
            results_csv=args.results_csv,
            output=args.output,
            marker_column=args.marker_column,
            r_squared_column=args.r_squared_column,
        )
        print(f"Wrote PLINK clumping association file: {output}")
    else:  # pragma: no cover - argparse should prevent this.
        raise ValueError(f"Unknown command: {args.command}")


if __name__ == "__main__":
    main()
