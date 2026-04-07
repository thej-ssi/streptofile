#!/usr/bin/env python3

from pathlib import Path
import polars as pl
import subprocess
import argparse
from importlib import resources
from wgs_data_manager.wgs_data_manager import WgsSample, WgsCollection


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input",
                        nargs="*",
                        help = "Input fasta file(s)",
                        type = Path)
    parser.add_argument("-o", "--output",
                        help = "Output directory.",
                        type=Path,
                        required = True)
    parser.add_argument("-d", "--database",
                        help = "Path to fasta file with emm allele sequences",
                        type=Path,
                        default = resources.files("streptofile") / "db" / "alltrimmed.tfa")
    parser.add_argument("--full_path",
                        help = "Print full path to fasta input in output tsv rather than just file name",
                        action= "store_true",
                        default=False)
    return parser.parse_args()

def run_emm_blast(assembly_file: Path, emm_allele_fasta: Path, output_file: Path) -> tuple(str) | None:
    output_file = Path(output_file)
    #if not output_file.exists():
    cmd = f'blastn -query {emm_allele_fasta} -subject {assembly_file} -qcov_hsp_perc 90 -out {output_file} -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send sseq evalue bitscore"'
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        return False
    else:
        return True
    #else:
    #    return None
    
def extract_emm_type(emm_blast_tsv: Path):
    emm_types_in_emm_plus_mrp_operons = []  # to update
    mrp_types_in_emm_plus_mrp_operons = [
        "134",
        "156",
        "159",
        "164",
        "174",
        "205",
    ]  # to update

    emm_blast_tsv = Path(emm_blast_tsv)
    emm_typing_results = {"EMM_type": "-", "ENN_type": "-", "MRP_type": "-"}

    if not emm_blast_tsv.exists():
        emm_typing_results["emm_typing_notes"] = "No blast output found for EMM genes"
        return emm_typing_results

    try:
        blast_df = pl.read_csv(
            emm_blast_tsv,
            separator="\t",
            has_header=False,
            new_columns=[
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "qlen",
                "qstart",
                "qend",
                "sstart",
                "send",
                "sseq",
                "evalue",
                "bitscore",
            ],
        )
    except pl.exceptions.NoDataError:
        emm_typing_results["emm_typing_notes"] = "Empty blast output, no EMM genes detected"
        return emm_typing_results

    notes = []

    blast_df = (
        blast_df
        .with_columns([
            (pl.col("length") / pl.col("qlen") * 100).alias("plen"),
            pl.when(pl.col("sstart") < pl.col("send"))
            .then(((pl.col("sstart") - pl.col("qstart") + 1) / 100).round(0))
            .otherwise(((pl.col("send") - pl.col("qstart") + 1) / 100).round(0))
            .cast(pl.Int64)
            .alias("extended_sstart"),
        ])
        .filter(pl.col("bitscore") > 200)
    )

    blast_df_unique = (
        blast_df
        .sort("bitscore", descending=True)
        .group_by("extended_sstart")
        .first()
    )

    if blast_df_unique.height == 0:
        notes.append("No blast hits found for EMM genes")

    elif blast_df_unique["sseqid"].n_unique() == 1:
        if blast_df_unique.height == 1:
            row0 = blast_df_unique.row(0, named=True)
            emm_typing_results["EMM_type"] = "EMM" + row0["qseqid"][3:]
            if row0["length"] < row0["qlen"] or row0["pident"] < 100:
                emm_typing_results["EMM_type"] += "*"
                notes.append(
                    f"EMM{row0['qseqid'][3:]} with {round(row0['pident'], 2)} "
                    f"and length {row0['length']}/{row0['qlen']}"
                )
        else:
            first_row = blast_df_unique.row(0, named=True)
            if first_row["sstart"] < first_row["send"]:
                blast_df_unique = blast_df_unique.sort("sstart")
            else:
                blast_df_unique = blast_df_unique.sort("sstart", descending=True)

            rows = list(blast_df_unique.iter_rows(named=True))

            if blast_df_unique.height == 2:
                row0, row1 = rows

                emm_typing_results["EMM_type"] = "EMM" + row0["qseqid"][3:]
                if row0["length"] < row0["qlen"] or row0["pident"] < 100:
                    emm_typing_results["EMM_type"] += "*"
                    notes.append(
                        f"EMM{row0['qseqid'][3:]} with pident {round(row0['pident'], 2)} "
                        f"and length {row0['length']}/{row0['qlen']}"
                    )

                emm_typing_results["ENN_type"] = "EMM" + row1["qseqid"][3:]
                if row1["length"] < row1["qlen"] or row1["pident"] < 100:
                    emm_typing_results["ENN_type"] += "*"
                    notes.append(
                        f"ENN{row1['qseqid'][3:]} with pident {round(row1['pident'], 2)} "
                        f"and length {row1['length']}/{row1['qlen']}"
                    )

                emm_maintype = row0["qseqid"][3:].split(".")[0]
                mrp_maintype = row1["qseqid"][3:].split(".")[0]

                if (
                    mrp_maintype in emm_types_in_emm_plus_mrp_operons
                    or emm_maintype in mrp_types_in_emm_plus_mrp_operons
                ):
                    emm_typing_results["MRP_type"] = "EMM" + emm_typing_results["EMM_type"][3:]
                    emm_typing_results["EMM_type"] = "EMM" + emm_typing_results["ENN_type"][3:]
                    emm_typing_results["ENN_type"] = "-"
                    notes.append("EMM redesignated due to known MRP+EMM operon")

            elif blast_df_unique.height == 3:
                row0, row1, row2 = rows

                emm_typing_results["MRP_type"] = "EMM" + row0["qseqid"][3:]
                if row0["length"] < row0["qlen"] or row0["pident"] < 100:
                    emm_typing_results["MRP_type"] += "*"
                    notes.append(
                        f"MRP{row0['qseqid'][3:]} with pident {round(row0['pident'], 2)} "
                        f"and length {row0['length']}/{row0['qlen']}"
                    )

                emm_typing_results["EMM_type"] = "EMM" + row1["qseqid"][3:]
                if row1["length"] < row1["qlen"] or row1["pident"] < 100:
                    emm_typing_results["EMM_type"] += "*"
                    notes.append(
                        f"EMM{row1['qseqid'][3:]} with pident {round(row1['pident'], 2)} "
                        f"and length {row1['length']}/{row1['qlen']}"
                    )

                emm_typing_results["ENN_type"] = "EMM" + row2["qseqid"][3:]
                if row2["length"] < row2["qlen"] or row2["pident"] < 100:
                    emm_typing_results["ENN_type"] += "*"
                    notes.append(
                        f"ENN{row2['qseqid'][3:]} with pident {round(row2['pident'], 2)} "
                        f"and length {row2['length']}/{row2['qlen']}"
                    )

    else:
        emm_genes = []
        for row in blast_df_unique.iter_rows(named=True):
            if row["length"] < row["qlen"] or row["pident"] < 100:
                emm_genes.append(row["qseqid"][3:] + "*")
            else:
                emm_genes.append(row["qseqid"][3:])
        notes.append(
            "EMM and EMM-like genes found on multiple contigs. Alleles found: "
            + "/".join(emm_genes)
        )

    emm_typing_results["emm_typing_notes"] = ", ".join(notes)
    return emm_typing_results

def emm_typer(assembly_file: Path,
              output_folder: Path,
              emm_allele_fasta: Path) -> dict:
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok= True)
    blast_output_tsv = output_folder.joinpath("emm_blast.tsv")
    blast_complete = run_emm_blast(assembly_file = assembly_file,
                                   emm_allele_fasta = emm_allele_fasta,
                                   output_file=blast_output_tsv)
    if blast_complete:
        emm_results_dict = extract_emm_type(emm_blast_tsv = blast_output_tsv)
    else:
        emm_results_dict = {"EMM_type": "-", "ENN_type": "-", "MRP_type": "-","emm_typing_notes": "Failed to run blast"}
    return emm_results_dict

def print_emm_results(emm_results_dicts: list[dict],
                      output_file: Path) -> None:
    emm_results_df = pl.from_dict(emm_results_dicts)
    emm_results_df.write_csv(file = output_file,
                             separator = "\t")
    return None


def main():
    args = parse_args()
    emm_result_dicts = []
    for assembly_file in args.input:
        sample_output_dir = args.output / assembly_file.stem
        if args.full_path:
            emm_results = {"sample": str(assembly_file)}
        else:
            emm_results = {"sample": assembly_file.name}
        emm_results.update(emm_typer(assembly_file=assembly_file, output_folder=sample_output_dir, emm_allele_fasta=args.database))
        emm_result_dicts.append(emm_results)
    print(emm_result_dicts)
    

if __name__ == "__main__":
    main()
