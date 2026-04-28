
from importlib import resources
from streptofile import emm_typer, mlstyper, virulence_profiler
from pathlib import Path

samplesheet_tsv = resources.files("streptofile") / "test_output" / "test_output" / "GCA_034705955.1_PDT001983156.1_genomic" / "emm_blast.tsv"
assembly_dir = resources.files("wgs_data_manager") / "test" / "assembly_files"
paired_end_reads_dir = resources.files("wgs_data_manager") / "test" / "paired_end_read_files"
test_output_dir = resources.files("streptofile") / "test_output"



def get_test_input(test_output_dir) -> tuple[list[Path], list[Path], list[Path]]:
    emm_blast_outputs = []
    mlst_blast_outputs = []
    virulence_blast_outputs = []
    return None



def test_mlst_1() -> None:
    return None

def test_emm_1() -> None:
    emm_results = emm_typer.extract_emm_type(test_output_dir / "GCF_000011765.3_ASM1176v2_genomic" / "emm_blast.tsv")
    print(emm_results.row(0))
    assert(emm_results.row(0) == ('1.0', '-', '-', '1.0', 'All exact allele matches'))

def test_emm_2() -> None:
    emm_results = emm_typer.extract_emm_type(test_output_dir / "GCA_034705955.1_PDT001983156.1_genomic" / "emm_blast.tsv")
    assert(emm_results.row(0) == ('11.0', '202.1', '141.2', '141.2,11.0,202.1', 'All exact allele matches'))


def test_emm_3() -> None:
    emm_results = emm_typer.extract_emm_type(test_output_dir / "GCA_037112455.1_PDT002108051.1_genomic" / "emm_blast.tsv")
    print(emm_results.row(0))
    assert(emm_results.row(0) == ('81.0', '164.3', '156.4*', '156.4,81.0,164.3', 'MRP156.4__180__99.44'))

def test_emm_4() -> None:
    emm_results = emm_typer.extract_emm_type(test_output_dir / "GCF_001535505.1_ASM153550v1_genomic" / "emm_blast.tsv")
    print(emm_results.row(0))
    assert(emm_results.row(0) == ('4.0', '-', '156.0', '156.0,4.0', 'EMM redesignated due to known MRP+EMM operon'))

def test_emm_5() -> None:
    emm_results = emm_typer.extract_emm_type(test_output_dir / "GCF_002236855.1_ASM223685v1_genomic" / "emm_blast.tsv")
    print(emm_results.row(0))
    assert(emm_results.row(0) == ('-', '-', '-', '156.4,111.2,111.2,29.32', 'More than three EMM genes found, 156.4__177__98.87, 111.2__180__99.44, 111.2__165__99.39, 29.32__167__98.2'))


test_emm_1()
test_emm_2()
test_emm_3()
test_emm_4()
test_emm_5()