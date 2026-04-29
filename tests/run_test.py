
from importlib import resources
from streptofile import emm_typer, mlstyper, virulence_profiler
from pathlib import Path

test_output_dir = resources.files("streptofile") / "test_output"
mlst_profiles_tsv = resources.files("streptofile") / "db" / "mlst" / "profiles.tsv"


#### Tests for proper parsing of EMM blast output

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

def test_emm_6() -> None:
    emm_results = emm_typer.extract_emm_type(test_output_dir / "GCA_003240915.2_ASM324091v2_genomic" / "emm_blast.tsv")
    print(emm_results.row(0))
    assert(emm_results.row(0) == ('-', '-', '-', '-', 'EMM blast output empty'))


#### Tests for proper parsing of MLST blast output

def test_mlst_1() -> None:
    mlst_results = mlstyper.extract_mlst_type(mlst_blast_tsv=test_output_dir / "GCF_000011765.3_ASM1176v2_genomic" / "mlst_blast.tsv",
                                              profiles_tsv=mlst_profiles_tsv)
    print(mlst_results.row(0))
    assert(mlst_results.row(0) == ('28', '4', '3', '4', '4', '4', '2', '4', 'All loci matched exactly'))


def test_mlst_2() -> None:
    mlst_results = mlstyper.extract_mlst_type(mlst_blast_tsv=test_output_dir / "GCA_034705955.1_PDT001983156.1_genomic" / "mlst_blast.tsv",
                                              profiles_tsv=mlst_profiles_tsv)
    print(mlst_results.row(0))
    assert(mlst_results.row(0) == ('ND', '3', '4', '6', '6', '14', '5', '-', 'Missing MLST genes: yqiL'))

def test_mlst_3() -> None:
    mlst_results = mlstyper.extract_mlst_type(mlst_blast_tsv=test_output_dir / "GCA_037112455.1_PDT002108051.1_genomic" / "mlst_blast.tsv",
                                              profiles_tsv=mlst_profiles_tsv)
    print(mlst_results.row(0))
    assert(mlst_results.row(0) == ('624', '38', '2', '2', '23', '39', '2', '1', 'All loci matched exactly'))

def test_mlst_4() -> None:
    mlst_results = mlstyper.extract_mlst_type(mlst_blast_tsv=test_output_dir / "GCF_001535505.1_ASM153550v1_genomic" / "mlst_blast.tsv",
                                              profiles_tsv=mlst_profiles_tsv)
    print(mlst_results.row(0))
    assert(mlst_results.row(0) == ('39', '5', '11', '8', '5', '15', '2', '1', 'All loci matched exactly'))

def test_mlst_5() -> None:
    mlst_results = mlstyper.extract_mlst_type(mlst_blast_tsv=test_output_dir / "GCF_002236855.1_ASM223685v1_genomic" / "mlst_blast.tsv",
                                              profiles_tsv=mlst_profiles_tsv)
    print(mlst_results.row(0))
    assert(mlst_results.row(0) == ('1065', '2', '2', '125', '3', '18', '24', '2', 'All loci matched exactly'))

def test_mlst_6() -> None:
    mlst_results = mlstyper.extract_mlst_type(mlst_blast_tsv=test_output_dir / "GCA_003240915.2_ASM324091v2_genomic/" / "mlst_blast.tsv",
                                              profiles_tsv=mlst_profiles_tsv)
    print(mlst_results.row(0))
    assert(mlst_results.row(0) == ('ND', '-', '-', '-', '-', '-', '-', '-', 'MLST blast output empty'))


test_emm_1()
test_emm_2()
test_emm_3()
test_emm_4()
test_emm_5()
test_emm_6()
test_mlst_1()
test_mlst_2()
test_mlst_3()
test_mlst_4()
test_mlst_5()
test_mlst_6()