import pandas as pd
import os


# Run coloc on a pair of extracted regions
rule run_coloc:
    input:
        input_signals=config["DATA_DIR"] + "/input_signals.csv",
        input_pairs=config["OUTPUTS_DIR"] + "/pairs_to_test.csv",
        input_joined_sumstats=config["OUTPUTS_DIR"]
        + "/extracted_regions/{pair_id}_sumstats.txt",
    params:
        pair_id="{pair_id}",
    log:
        config["OUTPUTS_DIR"] + "/logs/run_coloc_{pair_id}.log",
    output:
        output_results=config["OUTPUTS_DIR"] + "/coloc_results/{pair_id}.txt",
        output_cs=config["OUTPUTS_DIR"] + "/coloc_credible_sets/{pair_id}_cs95.txt",
    script:
        "../scripts/run_coloc.R"


rule aggregate_coloc_results:
    input:
        pairs_file=config["OUTPUTS_DIR"] + "/pairs_to_test.csv",
        coloc_results=lambda wildcards: expand(
            config["OUTPUTS_DIR"] + "/coloc_results/{pair_id}.txt",
            pair_id=get_pairs_to_test(wildcards),
        ),
    log:
        config["OUTPUTS_DIR"] + "/logs/aggregate_coloc_results.log",
    output:
        aggregated_results=config["OUTPUTS_DIR"] + "/all_coloc_results.txt",
    shell:
        "awk 'FNR == 2' {input.coloc_results} > {output.aggregated_results}"


# Checkpoint for defining colocalized signals
checkpoint define_colocalized_signals:
    input:
        input_coloc_results=config["OUTPUTS_DIR"] + "/all_coloc_results.txt",
        input_pairs=config["OUTPUTS_DIR"] + "/pairs_to_test.csv",
    params:
        dense_nvar_threshold=config["DENSE_NVAR_THRESHOLD"],
        signal_dropout_threshold=config["SIGNAL_DROPOUT_THRESHOLD"],
        coloc_h4_threshold=config["COLOC_H4_THRESHOLD"],
        sensitivity_p12_threshold=config["SENSITIVITY_P12_THRESHOLD"],
        min_n=config["MIN_N"],
        min_dense=config["MIN_DENSE"],
        max_dropout=config["MAX_DROPOUT"],
        min_coloc=config["MIN_COLOC"],
        min_robust=config["MIN_ROBUST"],
    log:
        config["OUTPUTS_DIR"] + "/logs/define_colocalized_signals.log",
    output:
        output_coloc_signals=config["OUTPUTS_DIR"] + "/colocalized_signals.txt",
        output_coloc_codes=config["OUTPUTS_DIR"] + "/colocalized_signal_codes.txt",
    script:
        "../scripts/define_colocalized_signals.R"


rule join_credible_sets:
    input:
        input_coloc_codes=config["OUTPUTS_DIR"] + "/colocalized_signal_codes.txt",
        input_cs_dir=config["OUTPUTS_DIR"] + "/coloc_credible_sets/",
    params:
        signal_code="{signal_code}",
    log:
        config["OUTPUTS_DIR"] + "/logs/{signal_code}_join_credible_sets.log",
    output:
        output_joined_cs=config["OUTPUTS_DIR"]
        + "/joined_credible_sets/{signal_code}_cs95.txt",
    script:
        "../scripts/join_credible_sets.R"


rule all:
    input:
        config["OUTPUTS_DIR"] + "/all_coloc_results.txt",
        config["OUTPUTS_DIR"] + "/colocalized_signals.txt",
        config["OUTPUTS_DIR"] + "/colocalized_signal_codes.txt",
        lambda wildcards: expand(
            config["OUTPUTS_DIR"] + "/joined_credible_sets/{signal_code}_cs95.txt",
            signal_code=get_coloc_ids_from_checkpoint(wildcards),
        ),
