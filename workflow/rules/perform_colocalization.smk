# Run coloc on a pair of extracted regions
rule run_coloc:
    input:
        input_signals=config["DATA_DIR"] + "/input_signals.csv",
        input_pairs=config["OUTPUTS_DIR"] + "/pairs_to_test.csv",
        input_joined_sumstats=config["OUTPUTS_DIR"]
        + "/extracted_regions/{pair_id}_sumstats.txt",
    params:
        pair_id="{pair_id}",
    output:
        output_results=config["OUTPUTS_DIR"] + "/coloc_results/{pair_id}.txt",
        output_cs=config["OUTPUTS_DIR"] + "/coloc_credible_sets/{pair_id}_cs95.txt",
    script:
        "../scripts/run_coloc.R"


rule aggregate_coloc_results:
    input:
        expand(config["OUTPUTS_DIR"] + "/coloc_results/{pair_id}.txt", pair_id=PAIRS),
    output:
        aggregated_results=config["OUTPUTS_DIR"] + "/all_coloc_results.txt",
    shell:
        "awk 'FNR == 2' {input} > {output.aggregated_results}"


rule define_colocalized_signals:
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
    output:
        output_joined_cs=config["OUTPUTS_DIR"]
        + "/joined_credible_sets/{signal_code}_cs95.txt",
    script:
        "../scripts/join_credible_sets.R"
