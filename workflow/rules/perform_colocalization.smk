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
