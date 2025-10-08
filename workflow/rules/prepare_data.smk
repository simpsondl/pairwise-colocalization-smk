rule generate_analyses:
    input:
        input_signals=config["DATA_DIR"] + "/input_signals.csv",
    output:
        output_pairs=config["OUTPUTS_DIR"] + "/pairwise_colocalization_to_perform.csv",
    script:
        "../scripts/generate_analyses.R"
