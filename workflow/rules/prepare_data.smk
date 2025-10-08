rule generate_analyses:
    input:
        input_signals="{DATA_DIR}/input_signals.csv"
    output:
        output_pairs="{OUTPUTS_DIR}/pairwise_colocalization_to_perform.csv"
    script:
        "{SCRIPTS_DIR}/generate_analyses.R"
