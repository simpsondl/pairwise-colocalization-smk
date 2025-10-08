rule generate_analyses:
    input:
        input_signals=config["DATA_DIR"] + "/input_signals.csv",
    output:
        output_pairs=config["OUTPUTS_DIR"] + "/pairwise_colocalizations_to_check.csv",
    script:
        "../scripts/generate_analyses.R"


# rule create_joined_sumstats:
#     input:
#         input_sumstats_appended="preprocessed.sorted.tsv",
#         input_pairs=rules.generate_analyses.output.output_pairs,
#     output:
#         output_joined_sumstats=config["OUTPUTS_DIR"]
#         + "/joined_sumstats_{pairwise_analysis}.tsv",
