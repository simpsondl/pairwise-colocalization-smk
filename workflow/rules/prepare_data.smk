# First rule generates all pairs to test
rule find_pairs:
    input:
        input_signals=config["DATA_DIR"] + "/input_signals.csv",
    output:
        output_pairs=config["OUTPUTS_DIR"] + "/pairs_to_test.csv",
    params:
        distance=config["DISTANCE_LIMIT"],
    script:
        "../scripts/generate_analyses.R"


# Extract regions for a single pair
rule extract_regions:
    input:
        pairs=config["OUTPUTS_DIR"] + "/pairs_to_test.csv",
        sumstats_dir=config["SUMSTATS_DIR"],
    output:
        output_joined_sumstats=config["OUTPUTS_DIR"]
        + "/extracted_regions/{pair_id}_sumstats.txt",
    params:
        pair_id="{pair_id}",
    script:
        "../scripts/extract_regions.R"
