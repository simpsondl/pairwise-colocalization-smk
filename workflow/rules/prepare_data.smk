# First rule generates all pairs to test - now a checkpoint
checkpoint find_pairs:
    input:
        input_signals=config["DATA_DIR"] + "/input_signals.csv",
    output:
        output_pairs=config["OUTPUTS_DIR"] + "/pairs_to_test.csv",
        output_exact_matches=config["OUTPUTS_DIR"] + "/exact_variant_matches.csv",
    params:
        distance=config["DISTANCE_LIMIT"],
    log:
        config["OUTPUTS_DIR"] + "/logs/find_pairs.log",
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
    log:
        config["OUTPUTS_DIR"] + "/logs/extract_regions_{pair_id}.log",
    script:
        "../scripts/extract_regions.R"
