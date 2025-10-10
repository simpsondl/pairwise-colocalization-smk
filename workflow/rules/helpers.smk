import pandas as pd
import os


# Read the pairs file to get list of pairs to process
if os.path.exists(config["OUTPUTS_DIR"] + "/pairs_to_test.csv"):
    PAIRS = pd.read_csv(config["OUTPUTS_DIR"] + "/pairs_to_test.csv")[
        "pair_index"
    ].tolist()
else:
    PAIRS = []

# Get list of credible sets to join
if os.path.exists(config["OUTPUTS_DIR"] + "/colocalized_signal_codes.txt"):
    COLOCALIZED1 = pd.read_csv(config["OUTPUTS_DIR"] + "/colocalized_signal_codes.txt")[
        "analysis1"
    ].tolist()
    COLOCALIZED2 = pd.read_csv(config["OUTPUTS_DIR"] + "/colocalized_signal_codes.txt")[
        "analysis2"
    ].tolist()
    COLOCALIZED = list(set(COLOCALIZED1 + COLOCALIZED2))
    COLOC_IDS = pd.read_csv(config["OUTPUTS_DIR"] + "/colocalized_signal_codes.txt")[
        "coloc_id"
    ].tolist()
else:
    COLOCALIZED = []
