import pandas as pd
import os


# Read the pairs file to get list of pairs to process
if os.path.exists(config["OUTPUTS_DIR"] + "/pairs_to_test.csv"):
    PAIRS = pd.read_csv(config["OUTPUTS_DIR"] + "/pairs_to_test.csv")[
        "pair_index"
    ].tolist()
else:
    PAIRS = []
