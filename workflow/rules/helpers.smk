import pandas as pd
import os


# Input functions that read checkpoint outputs at runtime
def get_pairs_to_test(wildcards):
    """Get list of pair IDs from the checkpoint output"""
    checkpoint_output = checkpoints.find_pairs.get(**wildcards).output.output_pairs
    pairs_df = pd.read_csv(checkpoint_output)
    return pairs_df["pair_index"].tolist()


def get_coloc_ids(wildcards):
    """Get list of colocalized signal IDs"""
    # This depends on the define_colocalized_signals rule completing
    coloc_file = config["OUTPUTS_DIR"] + "/colocalized_signal_codes.txt"
    if os.path.exists(coloc_file):
        coloc_df = pd.read_csv(coloc_file)
        return coloc_df["coloc_id"].tolist()
    else:
        return []


def get_coloc_ids_from_checkpoint(wildcards):
    """Get colocalized signal IDs from checkpoint output"""
    checkpoint_output = checkpoints.define_colocalized_signals.get(
        **wildcards
    ).output.output_coloc_codes
    coloc_df = pd.read_csv(checkpoint_output)
    return coloc_df["coloc_id"].tolist()
