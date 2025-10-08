import pandas as pd

# SAMPLE_DATA = {}


# def load_sample_data(wildcards):
#     """
#     Reads the CSV output by the checkpoint into a global dictionary
#     and returns the list of all sample IDs (for use in rule all).
#     """
#     global SAMPLE_DATA

#     # 1. Wait for the checkpoint to finish
#     checkpoint_output = checkpoints.generate_data_list.get(**wildcards)

#     # The actual input CSV is at: checkpoint_output.input[0]
#     csv_file = checkpoint_output.input[0]

#     # 2. Read the CSV using Pandas
#     # Use 'ID' as the index (this is the key for easy lookup later)
#     df = pd.read_csv(csv_file).set_index("ID", drop=False)

#     # 3. Store the data in the global dictionary for later rule access
#     # We convert it to dict for fast, easy Snakemake access
#     SAMPLE_DATA = df.T.to_dict()

#     # 4. Return the list of IDs (from the index) to dynamically expand rules
#     return list(df.index)


# def input_for_data_check(wildcards):
#     # This function is run during Snakemake's DAG construction.
#     # It forces the execution of the full data loading logic (which populates SAMPLE_DATA).

#     # We call load_sample_data() but discard the list of IDs it returns,
#     # because we only need the side effect of populating SAMPLE_DATA.
#     # We MUST return the file path that is the dependency for the checkpoint.

#     load_sample_data(wildcards)  # This populates the global dictionary
#     # The rule needs a file to depend on, which is the checkpoint's input file
#     return "params/.data_ready.txt"
