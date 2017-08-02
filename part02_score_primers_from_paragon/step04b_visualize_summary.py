import pandas as pd
from bokeh.charts import Scatter, Histogram, output_file, save
from bokeh.layouts import column

def visualize_summary():
    output_file("summary.html")

    df = pd.read_csv("summary.csv")

    p = Scatter(df, x='insert_length', y='score', color="max_ambiguities(l|r)",
                title="Score vs Insert_Length Colored by Max_Ambiguities(L|R)")

    hl = Histogram(df, values="insert_length", title="Insert_Length Distribution Colored by Max_Ambiguities(L|R)",
                   color="max_ambiguities(l|r)")

    hs = Histogram(df, values='score', title="Score Distribution Colored by Max_Ambiguities(L|R)",
                   color="max_ambiguities(l|r)", bins=7)

    save(column(p, hl, hs))

if __name__ == '__main__':
    import argparse
    import pytoml

    parser = argparse.ArgumentParser()
    parser.add_argument('--configPath', help='configPath', default='../config.toml')
    args = parser.parse_args()

    # load config file
    with open(args.configPath) as toml_data:
        config = pytoml.load(toml_data)
