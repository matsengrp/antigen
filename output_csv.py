"""
Output case counts per day and deme
using (csv of) out.timeseries from antigen
"""
import pandas as pd
import yaml
import argparse


COLUMN_NAMES = ['date', 'location', 'cases']


def parse_args():
    parser = argparse.ArgumentParser()

    # Input and Output file names
    parser.add_argument('-i', '--input', nargs='?', default="example/out_timeseries.csv",
                        help="Input csv file name used for outputting cases per day, deme, and variant ")
    parser.add_argument('-o', '--output', nargs='?', default="demes.timeseries", help="Output timeseries file name")

    return parser.parse_args()


def count_demes(args):
    # Read CSV
    data = pd.read_csv(args.input)

    # Create DataFrame
    df = pd.DataFrame(columns=COLUMN_NAMES)

    print(df)

    # Parse parameters.yml
    with open('parameters.yml') as f:
        parameters_data = yaml.load(f, Loader=yaml.FullLoader)

    # Get demeNames from parameters.yml
    deme_names = parameters_data.get('demeNames')

    # Iterate over each deme
    for deme in deme_names:
        deme_index = deme_names.index(deme)
        # Iterate over each row in out.timeseries
        for index, row in data.iterrows():
            # Make dataframe for current deme in current row
            next = [[row[0], deme, row[20 + 10 * deme_index]]]
            df_next = pd.DataFrame(next, columns=COLUMN_NAMES)
            df = pd.concat([df, df_next], ignore_index=True)

    # Output file
    df.to_csv(args.output, index=False, sep='\t')


def main():
    # Parse arguments
    arguments = parse_args()

    # Create output file per demes
    count_demes(arguments)


if __name__ == '__main__':
    main()