# import libraries
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

def make_chart():
    data = pd.read_csv("mutation.csv")
    starting_sequence = str(data.iloc[0])

    nucleotide_counter = Counter(starting_sequence)
    print('A:', nucleotide_counter['A'])
    print('C:', nucleotide_counter['C'])
    print('G:', nucleotide_counter['G'])
    print('T:', nucleotide_counter['T'])

    mutations = data.iloc[1:]
    plt.figure()

    ax = mutations.value_counts().plot.bar()
    fig = ax.get_figure()
    fig.savefig('mutation_distribution.png')


def main():
    # Create bar chart
    make_chart()


if __name__ == '__main__':
    main()