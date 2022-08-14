# import libraries
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter

def make_chart():
    data = pd.read_csv("mutations.csv")
    starting_sequence = str(data.iloc[0])

    nucleotide_counter = Counter(starting_sequence)
    print('Starting sequence nucleotide count:')
    print('A:', nucleotide_counter['A'])
    print('C:', nucleotide_counter['C'])
    print('G:', nucleotide_counter['G'])
    print('T:', nucleotide_counter['T'])

    mutations = data.iloc[1:]
    plt.figure()

    mutation_counts = mutations.value_counts()
    transition_counts = mutation_counts[['AG', 'GA', 'CT', 'TC']].sum()
    transversion_counts = mutation_counts[['AC', 'AT', 'GC', 'GT', 'TA', 'TG', 'CA', 'CG']].sum() \

    print('Number of transitions:', transition_counts)
    print('Number of transversion:', transversion_counts)

    ratioString = 'Calculated transition transversion ratio {transition:.3f} : {transversion:.1f}'.format(transition=transition_counts/transversion_counts, transversion=1)

    ax = mutation_counts.plot.bar()
    fig = ax.get_figure()
    plt.suptitle('Nucleotide Mutations')
    plt.title(ratioString)
    plt.xlabel('Wild type to mutant nucleotide mutation')
    plt.ylabel('Count')
    fig.savefig('mutation_distribution.png', bbox_inches = "tight")


def main():
    # Create bar chart
    make_chart()


if __name__ == '__main__':
    main()