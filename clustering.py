# import libraries
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import seaborn as sns
import argparse


def parse_args():
    parser = argparse.ArgumentParser()

    # Input and Output file names
    parser.add_argument('-i', '--input', nargs='?', default="example/out.tips", help="Input csv file name used for "
                                                                                     "clustering")
    parser.add_argument('-o', '--output', nargs='?', default="mapClustersSklearn.png", help="Output png file name")

    # Arguments for KMeans()
    parser.add_argument('-c', '--clusters', nargs='?', type=int, default=11, help="KMeans number of clusters argument")
    parser.add_argument('-n', '--n_init', nargs='?', type=int, default=40, help="KMeans n_init argument")
    parser.add_argument('-m', '--max', nargs='?', type=int, default=100, help="KMeans max_iter argument")
    parser.add_argument('-r', '--random', nargs='?', type=int, default=None, help="KMeans random_state argument")

    return parser.parse_args()


def make_map(args):
    # Create dataframe from CSV
    data = pd.read_csv(args.input)

    # Only use columns ag1 and ag2
    data = data[["ag1", "ag2"]]

    kmeans = KMeans(n_clusters=args.clusters, n_init=args.n_init, max_iter=args.max, random_state=args.random).fit(data)
    labels = kmeans.labels_

    data["clusters"] = labels

    data.groupby(["clusters"]).mean()

    # Create clustering graph
    sns.set_palette("Paired")
    sns.lmplot(x="ag1", y="ag2",
               data=data,
               fit_reg=False,
               hue="clusters",
               scatter_kws={"marker": "D",
                            "s": 7})
    plt.title("Map Clusters")
    plt.xlabel("ag1")
    plt.ylabel("ag2")
    plt.savefig(args.output, bbox_inches="tight")


def main():
    # Parse arguments
    arguments = parse_args()
    # Create graph
    make_map(arguments)


if __name__ == '__main__':
    main()
