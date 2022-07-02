# import libraries
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import scipy.stats as stats

def make_map():
    for i in range(10): # update range() to change the number of graphs to produce
        data = pd.read_csv("valuesGammaDistribution/0_site" + str(i) + ".csv")

        data = data[["r", "theta"]]

        plt.figure()
        plt.title('Density plot = r')
        sns.distplot(data['r'])

        # ~Gamma
        x = np.arange(0, 5, 0.1)
        y = stats.gamma.pdf(x, a=1, scale=0.3)

        plt.plot(x, y)
        plt.show()


def main():
    # Create graph
    make_map()


if __name__ == '__main__':
    main()