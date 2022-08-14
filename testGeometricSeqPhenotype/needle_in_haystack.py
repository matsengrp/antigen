# unfinished
import pandas as pd
import numpy as np

def parse():
    max_e = float('-inf')
    max_e_virus = ""

    max_ne = float('-inf')
    max_ne_virus = ""

    visited = set()
    trunk_Mutations = np.array([0,0])
    branch_Mutations = np.array([0,0])
    out_f = 'out.branches'
    with open(out_f) as f:
        lines = f.readlines()
    counts_dict = {}
    for line in lines:
        elements = line.split(',')
        #print(elements)
        child = elements[0].split('"')[1]
        v_isTrunk = elements[2]
        v_eMutations = int(elements[10].strip())
        v_neMutations = int(elements[11].split('}')[0].strip())
        parent = elements[11].split('"')[1]
        vp_eMutations = int(elements[21].strip())
        vp_neMutations = int(elements[22].split('}')[0].strip())

        vp_isTrunk = elements[13]
        key = (vp_isTrunk, v_isTrunk)
        if key not in counts_dict:
            counts_dict[key] = 1
        else:
            counts_dict[key] += 1

        v_Mutations = np.array([v_eMutations, v_neMutations])
        vp_Mutations = np.array([vp_eMutations, vp_neMutations])
        diff = v_Mutations - vp_Mutations
        #print(vp_Mutations, v_Mutations)
        if key == ('1', '1'):
            if v_eMutations > max_e:
                max_e = v_eMutations
                max_e_virus = child
            if v_neMutations > max_ne:
                max_ne = v_neMutations
                max_ne_virus = child

            trunk_Mutations += diff
        else:
            if sum(diff) == 0:
                print(child, parent)
            branch_Mutations += diff
        visited.add(parent)
    print(counts_dict)
    print(trunk_Mutations)
    print(branch_Mutations)

def main():

    # Goal: count the number of epitope vs nonepitope mutations
    # summing them up based on whether they occur between two
    # trunk nodes or not

    # Initialize at zero, zero
    # v_Mutations = np.array([6,10])
    # vp_Mutations = np.array([5,10])
    # diff = v_Mutations - vp_Mutations
    # trunk_Mutations += diff
    # print(diff)
    # print(trunk_Mutations)

    parse()


if __name__ == '__main__':
    main()