""" A script to edit parameters from a yaml file to enable simulation parameter
sweeps.
"""
import os
import yaml
import argparse
from itertools import product

def parse_yaml(file_name):
    """ Store a .yml file as a dictionary.
    """
    with open(file_name) as params_file:
        params = yaml.load(params_file, Loader=yaml.FullLoader)
    return params

def dump_yaml(data, out_path):
    """ Dump a .yml file.
    """
    with open(out_path, 'w') as file:
        yaml.dump(data, file)
        print(f'Updates simulation parameters dumped to {out_path}.')

def product_dict(options_dict):
    """ Take in a dictionary of changes to make, and create a dictionary of the cartesian product.

    Parameters
    -----------
    options_dict : Dictionary
        A dictionary of desired hyperparameters to try in antigen sims.

    Returns
    --------
    dict_list : list[Dictionary, ...]
        A list of dictionaries that are analgous to one yaml configuration file.
    """
    keys = []
    vals = []
    for key, val in options_dict.items():
        if isinstance(val, list):
            vals.append(val)
        else:
            vals.append([val])
        keys.append(key)
    return [dict(zip(keys, items)) for items in product(*vals)]

def make_dirs(new_params_list, runs=1):
    """ Create the directories needed for output files.

    Take in a list of dictionaries with key-value pairs of parameter configs.
    For each entry in the list, create a directory in key_val_key_val_... format.
    If runs is greater than 1 --, in each directory, create sub-directories for each run.
    """
    path_list = []
    # Loop through list of dictinoaries
    for params in new_params_list:
        config_path = "" # Stem for new hyperparam configuration path.
        # Iterater over keys and values:
        for key, val in params.items():
            config_path += str(key) + '_' + str(val) + '_'
        config_path = config_path[:-1]
        path_list.append(config_path)
        if not os.path.exists(config_path):
            os.makedirs(config_path)
    return path_list

def create_yaml_files(original_params, path_list, params_list):
    """ Creates the relevant yaml files in each directory.
    """
    for path, param_set in zip(path_list, params_list):
        # Create a file called parameters.yml in each subdir.
        yml_path = path + "/parameters.yml"
        original_params.update(param_set)
        dump_yaml(original_params, yml_path)

def main():
    """ Run the main program. """
    parser = argparse.ArgumentParser(description='Update parameter file.')
    parser.add_argument('in_file', help='Path to source parameter.yml file.', type=str)
    parser.add_argument('edits_file', help='Path to YAML file holding desired parameter updates.', type=str)

    args = parser.parse_args()

    params = parse_yaml(args.in_file)
    new_params = parse_yaml(args.edits_file)
    new_params_list = list(product_dict(new_params))
    path_list = make_dirs(new_params_list)
    create_yaml_files(params, path_list, new_params_list)

if __name__ == "__main__":
    main()
