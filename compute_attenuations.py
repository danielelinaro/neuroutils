# -*- coding: utf-8 -*-
"""
Copyright (C) 2024 Daniele Linaro <danielelinaro@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
"""

import os
import sys
import json
from tqdm import tqdm
import numpy as np

from dlutils.utils import extract_mechanisms
from neuroutils.trees import SWCImpedanceTree

progname = os.path.basename(sys.argv[0])


def find_param(params, name):
    D = {}
    section_IDs = {'somatic': 1, 'axonal': 2, 'basal': 3, 'apical': 4}
    for param in params:
        if param['param_name'] == name:
            if param['sectionlist'] == 'all':
                return param['value']
            D[section_IDs[param['sectionlist']]] = param['value']
    return D


def usage(exit_code=None):
    print(f'usage: {progname} [-h | --help] [-o | --outfile <filename>]')
    prefix = '       ' + ' ' * (len(progname)+1)
    print(prefix + '[-f | --force] [--no-plot] [-F | --frequency <value>] SWC_file')
    if exit_code is not None:
        sys.exit(exit_code)


if __name__ == '__main__':
    
    i = 1
    N_args = len(sys.argv)
    outfile = None
    force = False
    with_plot = False
    F = None

    while i < N_args:
        arg = sys.argv[i]
        if arg in ('-h', '--help'):
            usage(0)
        elif arg in ('-o','--outfile'):
            i += 1
            outfile = sys.argv[i]
        elif arg in ('-F', '--frequency'):
            i += 1
            F = float(sys.argv[i])
        elif arg in ('-f', '--force'):
            force = True
        elif arg == '--no-plot':
            with_plot = False
        elif arg[0] == '-':
            print(f'{progname}: unknown option `{arg}`.')
            sys.exit(1)
        else:
            break
        i += 1

    if i == N_args:
        print(f'{progname}: you must specify a configuration file')
        sys.exit(1)
    if i == N_args-1:
        config_file = sys.argv[i]
    else:
        print(f'{progname}: arguments after project name are not allowed')
        sys.exit(1)

    if not os.path.isfile(config_file):
        print(f'{progname}: {config_file}: no such file.')
        sys.exit(1)
    
    config = json.load(open(config_file, 'r'))
    if F is None:
        F = config['F']
    cm = config['cm']
    rm = config['rm']
    if isinstance(rm, dict):
        rm = {int(k): v for k,v in rm.items()}
    ra = config['ra']
    if isinstance(ra, dict):
        ra = {int(k): v for k,v in ra.items()}
    swc_file = config['swc_file']

    if outfile is None:
        outfile = '{}_F={:.4f}.npz'.\
            format(os.path.splitext(os.path.basename(config_file))[0], F)
    if os.path.isfile(outfile) and not force:
        print(f'{progname}: {outfile}: file exists, use -f to overwrite.')
        sys.exit(1)

    tree = SWCImpedanceTree(swc_file, cm, rm, ra, root_point=1)
    tree.compute_impedances(F=F)
    tree.compute_attenuations()

    dtypes = [('ID',int), ('type',int),
              ('x',float), ('y',float), ('z',float), ('diam',float),
              ('distance_from_soma',float),
              ('A_out', float), ('A_in', float)]
    data = np.array([(node.ID, node.type,
                      node.x, node.y, node.z, node.diam,
                      node.distance,
                      np.abs(node.A_from_root), 1.0) for node in tree],
                    dtype=dtypes)
    
    dict_IDs = list(tree._nodes_dict.keys())
    assert len(data) == len(dict_IDs)
    assert all([a==b for a,b in zip(data['ID'], dict_IDs)])
    N_nodes = len(data)
    
    for i in tqdm(range(1, N_nodes), ascii=True, ncols=80):
        node = tree._nodes_dict[data['ID'][i]]
        SWCImpedanceTree.make_root(node)
        # !!! DON'T FORGET THE FOLLOWING LINE !!!
        # it is necessary to compute all tree branches
        tree.root = node
        tree.compute_impedances(F=F)
        tree.compute_attenuations()
        data['A_in'][i] = tree.compute_attenuation(node.ID, 1)
    
    np.savez_compressed(outfile, data=data, F=F, cm=cm, ra=ra, rm=rm, config=config)
    
    if with_plot:
        import matplotlib.pyplot as plt
        outfile = os.path.splitext(outfile) + '.pdf'
        
    
