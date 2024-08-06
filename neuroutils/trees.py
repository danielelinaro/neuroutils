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


from .nodes import Node, ImpedanceNode, SWCImpedanceNode, SWCNode
import numpy as np
from neuron import h

__all__ = ['Tree', 'BaseImpedanceTree', 'ImpedanceTree', 'SWCImpedanceTree', 'SWCTree']


class Tree (object):
    def __init__(self, root):
        self.root = root
        self._branches = []
        self._make_branches(self.root, self.branches)

    @property
    def root(self):
        return self._root
    @root.setter
    def root(self, r):
        self._root = r

    @property
    def branches(self):
        return self._branches
        
    def _gather_nodes(self, node, node_list):
        if not node is None:
            node_list.append(node)
            for child in node.children:
                self._gather_nodes(child, node_list)

    def __iter__(self):
        nodes = []
        self._gather_nodes(self.root, nodes)
        for n in nodes:
            yield n
    
    def find_node(self, node):
        for n in self:
            if n == node:
                return n
        return None

    def find_node_with_ID(self, ID):
        for n in self:
            if n.ID == ID:
                return n
        return None

    def _make_branches(self, node, branches):
        branch = []
        while len(node.children) == 1:
            branch.append(node)
            node = node.children[0]
        branch.append(node)
        branches.append(branch)
        for child in node.children:
            self._make_branches(child, branches)

    def find_connecting_path(self, ID_i, ID_j):
        node_i = self.find_node_with_ID(ID_i)
        if node_i is None:
            raise ValueError(f"ID '{ID_i}' not in tree")
        node_j = self.find_node_with_ID(ID_j)
        if node_j is None:
            raise ValueError(f"ID '{ID_j}' not in tree")
        def _find_path_to_root(node):
            path = [node]
            while node.parent is not None:
                idx = node.parent.children.index(node)
                node = node.parent
                path.append(node)
            return path[::-1]
        path_i = _find_path_to_root(node_i)
        try:
            idx = path_i.index(node_j)
            return path_i[idx:]
        except:
            pass
        path_j = _find_path_to_root(node_j)
        try:
            idx = path_j.index(node_i)
            return path_j[idx:]
        except:
            pass
        depth = {}
        for node in path_i:
            if len(node.children) > 1 and node in path_j:
                depth[path_j.index(node)] = node
        max_depth = max(list(depth.keys()))
        # the deepest common ancestor of node_i and node_j
        ancestor = depth[max_depth]
        path = path_i[path_i.index(ancestor)+1:][::-1] + path_j[path_j.index(ancestor):]
        return path


class BaseImpedanceTree (Tree):
    def __init__(self, root):
        super().__init__(root)

    def compute_impedances(self, F):
        self.root.compute_impedances(F)
        
    def compute_attenuations(self):
        self.root.compute_attenuations()

    def compute_attenuation(self, ID_start, ID_end, full_output=False):
        path = super().find_connecting_path(ID_start, ID_end)
        n_nodes = len(path)-1
        A_on_path = np.zeros(n_nodes, dtype=complex)
        for i in range(n_nodes):
            j = path[i].children.index(path[i+1])
            A_on_path[i] = path[i].A[j]
        A = np.abs(np.prod(A_on_path))
        if full_output:
            return A,A_on_path
        return A


_get_ith_segment = lambda sec,i: [seg for seg in sec][i]
_get_first_segment = lambda sec: _get_ith_segment(sec,0)
_get_last_segment = lambda sec: _get_ith_segment(sec,-1)

class ImpedanceTree (BaseImpedanceTree):
    def __init__(self, root_seg=None, root_node=None):
        if root_seg is not None:
            root = ImpedanceNode(root_seg)
        elif root_node is not None:
            root = root_node
        else:
            raise Exception('one of root_node or root_seg must be passed')
        super().__init__(root)
        self._make_branch_increasing_x(root)
        self._make_branch_decreasing_x(root)

    def _make_branch_increasing_x(self, node):
        branch = [ImpedanceNode(seg) for seg in node.sec if seg.x > node.seg.x]
        if len(branch) > 0:
            node.add_child(branch[0])
            for node_i,node_j in zip(branch[:-1],branch[1:]):
                node_i.add_child(node_j)
            parent_node = branch[-1]
        else:
            parent_node = node
        for child_sec in node.sec.children():
            child_node = ImpedanceNode(_get_first_segment(child_sec),
                                       parent=parent_node)
            self._make_branch_increasing_x(child_node)

    def _make_branch_decreasing_x(self, node):
        branch = [ImpedanceNode(seg) for seg in node.sec if seg.x < node.seg.x][::-1]
        if len(branch) > 0:
            node.add_child(branch[0])
            for node_i,node_j in zip(branch[:-1],branch[1:]):
                node_i.add_child(node_j)
            parent_node = branch[-1]
        else:
            parent_node = node
        sref = h.SectionRef(node.sec)
        if sref.has_parent():
            parent_sec = sref.parent
            parent_node = ImpedanceNode(_get_last_segment(parent_sec),
                                        parent=parent_node)
            for child_sec in parent_sec.children():
                if child_sec != node.sec:
                    child_node = ImpedanceNode(_get_first_segment(child_sec),
                                               parent=parent_node)
                    self._make_branch_increasing_x(child_node)
            self._make_branch_decreasing_x(parent_node)

    def compute_attenuation(self, end_seg, full_output=False):
        return super().compute_attenuation(ImpedanceNode.make_ID(self.root.seg),
                                           ImpedanceNode.make_ID(end_seg),
                                           full_output)


class SWCImpedanceTree (BaseImpedanceTree):
    def __init__(self, swc_file, cm, rm, ra):
        def make_dict(value):
            if isinstance(value, dict):
                return value
            if isinstance(value, (int,float)):
                return {k: value for k in (1,2,3,4)}
            raise ValueError('Argument must be a dictionary or a number')
        cm_dict = make_dict(cm)
        rm_dict = make_dict(rm)
        ra_dict = make_dict(ra)
        import pandas as pd
        from collections import OrderedDict
        col_names = 'ID','typ','x','y','z','diam','parent_ID'
        col_types = {'ID': np.int32, 'typ': np.int32, 'x': np.float32,
                     'y': np.float32, 'z': np.float32, 'diam': np.float32,
                     'parent_ID': np.int32}
        df = pd.read_table(swc_file, sep=' ', header=None, names=col_names, index_col='ID')
        nodes = OrderedDict()
        soma_idx, = np.where(df.loc[:,'typ'] == 1)
        if len(soma_idx) == 3:
            # 3-point soma
            coords = df.loc[[2,3],['x','y','z']].to_numpy()
            diams = df.loc[[2,3],'diam'].to_numpy()
            typ = 1
            nodes[1] = SWCImpedanceNode(1, coords, diams, cm_dict[typ],
                                        rm_dict[typ], ra_dict[typ],
                                        node_type=typ)
        else:
            raise NotImplementedError('Only morphologies with 3-point soma are supported')
        for ID,child in df.iterrows():
            if child.typ != 1:
                parent = df.loc[child.parent_ID]
                coords = [[child.x, child.y, child.z],[parent.x, parent.y, parent.z]]
                diams = [child.diam, parent.diam]
                nodes[ID] = SWCImpedanceNode(ID, coords, diams,
                                             cm_dict[child.typ], rm_dict[child.typ],
                                             ra_dict[child.typ], child.typ,
                                             parent=nodes[child.parent_ID])
        self._nodes_dict = nodes
        super().__init__(nodes[1])
        
    @property
    def nodes(self):
        return self._nodes_df



class SWCTree (Tree):
    def __init__(self, swc_file):
        from collections import OrderedDict
        data = np.loadtxt(swc_file)
        idx = (data[:,1] == 2) | (data[:,1] == 3) | (data[:,1] == 4)
        x = data[idx, 2]
        y = data[idx, 3]
        z = data[idx, 4]
        self.xy_ratio = (x.max() - x.min()) / (y.max() - y.min())
        self.bounds = np.array([[x.min(), x.max()], [y.min(), y.max()], [z.min(), z.max()]])
        nodes = OrderedDict()
        for row in data:
            node_id   = int(row[0])
            node_type = int(row[1])
            x, y, z,  = row[2:5]
            diam      = row[5]
            parent_id = int(row[6])
            parent = nodes[parent_id] if parent_id > 0 else None
            nodes[node_id] = SWCNode(x,y,z,diam,node_type,node_id,parent)
        self._nodes_dict = nodes
        keys = list(nodes.keys())
        super().__init__(nodes[keys[0]])

    @property
    def nodes(self):
        return self._nodes_dict

    def plot(self, type_ids=(1,2,3,4), scalebar_length=None, cmap=None, points=None, values=None,
             cbar_levels=None, cbar_ticks=10, cbar_orientation='vertical', cbar_label='', ax=None,
             bounds=None, diam_coeff=1, cbar_ticks_fun=lambda x: x):
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        from matplotlib.collections import LineCollection
        from matplotlib import cm
    
        if ax is None:
            ax = plt.gca()

        if points is None or values is None:
            uniform_color_branches = True
            if cmap is None:
                color_fun = lambda i: [
                    [0,0,0],    # soma
                    [.2,.2,.2], # axon
                    [.7,0,.7],  # basal
                    [0,.7,0]    # apical
                ][i-1]
            elif isinstance(cmap, dict):
                color_fun = lambda key: cmap[key]
            else:
                color_fun = cmap
        else:
            from scipy.interpolate import NearestNDInterpolator
            uniform_color_branches = False
            interp = NearestNDInterpolator(points, values)
            norm = colors.Normalize(vmin = values.min(), vmax = values.max())

        for branch in self.branches:
            if branch[0].type not in type_ids:
                continue
            if branch[0].parent is not None:
                node = branch[0].parent
                xyzd = np.concatenate((np.array([node.x, node.y, node.z, node.diam * diam_coeff], ndmin=2),
                                       [[node.x, node.y, node.z, node.diam * diam_coeff] for node in branch]))
            else:
                xyzd = np.array([[node.x, node.y, node.z, node.diam * diam_coeff] for node in branch])

            if branch[0].parent is not None and branch[0].parent.type == 1:
                xyzd[0,-1] = xyzd[1,-1]

            xy = xyzd[:,:2].reshape(-1, 1, 2)
            segments = np.concatenate([xy[:-1], xy[1:]], axis=1)
            if uniform_color_branches:
                lc = LineCollection(segments, linewidths=xyzd[:,-1]/2, colors=color_fun(branch[0].type))
            else:
                lc = LineCollection(segments, linewidths=xyzd[:,-1]/2, cmap=cmap, norm=norm)
                lc.set_array(interp(xyzd[:,:3]))
            line = ax.add_collection(lc)

            if bounds is not None:
                ax.set_xlim(bounds[0])
                ax.set_ylim(bounds[1])
            else:
                ax.set_xlim(self.bounds[0])
                ax.set_ylim(self.bounds[1])
                ax.axis('equal')

        if scalebar_length is not None:
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            x = xlim[0] / 1.5
            y = (ylim[1] - scalebar_length) / 2
            ax.plot(x + np.zeros(2), y + np.array([0, scalebar_length]), 'k', lw=2)
            ax.text(x - np.diff(xlim)/15, y + scalebar_length/2, r'{} $\mu$m'.format(scalebar_length), fontsize=12, \
                    horizontalalignment='center', verticalalignment='center', rotation=90)

        if not uniform_color_branches and cbar_levels is not None:
            if np.isscalar(cbar_ticks):
                ticks = np.linspace(values.min(), values.max(), cbar_ticks)
            else:
                ticks = cbar_ticks
            cbar = plt.colorbar(line, ax=ax, fraction=0.1, shrink=0.5, aspect=30, ticks=ticks, orientation=cbar_orientation)
            cbar.ax.set_yticklabels(cbar_ticks_fun(ticks))
            if cbar_orientation == 'vertical':
                cbar.ax.set_ylabel(cbar_label)
            else:
                cbar.ax.set_xlabel(cbar_label)

