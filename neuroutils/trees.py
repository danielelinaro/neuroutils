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


from .nodes import Node, ImpedanceNode
import numpy as np

__all__ = ['Tree', 'ImpedanceTree']


class Tree (object):
    def __init__(self, root):
        self.root = root
        
    @property
    def root(self):
        return self._root
    @root.setter
    def root(self, r):
        self._root = r

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

    def find_connecting_path(self, ID_i, ID_j):
        node_i = self.find_node_with_ID(ID_i)
        if node_i is None:
            raise ValueError(f"ID '{ID_i}' not in tree")
        node_j = self.find_node_with_ID(ID_j)
        if node_j is None:
            raise ValueError(f"ID '{ID_j}' not in tree")
        path = [node_j]
        node = node_j
        while node.parent is not None:
            idx = node.parent.children.index(node)
            node = node.parent
            path.append(node)
            if node == node_i:
                break
        if path[-1] != node_i:
            return None
        return path[::-1]


class ImpedanceTree (Tree):
    def __init__(self, root_node=None, root_sec=None, root_seg=None):
        if root_node is not None:
            root = self._make_branch(root_node.seg.sec)
        elif root_sec is not None:
            root = self._make_branch(root_sec)
        elif root_seg is not None:
            root = self._make_branch(root_seg.sec)
        else:
            raise Exception('one of root_node, root_sec or root_seg must be passed')
        super().__init__(root)

    def compute_impedances(self, F):
        self.root.compute_impedances(F)
        
    def compute_attenuations(self):
        self.root.compute_attenuations()

    def _make_branch(self, sec):
        branch = [ImpedanceNode(seg) for seg in sec]
        n_nodes = len(branch)
        for i in range(n_nodes-1):
            branch[i].add_child(branch[i+1])
        for child in sec.children():
            child_node = self._make_branch(child)
            branch[-1].add_child(child_node)
        return branch[0]

    def compute_attenuation(self, seg_i, seg_j, full_output=False):
        path = super().find_connecting_path(ImpedanceNode.make_ID(seg_i),
                                            ImpedanceNode.make_ID(seg_j))
        n_nodes = len(path)-1
        A_on_path = np.zeros(n_nodes, dtype=complex)
        for i in range(n_nodes):
            j = path[i].children.index(path[i+1])
            A_on_path[i] = path[i].A[j]
        A = np.abs(np.prod(A_on_path))
        if full_output:
            return A,A_on_path
        return A

