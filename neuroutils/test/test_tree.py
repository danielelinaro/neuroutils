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


import unittest
from ..nodes import Node
from ..trees import Tree

def make_binary_tree(node, num_levels):
    num_levels -= 1
    if num_levels == 0:
        return
    node.children = [Node(ID=node.ID+suffix) for suffix in 'LR']
    for child in node.children:
        make_binary_tree(child, num_levels)

class TestTree(unittest.TestCase):

    def test_find_node(self):
        root = Node(ID='')
        num_levels = 4
        make_binary_tree(root, num_levels)
        tree = Tree(root)
        node = tree.find_node(Node(ID='LR'))
        self.assertIsNotNone(node, "Node not in tree")

    def test_find_path(self):
        root = Node(ID='')
        num_levels = 4
        make_binary_tree(root, num_levels)
        tree = Tree(root)
        ID_i,ID_j = '','RLR'
        path = tree.find_connecting_path(ID_i, ID_j)
        self.assertIsNotNone(path, f"No path between nodes with IDs '{ID_i}' and '{ID_j}'.")

