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


__all__ = ['Node','Tree']


class Node (object):
    def __init__(self, ID, parent=None, children=[]):
        self.ID = ID
        self.parent = parent
        self.children = children

    @property
    def ID(self):
        return self._ID
    @ID.setter
    def ID(self, value):
        self._ID = value

    @property
    def parent(self):
        return self._parent
    @parent.setter
    def parent(self, p):
        self._parent = p

    @property
    def children(self):
        return self._children
    @children.setter
    def children(self, children):
        if not isinstance(children, list):
            raise ValueError('children must be a list')
        if not all([isinstance(child, Node) for child in children]):
            raise ValueError('all children must be nodes')
        self._children = children
        for child in self.children:
            child.parent = self
    def add_child(self, child):
        if not isinstance(child, Node):
            raise ValueError('child must be a node')
        self.children.append(child)
        child.parent = self

    def __str__(self):
        return f"'{self.ID}'"

    def __eq__(self, other):
        return other.ID == self.ID


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

