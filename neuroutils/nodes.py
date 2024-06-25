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


__all__ = ['Node','ImpedanceNode']

import numpy as np

class Node (object):
    def __init__(self, ID, parent=None, children=None):
        self.ID = ID
        self.parent = parent
        self.children = children if children is not None else []

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


class ImpedanceNode (Node):
    def __init__(self, seg, parent=None, children=None):
        self.seg = seg
        sec,x = seg.sec, seg.x
        # the section that contains this segment
        self._sec = sec
        # the position in the containing section
        self._x = x
        super().__init__(ID=ImpedanceNode.make_ID(seg),
                         parent=parent,
                         children=children)
        # the length of this segment
        self._L    = sec.L / sec.nseg * 1e-4  # [cm]
        self._diam = seg.diam * 1e-4          # [cm]
        # membrane resistance per unit area
        self._rm   = 1/seg.g_pas              # [Ohm cm2]
        # membrane capacitance per unit area
        self._cm = seg.cm*1e-6                # [F/cm2]
        # membrane time constant
        self._taum = self._rm * self._cm      # [s]
        # axial resistivity
        self._ra   = sec.Ra                   # [Ohm cm]
        # half the total axial resistance
        self._Ra   = 2*self._ra*self._L/(np.pi*self._diam**2) # [Ohm]
        self._A = np.zeros(len(self.children), dtype=complex)

    def make_ID(seg):
        return '{}-{:.4f}'.format(seg.sec.name(), seg.x)

    @property
    def Ra(self):
        return self._Ra
    @property
    def Zm(self):
        return self._Zm
    @property
    def Zp(self):
        return self._Zp
    @property
    def Zload(self):
        return self._Zload
    @property
    def A(self):
        return self._A
    
    def compute_impedances(self, F):
        self._Zm = self._rm / (np.pi*self._diam*self._L*(1+1j*2*np.pi*F*self._taum))
        n_children = len(self.children)
        if n_children == 0:
            self._Zp = self._Zm
        else:
            Z = np.zeros(n_children, dtype=complex)
            for i,child in enumerate(self.children):
                child.compute_impedances(F)
                Z[i] = child.Zload
            # the parallel of all children impedances
            Z = 1/np.sum(1/Z)
            self._Zp = (Z+self.Ra)*self.Zm/(Z+self.Ra+self.Zm)
        self._Zload = self.Ra + self.Zp
        
    def compute_attenuations(self):
        n_children = len(self.children)
        self._A = np.zeros(n_children, dtype=complex)
        for i,child in enumerate(self.children):
            self._A[i] = 1 + (self.Ra + child.Ra) / child.Zp
            child.compute_attenuations()

    def __eq__(self, other):
        return other.seg == self.seg

