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


import numpy as np

__all__ = ['Node', 'BaseImpedanceNode', 'ImpedanceNode', 'SWCImpedanceNode', 'SWCNode']


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
        if hasattr(self,'_parent'):
            if self._parent == p:
                return
            if self._parent is not None:
                self._parent.remove_child(self)
        self._parent = p
        if p is not None and self not in p.children:
            p.add_child(self)

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
    def remove_child(self, node):
        idx = self.children.index(node)
        self.children.pop(idx)

    def __str__(self):
        return f"'{self.ID}'"

    def __eq__(self, other):
        return other is not None and other.ID == self.ID


class BaseImpedanceNode (Node):
    def __init__(self, ID, L, diam, cm, rm, ra, parent=None, children=None):
        super().__init__(ID=ID, parent=parent, children=children)
        self._L    = L
        self._diam = diam
        # membrane capacitance per unit area
        self._cm   = cm
        # membrane resistance per unit area
        self._rm   = rm
        # axial resistivity
        self._ra   = ra
        # membrane time constant
        self._taum = self._rm * self._cm      # [s]
        # DC length constant
        self._lambda_DC = 0.5*np.sqrt(self._diam*self._rm/self._ra) # [cm]
        # DC input resistance
        self._r_inf = 2*np.sqrt(self._ra*self._rm) / \
                      (np.pi*self._diam**(3/2)) # [Ohm]
        self._A = np.zeros(len(self.children), dtype=complex)

    @property
    def L(self):
        return self._L
    @property
    def diam(self):
        return self._diam
    @property
    def cm(self):
        return self._cm
    @property
    def rm(self):
        return self._rm
    @property
    def ra(self):
        return self._ra
    @property
    def taum(self):
        return self._taum
    @property
    def Za(self):
        return self._Za
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
        w = 2*np.pi*F
        sqrt = np.sqrt(1 + 1j*w*self.taum)
        arg = self.L/self._lambda_DC * sqrt
        den = sqrt * np.sinh(arg)
        self._Zm = self._r_inf / den
        self._Za = self._r_inf * (np.cosh(arg)-1) / den
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
            self._Zp = (Z+self.Za)*self.Zm/(Z+self.Za+self.Zm)
        self._Zload = self.Za + self.Zp
        
    def compute_attenuations(self):
        n_children = len(self.children)
        self._A = np.zeros(n_children, dtype=complex)
        for i,child in enumerate(self.children):
            self._A[i] = 1 + (self.Za + child.Za) / child.Zp
            child.compute_attenuations()


class ImpedanceNode (BaseImpedanceNode):
    def __init__(self, seg, parent=None, children=None):
        self._seg = seg
        sec,x = seg.sec, seg.x
        # the section that contains this segment
        self._sec = sec
        # the position in the containing section
        self._x = x
        super().__init__(ID=ImpedanceNode.make_ID(seg),
                         L=sec.L/sec.nseg*1e-4, # [cm]
                         diam=seg.diam*1e-4,    # [cm]
                         cm=seg.cm*1e-6,        # [F/cm2]
                         rm=1/seg.g_pas,        # [Ohm.cm2]
                         ra=sec.Ra,             # [Ohm.cm]
                         parent=parent,
                         children=children)

    def make_ID(seg):
        return '{}-{:.4f}'.format(seg.sec.name(), seg.x)

    @property
    def seg(self):
        return self._seg
    @property
    def sec(self):
        return self._sec

    def __eq__(self, other):
        return other is not None and other.seg == self.seg


class SWCImpedanceNode (BaseImpedanceNode):
    def __init__(self, ID, coords, diams, cm, rm, ra, node_type, parent=None, children=None):
        L,diam,r1,r2,h,g,S = SWCImpedanceNode.compute_equivalent_cylinder_pars(coords, diams)
        self._r1,self._r2 = r1,r2
        self._h,self._g,self._S = h,g,S
        self._x,self._y,self._z = np.mean(coords, axis=0)
        self._xyz = np.array([self._x, self._y, self._z])
        self._node_type = node_type
        super().__init__(ID,
                         L*1e-4,      # [cm]
                         #TODO: figure out why the multiplication by 2 below is necessary
                         2*diam*1e-4, # [cm]
                         cm*1e-6,     # [F/cm2]
                         rm,          # [Ohm.cm2]
                         ra,          # [Ohm.cm]
                         parent=parent,
                         children=children)

    def compute_coords(self, units='um'):
        for child,A in zip(self.children,self.A):
            v = child._xyz - self._xyz            # direction vector
            if units == 'lambda':
                coeff = 1 / (self._lambda_DC * 1e4)  # [um]
            elif units == 'um':
                coeff = 1
            else:
                raise ValueError(f"Unknown units '{units}'")
            child._XYZ = self._XYZ + np.abs(A) * v * coeff
            child.compute_coords(units)

    @classmethod
    def compute_equivalent_cylinder_pars(cls, coords, diams):
        # radius of the larger base
        r1 = np.max(diams) / 2
        # radius of the smaller base
        r2 = np.min(diams) / 2
        # distance between bases (i.e., height of the truncated cone)
        h = np.sqrt(np.sum(np.diff(np.array(coords),axis=0)**2))
        # slanted height of the truncated cone
        g = np.sqrt(h**2 + (r1-r2)**2)
        S = np.pi * (r1+r2) * g
        L,diam = h,S/(np.pi*h)
        return L,diam,r1,r2,h,g,S


class SWCNode (Node):
    def __init__(self, x, y, z, diam, node_type, ID, parent=None, children=None):
        super().__init__(ID=ID,
                         parent=parent,
                         children=children)
        self._x,self._y,self._z = x,y,z
        self._xyz = np.array([x,y,z])
        self._diam = diam
        self._node_type = node_type

    @property
    def type(self):
        return self._node_type
    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, value):
        self._x = value
        self._xyz[0] = value
    @property
    def y(self):
        return self._y
    @y.setter
    def y(self, value):
        self._y = value
        self._xyz[1] = value
    @property
    def z(self):
        return self._z
    @z.setter
    def z(self, value):
        self._z = value
        self._xyz[2] = value
    @property
    def diam(self):
        return self._diam
    @diam.setter
    def diam(self, value):
        if value <= 0:
            raise Exception('Diameter must be > 0')
        self._diam = value
    @property
    def xyz(self):
        return self._xyz
    @xyz.setter
    def xyz(self, value):
        self._xyz = value
        self._x,self._y,self._z = value

