
from .tree import Node, Tree

__all__ = ['ImpedanceNode', 'ImpedanceTree']


class ImpedanceNode (Node):
    def __init__(self, seg, parent=None, children=[]):
        self.seg = seg
        sec,x = seg.sec, seg.x
        # the section that contains this segment
        self._sec = sec
        # the position in the containing section
        self._x = x
        super().__init__(ID=ImpedanceNode.make_ID(seg),
                         parent=parent,
                         children=children)
        # attenuation
        self._A = np.zeros(len(self.children), dtype=complex)
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
        self.compute_impedances(0)

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
        super.__init__(root)

    def compute_impedances(self, F):
        self.root.compute_impedances(F)
        
    def compute_attenuations(self):
        self.root.compute_attenuations()

    def _make_branch(self, sec):
        branch = [Node(seg) for seg in sec]
        n_nodes = len(branch)
        for i in range(n_nodes-1):
            branch[i].add_child(branch[i+1])
        for child in sec.children():
            child_node = self._make_branch(child)
            branch[-1].add_child(child_node)
        return branch[0]

    def compute_attenuation(self, seg_i, seg_j):
        path = super.find_connecting_path(ImpedanceNode.make_ID(seg_i),
                                          ImpedanceNode.make_ID(seg_j))
        A = []
        import ipdb
        ipdb.set_trace()

    def find_connecting_path(self, seg_i, seg_j, with_attenuation=True):
        node_i = self.find_segment(seg_i)
        if node_i is None:
            raise Exception('seg_i not in tree')
        node_j = self.find_segment(seg_j)
        if node_j is None:
            raise Exception('seg_j not in tree')
        path = [node_j]
        A = []
        node = node_j
        while node.parent is not None:
            idx = node.parent.children.index(node)
            A.append(node.parent.A[idx])
            node = node.parent
            path.append(node)
            if node == node_i:
                break
        if path[-1] != node_i:
            raise Exception('Cannot find path between seg_i and seg_j')
        if with_attenuation:
            return path[::-1],np.array(A[::-1])
        return path[::-1]

