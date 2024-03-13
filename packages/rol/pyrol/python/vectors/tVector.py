from mpi4py import MPI
from PyROL.PyROL import ROL
from PyTrilinos2.PyTrilinos2 import Teuchos
from PyTrilinos2.getTpetraTypeName import *

class tVector(ROL.Vector_double_t):
    def __init__(self, dimension=1, default_value=0., map=None, comm=None, scalar_type=getDefaultScalarType(), local_ordinal_type=getDefaultLocalOrdinalType(), global_ordinal_type=getDefaultGlobalOrdinalType(), node_type=getDefaultNodeType()):
        self.scalar_type = scalar_type
        self.local_ordinal_type = local_ordinal_type
        self.global_ordinal_type = global_ordinal_type
        self.node_type = node_type
        self.mapType = getTypeName('Map', local_ordinal_type=self.local_ordinal_type, global_ordinal_type=self.global_ordinal_type, node_type=self.node_type)
        self.vectorType = getTypeName('Vector', scalar_type=self.scalar_type, local_ordinal_type=self.local_ordinal_type, global_ordinal_type=self.global_ordinal_type, node_type=self.node_type)
        if map is None:
            if comm is None:
                comm = Teuchos.getTeuchosComm(MPI.COMM_WORLD)
            map = self.mapType(dimension, 0, comm)
        self.tvector = self.vectorType(map, False)
        self.tvector.putScalar(default_value)
        super().__init__()

    def plus(self, b):
        self.tvector.update(1., b.tvector, 1.)

    def scale(self, scale_factor):
        self.tvector.scale(scale_factor)

    def dot(self, b):
        return self.tvector.dot(b.tvector)

    def norm(self):
        return self.tvector.norm2()

    def clone(self):
        tmp = tVector(map=self.tvector.getMap(), scalar_type=self.scalar_type, local_ordinal_type=self.local_ordinal_type, global_ordinal_type=self.global_ordinal_type, node_type=self.node_type)
        return tmp

    def axpy(self, scale_factor, x):
        self.tvector.update(scale_factor, x.tvector, 1.)

    def dimension(self):
        return self.tvector.getMap().getGlobalNumElements()

    def setScalar(self, new_value):
        self.tvector.putScalar(new_value)

    def __getitem__(self, index):
        if isinstance( index, int ):
            map = self.tvector.getMap()
            if map.isNodeGlobalElement(index):
                local_index = map.getLocalElement(index)
                view = self.tvector.getLocalViewHost()
                return view[local_index]
        if isinstance( index, slice ):
            map = self.tvector.getMap()
            view = self.tvector.getLocalViewHost()
            global_indices = range(*index.indices(self.dimension()))
            local_indices = np.empty(np.size(global_indices), dtype=int)
            for i in range(0, len(global_indices)):
                if map.isNodeGlobalElement(global_indices[i]):
                    local_indices[i] = map.getLocalElement(global_indices[i])
                else:
                    local_indices[i] = 0
            return view[local_indices]

    def __setitem__(self, index, val):
        if isinstance( index, int ):
            map = self.tvector.getMap()
            if map.isNodeGlobalElement(index):
                self.tvector.replaceGlobalValue(index, val)
        if isinstance( index, slice ):
            map = self.tvector.getMap()
            global_indices = range(*index.indices(self.dimension()))
            for i in range(0, len(global_indices)):
                if map.isNodeGlobalElement(global_indices[i]):
                    if len(val) > 1:
                        self.tvector.replaceGlobalValue(global_indices[i], val[i])
                    else:
                        self.tvector.replaceGlobalValue(global_indices[i], val)

    def reduce(self, op):
        reductionType = op.reductionType()
        raise NotImplementedError(reductionType)
