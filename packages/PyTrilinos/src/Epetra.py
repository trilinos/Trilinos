from UserArray            import *
from PyTrilinos.RawEpetra import *

class Vector(UserArray,NumPyVector):

    def __init__(self, *args):
        NumPyVector.__init__(self, *args)
        UserArray.__init__(self,self.getArray(),'d',0,1)

    def __str__(self):
        return str(self.array)
