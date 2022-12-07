__version__ = '0.1.0'
def version():
    return 'PyROL version: ' + __version__


import importlib
from . getTypeName import getTypeName, getDefaultScalarType, ROL_classes

def getWrapper(classname):
    def wrapper(scalarType=getDefaultScalarType()):
        actualClass = getTypeName(classname, scalarType)
        return actualClass

    return wrapper

defaultScalarType = getDefaultScalarType()
for classnameLong, _ in ROL_classes:
    pos = classnameLong.find('_'+defaultScalarType+'_t')
    if pos <= 0:
        continue
    classname = classnameLong[:pos]


    locals().update({classname: getTypeName(classname, defaultScalarType)})
    locals().update({classname+'_forScalar': getWrapper(classname)})
