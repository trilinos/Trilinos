__version__ = '0.1.0'
def version():
    return 'PyROL version: ' + __version__


import importlib
from . getTypeName import getTypeName, getDefaultScalarType, ROL_classes, ROL_members

def getWrapper(classname):
    def wrapper(scalarType=getDefaultScalarType()):
        actualClass = getTypeName(classname, scalarType)
        return actualClass

    return wrapper

defaultScalarType = getDefaultScalarType()
for classnameLong in ROL_members:
    class_obj, isClass = ROL_members[classnameLong]
    pos = classnameLong.find('_'+defaultScalarType+'_t')
    if pos <= 0:
        locals().update({classnameLong: class_obj})
        continue
    classname = classnameLong[:pos]


    locals().update({classname: getTypeName(classname, defaultScalarType)})
    locals().update({classname+'_forScalar': getWrapper(classname)})
