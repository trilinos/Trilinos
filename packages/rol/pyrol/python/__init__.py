import importlib

from . getTypeName import getTypeName, getDefaultScalarType, ROL_classes, ROL_members

__version__ = '0.1.0'
def version():
    return 'PyROL version: ' + __version__


def getWrapper(classname):
    def wrapper(scalarType=getDefaultScalarType()):
        actualClass = getTypeName(classname, scalarType)
        return actualClass
    return wrapper


defaultScalarType = getDefaultScalarType()

supported_objects = {"Bounds", "Constraint", "LinearOperator",
                     "LinearConstraint", "Objective", "Problem", "Solver",
                     "Vector", "getCout", "getParametersFromXmlFile"}

unsupported = importlib.import_module('.unsupported', 'pyrol')

for classnameLong in ROL_members:
    class_obj, _ = ROL_members[classnameLong]
    pos = classnameLong.find('_'+defaultScalarType+'_t')
    if pos <= 0:
        if classnameLong in supported_objects:
            locals().update({classnameLong: class_obj})
        else:
            setattr(unsupported, classnameLong, class_obj)
        continue

    classname = classnameLong[:pos]

    if classname in supported_objects:
        locals().update({classname: getTypeName(classname, defaultScalarType)})
        locals().update({classname+'_forScalar': getWrapper(classname)})
    else:
        setattr(unsupported, classname, getTypeName(classname, defaultScalarType))
        setattr(unsupported, classname + '_forScalar', getWrapper(classname))

del importlib, getTypeName, getDefaultScalarType, ROL_classes, ROL_members
del class_obj, classname, classnameLong, getWrapper, pos
