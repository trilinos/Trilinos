"""Python interface to the Rapid Optimization Library (ROL).

The commonly used ROL classes are available here for the default scalar
type, e.g. pyrol.Objective, and pyrol.<Name>_forScalar(scalarType) returns
the class for another scalar type. All other classes are in
pyrol.unsupported.
"""
import importlib

from . getTypeName import getTypeName, getDefaultScalarType, ROL_classes, ROL_members
from pyrol.pyrol.Teuchos import ParameterList
from ._version import __version__


def version():
    return 'PyROL version: ' + __version__


# getTypeName is bound as a default argument because it is deleted from this
# module's namespace at the end of this file.
def getWrapper(classname, getTypeName=getTypeName):
    def wrapper(scalarType=getDefaultScalarType()):
        return getTypeName(classname, scalarType)
    return wrapper


defaultScalarType = getDefaultScalarType()

supported_objects = {"Bounds",
                     "Constraint",
                     "getCout",
                     "getParametersFromXmlFile",
                     "getParametersFromYamlFile",
                     "LinearConstraint",
                     "LinearOperator",
                     "Objective",
                     "Problem",
                     "Solver",
                     "Vector"}

unsupported = importlib.import_module('.unsupported', 'pyrol')

# Export every ROL member: here if it is in supported_objects, otherwise in
# pyrol.unsupported. Instances for the default scalar type are exported under
# their short name, e.g. Objective_double_t as Objective and
# Objective_forScalar.
for classnameLong in ROL_members:
    class_obj, _ = ROL_members[classnameLong]
    pos = classnameLong.find('_' + defaultScalarType + '_t')
    if pos <= 0:
        if classnameLong in supported_objects:
            globals()[classnameLong] = class_obj
        else:
            setattr(unsupported, classnameLong, class_obj)
        continue

    classname = classnameLong[:pos]

    if classname in supported_objects:
        globals()[classname] = getTypeName(classname, defaultScalarType)
        globals()[classname + '_forScalar'] = getWrapper(classname)
    else:
        setattr(unsupported, classname, getTypeName(classname, defaultScalarType))
        setattr(unsupported, classname + '_forScalar', getWrapper(classname))

del importlib, getTypeName, getDefaultScalarType, ROL_classes, ROL_members
del class_obj, classname, classnameLong, getWrapper, pos
