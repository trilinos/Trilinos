from PyROL.PyROL import ROL
import inspect
import sys

ROL_classes = [(cls_name , cls_obj) for cls_name, cls_obj in inspect.getmembers(sys.modules['PyROL.PyROL.ROL']) if inspect.isclass(cls_obj)]

def getDefaultScalarType():
    return "double"

def getTypeName(class_name, scalar_type=getDefaultScalarType()):
    class_name_scalar = class_name.lower()+"_"+scalar_type.lower()
    for i in range(0, len(ROL_classes)):
        if ROL_classes[i][0].lower().find(class_name_scalar) == 0:
            return ROL_classes[i][1]
    print("Warning: Unknown type, the function returns None.")
    return None
