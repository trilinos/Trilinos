from PyROL.PyROL import ROL
import inspect
import sys

ROL_classe_names = [cls_name for cls_name, cls_obj in inspect.getmembers(sys.modules['PyROL.PyROL.ROL']) if inspect.isclass(cls_obj)]
ROL_classes = [cls_obj for cls_name, cls_obj in inspect.getmembers(sys.modules['PyROL.PyROL.ROL']) if inspect.isclass(cls_obj)]

def getDefaultScalarType():
    return "double"

def getTypeName(class_name, scalar_type=getDefaultScalarType()):
    class_name_scalar = class_name.lower()+"_"+scalar_type.lower()
    for i in range(0, len(ROL_classe_names)):
        if ROL_classe_names[i].lower().find(class_name_scalar) == 0:
            return ROL_classes[i]
    print("Warning: Unknown type, the function returns None.")
    return None
