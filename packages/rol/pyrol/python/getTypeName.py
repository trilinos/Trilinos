from pyrol.pyrol import ROL
import inspect
import sys

trackedTypes = []

def track(self, *args, **kwargs):
    for arg in args:
        if isinstance(arg, tuple(trackedTypes)):
            self._tracked_constructor_args.append(arg)
    for key in kwargs:
        val = kwargs[key]
        if isinstance(val, tuple(trackedTypes)):
            self._tracked_constructor_args.append(val)

def tracking_method(cls_method):
    def new_method(self, *args, **kwargs):
        track(self, *args, **kwargs)
        return cls_method(self, *args, **kwargs)
    return new_method

def withTrackingConstructor(cls):
    class newCls(cls):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

            self._tracked_constructor_args = []
            
            track(self, *args, **kwargs)
            method_names = [m for m in dir(cls) if callable(getattr(cls, m))]

            for name in method_names:
                if name.startswith('add'):
                    setattr(newCls, name, tracking_method(getattr(cls, name)))

    return newCls


ROL_members = {}
for cls_name, cls_obj in inspect.getmembers(sys.modules['pyrol.pyrol.ROL']):
    if inspect.isclass(cls_obj):
        cls_obj = withTrackingConstructor(cls_obj)
        trackedTypes.append(cls_obj)
        setattr(sys.modules['pyrol.pyrol.ROL'], cls_name, cls_obj)
    ROL_members[cls_name] = (cls_obj, inspect.isclass(cls_obj))


ROL_classes = [(cls_name , cls_obj) for cls_name, cls_obj in inspect.getmembers(sys.modules['pyrol.pyrol.ROL']) if inspect.isclass(cls_obj)]

def getDefaultScalarType():
    return "double"

def getTypeName(class_name, scalar_type=getDefaultScalarType()):
    class_name_scalar = class_name.lower()+"_"+scalar_type.lower()
    for i in range(0, len(ROL_classes)):
        if ROL_classes[i][0].lower().find(class_name_scalar) == 0:
            return ROL_classes[i][1]
    print("Warning: Unknown type \"{}\", the function returns None.".format(class_name))
    return None
