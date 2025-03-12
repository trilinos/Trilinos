# @HEADER
# *****************************************************************************
#          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
#
# Copyright 2022 NTESS and the PyTrilinos2 contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

__version__ = '0.1.0'
def version():
    return 'PyTrilinos2 version: ' + __version__


from . PyTrilinos2 import *


def setupInstantiationHelpers():
    import inspect
    import sys
    from copy import copy

    SLGN = ('Scalar', 'LocalOrdinal', 'GlobalOrdinal', 'Node')
    LGN = ('LocalOrdinal', 'GlobalOrdinal', 'Node')

    templatedClasses = {
        "PyTrilinos2.PyTrilinos2.Tpetra": [("MultiVector", SLGN),
                                           ("FEMultiVector", SLGN),
                                           ("BlockMultiVector", SLGN),
                                           ("Export", LGN),
                                           ("Import", LGN),
                                           ("Vector", SLGN),
                                           ("Map", LGN),
                                           ("Operator", SLGN),
                                           ("RowMatrix", SLGN),
                                           ("CrsGraph", LGN),
                                           ("FECrsGraph", LGN),
                                           ("RowGraph", LGN),
                                           ("CrsMatrix", SLGN),
                                           ("FECrsMatrix", SLGN)]
    }

    replacements = [('long_long', 'long long'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Serial_Kokkos_HostSpace_t', 'serial'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Serial_t', 'serial'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_OpenMP_Kokkos_HostSpace', 'openmp'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_OpenMP', 'openmp'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Cuda_Kokkos_CudaSpace', 'cuda'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Cuda_Kokkos_CudaUVMSpace', 'cuda_uvm'),
                    ]



    def getInstantiation(class_name, templateParams, available_instantiations):

        defaults = {k: v for k, v in zip(templateParams, list(available_instantiations.keys())[0])}

        def instantiator(**kwargs):
            k = copy(defaults)
            k.update(kwargs)
            args = tuple([k[templateParam] for templateParam in templateParams])
            if args in available_instantiations:
                instantiation = available_instantiations[*args][1]
                return instantiation
            else:
                raise NotImplementedError("No instantiation of \"{}\" available for {}. Available instatiations: {}".format(class_name, args, list(available_instantiations.keys())))

        return instantiator


    for pkg in templatedClasses:
        classes = [(cls_name , cls_obj) for cls_name, cls_obj in inspect.getmembers(sys.modules[pkg]) if inspect.isclass(cls_obj)]
        for class_name, templateParams in templatedClasses[pkg]:
            numTemplates = len(templateParams)
            available_instantiations = {}
            matching_instantiations = [c for c in classes if c[0].find(class_name+'_') == 0]
            for c in matching_instantiations:
                d = c[0]
                for r in replacements:
                    d = d.replace(*r)
                components = d.split('_')
                available_instantiations[tuple(components[1:])] = c
            setattr(sys.modules[pkg], class_name, getInstantiation(class_name, templateParams, available_instantiations))


setupInstantiationHelpers()
