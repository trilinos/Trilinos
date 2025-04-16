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
    S = ('Scalar', )

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
                                           ("FECrsMatrix", SLGN)],
        "PyTrilinos2.PyTrilinos2.Stratimikos": [("LinearSolverBuilder", S)],
        "PyTrilinos2.PyTrilinos2.Teuchos": [("ScalarTraits", S)],
        "PyTrilinos2.PyTrilinos2.Thyra": [("TpetraEuclideanScalarProd", SLGN),
                                          ("TpetraLinearOp", SLGN),
                                          ("TpetraMultiVector", SLGN),
                                          ("TpetraOperatorVectorExtraction", SLGN),
                                          ("TpetraVectorSpace", SLGN),
                                          ("TpetraVector", SLGN),
                                          ("XpetraLinearOp", SLGN),]
    }

    replacements = [('_t', ''),
                    ('unsigned_int', 'unsigned int'),
                    ('unsigned_short', 'unsigned short'),
                    ('unsigned_long', 'unsigned long'),
                    ('long_long', 'long long'),
                    ('unsigned_long_long', 'unsigned long long'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Serial_Kokkos_HostSpace', 'serial'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Serial', 'serial'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_OpenMP_Kokkos_HostSpace', 'openmp'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_OpenMP', 'openmp'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Cuda_Kokkos_CudaSpace', 'cuda'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Cuda_Kokkos_CudaUVMSpace', 'cuda_uvm'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_HIP_Kokkos_HIPSpace', 'hip'),
                    ('Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_HIP_Kokkos_HIPManagedSpace', 'hip_managed'),
                    ]


    class Instantiator:
        def __init__(self, class_name, templateParams, available_instantiations):
            self.class_name = class_name
            self.templateParams = templateParams
            self.available_instantiations = available_instantiations
            self.defaults = {k: v for k, v in zip(templateParams, list(available_instantiations.keys())[0])}

        @property
        def templatedClasses(self):
            return tuple([inst[1] for inst in self.available_instantiations.values()])

        def __call__(self, **kwargs):
            k = copy(self.defaults)
            k.update(kwargs)
            args = tuple([k[templateParam] for templateParam in self.templateParams])
            if args in self.available_instantiations:
                instantiation = self.available_instantiations[*args][1]
                return instantiation
            else:
                raise NotImplementedError("No instantiation of \"{}\" available for {}. Available instatiations: {}".format(self.class_name, args, list(self.available_instantiations.keys())))

        def __getitem__(self, args):
            if isinstance(args, str):
                args = (args, )
            kwargs = {k: v for k, v in zip(self.templateParams, args)}
            return self(**kwargs)

        def __repr__(self):
            return "{}{}".format(self.class_name, self.templateParams)

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
            setattr(sys.modules[pkg], class_name, Instantiator(class_name, templateParams, available_instantiations))


setupInstantiationHelpers()
