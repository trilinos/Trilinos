# @HEADER
# *****************************************************************************
#          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
#
# Copyright 2022 NTESS and the PyTrilinos2 contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

import glob
import os
import sys
import numpy as np

def get_list_of_ETI_files_to_include(list_all_ETI_files, list_all_classes_to_ETI):
    list_ETI_files = []
    for ETI_file in list_all_ETI_files:
        for ETI_class in list_all_classes_to_ETI:
            if ETI_file.startswith(ETI_class):
                list_ETI_files.append(ETI_file)
                break
    return list_ETI_files


def write_ETI_include_file(source_dir, filename, list_ETI_files):
    with open(source_dir+'/'+filename, 'w') as fh:
        fh.write('#ifndef PYTRILINOS2_TPETRA_ETI\n')
        fh.write('#define PYTRILINOS2_TPETRA_ETI\n\n')

        for ETI_file in list_ETI_files:
            fh.write('#include <'+ETI_file+'>\n')

        fh.write('\n#endif // PYTRILINOS2_TPETRA_ETI\n')


def write_ETI_getTpetraTypeName_file(source_dir, filename, list_ETI_files):
    with open(source_dir+'/'+filename, 'w') as fh:
        fh.write('from PyTrilinos2.PyTrilinos2 import Tpetra\n\n')

        all_scalar_types = []
        all_local_ordinal_types = []
        all_global_ordinal_types = []
        all_node_types = []

        for ETI_file in list_ETI_files:
            tmp = ETI_file[:ETI_file.index('.')].split("_")

            if tmp[0].lower() == 'tpetra':
                scalar_type = tmp[2].lower()
                scalar_type_internal = scalar_type
                if scalar_type == 'long':
                    scalar_type = 'long long'

                node_type = tmp[-1].lower()
                global_ordinal_type = tmp[-2].lower()
                if global_ordinal_type == 'long':
                    global_ordinal_type = 'long long'
                    local_ordinal_type = tmp[-4].lower()
                else:
                    local_ordinal_type = tmp[-3].lower()
                if local_ordinal_type == 'long':
                    local_ordinal_type = 'long long'

                if not scalar_type in all_scalar_types:
                    all_scalar_types = np.append(all_scalar_types, scalar_type)
                if not local_ordinal_type in all_local_ordinal_types:
                    all_local_ordinal_types = np.append(all_local_ordinal_types, local_ordinal_type)
                if not global_ordinal_type in all_global_ordinal_types:
                    all_global_ordinal_types = np.append(all_global_ordinal_types, global_ordinal_type)
                if not node_type in all_node_types:
                    all_node_types = np.append(all_node_types, node_type)

        default_scalar_type = all_scalar_types[0]
        if 'double' in all_scalar_types:
            default_scalar_type = 'double'

        default_local_ordinal_type = all_local_ordinal_types[0]
        if 'int' in all_local_ordinal_types:
            default_local_ordinal_type = 'int'

        default_global_ordinal_type = all_global_ordinal_types[0]
        if 'long long' in all_global_ordinal_types:
            default_global_ordinal_type = 'long long'

        default_node_type = all_node_types[0]
        if 'cuda' in all_node_types:
            default_node_type = 'cuda'
        elif 'cudaUVM' in all_node_types:
            default_node_type = 'cudaUVM'
        elif 'openmp' in all_node_types:
            default_node_type = 'openmp'
        elif 'threads' in all_node_types:
            default_node_type = 'threads'

        fh.write('def getDefaultScalarType():\n')
        fh.write('\treturn "'+default_scalar_type+'"\n\n')

        fh.write('def getDefaultLocalOrdinalType():\n')
        fh.write('\treturn "'+default_local_ordinal_type+'"\n\n')

        fh.write('def getDefaultGlobalOrdinalType():\n')
        fh.write('\treturn "'+default_global_ordinal_type+'"\n\n')

        fh.write('def getDefaultNodeType():\n')
        fh.write('\treturn "'+default_node_type+'"\n\n')

        fh.write('def getTypeName(class_name, scalar_type=getDefaultScalarType(), local_ordinal_type=getDefaultLocalOrdinalType(), global_ordinal_type=getDefaultGlobalOrdinalType(), node_type=getDefaultNodeType()):\n')

        related_classes = {}

        for ETI_file in list_ETI_files:
            tmp = ETI_file[:ETI_file.index('.')].split("_")

            if tmp[0].lower() == 'tpetra':

                class_name = tmp[1].lower()
                class_name_internal = class_name
                if class_name == 'vector':
                    class_name_internal = 'Vector'
                if class_name == 'multivector':
                    class_name_internal = 'MultiVector'
                if class_name == 'crsmatrix':
                    class_name_internal = 'CrsMatrix'
                scalar_type = tmp[2].lower()
                scalar_type_internal = scalar_type
                if scalar_type == 'long':
                    scalar_type = 'long long'
                    scalar_type_internal = 'long_long'

                node_type = tmp[-1].lower()
                global_ordinal_type = tmp[-2].lower()
                global_ordinal_type_internal = global_ordinal_type
                if global_ordinal_type == 'long':
                    global_ordinal_type = 'long long'
                    global_ordinal_type_internal = 'long_long'
                    local_ordinal_type = tmp[-4].lower()
                else:
                    local_ordinal_type = tmp[-3].lower()
                local_ordinal_type_internal = local_ordinal_type
                if local_ordinal_type == 'long':
                    local_ordinal_type = 'long long'
                    local_ordinal_type_internal = 'long_long'

                if node_type == 'serial':
                    node_type_internal = 'Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Serial_Kokkos_HostSpace'
                if node_type == 'threads':
                    node_type_internal = 'Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Threads_Kokkos_HostSpace'
                if node_type == 'openmp':
                    node_type_internal = 'Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_OpenMP_Kokkos_HostSpace'
                if node_type == 'cuda':
                    node_type_internal = 'Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Cuda_Kokkos_CudaSpace'
                if node_type == 'cudaUVM':
                    node_type_internal = 'Tpetra_KokkosCompat_KokkosDeviceWrapperNode_Kokkos_Cuda_Kokkos_CudaUVMSpace'

                fh.write('\tif class_name.lower() == "'+class_name+'" and scalar_type.lower() == "'+scalar_type+'" and local_ordinal_type.lower() == "'+local_ordinal_type+'" and global_ordinal_type.lower() == "'+global_ordinal_type+'" and node_type.lower() == "'+node_type+'":\n')
                fh.write('\t\treturn Tpetra.'+class_name_internal+'_'+scalar_type_internal+'_'+local_ordinal_type_internal+'_'+global_ordinal_type_internal+'_'+node_type_internal+'_t\n')
                if class_name == 'vector':
                    related_class_name = 'Map_'+local_ordinal_type_internal+'_'+global_ordinal_type_internal+'_'+node_type_internal+'_t'
                    if not related_class_name in related_classes:
                        # Need to add the Map
                        fh.write('\tif class_name.lower() == "map" and local_ordinal_type.lower() == "'+local_ordinal_type+'" and global_ordinal_type.lower() == "'+global_ordinal_type+'" and node_type.lower() == "'+node_type+'":\n')
                        fh.write('\t\treturn Tpetra.'+related_class_name+'\n')
                        related_classes[related_class_name] = True
                    related_class_name = 'Export_'+local_ordinal_type_internal+'_'+global_ordinal_type_internal+'_'+node_type_internal+'_t'
                    if not related_class_name in related_classes:
                        # Need to add the Export
                        fh.write('\tif class_name.lower() == "export" and local_ordinal_type.lower() == "'+local_ordinal_type+'" and global_ordinal_type.lower() == "'+global_ordinal_type+'" and node_type.lower() == "'+node_type+'":\n')
                        fh.write('\t\treturn Tpetra.'+related_class_name+'\n')
                        related_classes[related_class_name] = True
                    related_class_name = 'Import_'+local_ordinal_type_internal+'_'+global_ordinal_type_internal+'_'+node_type_internal+'_t'
                    if not related_class_name in related_classes:
                        # Need to add the Import
                        fh.write('\tif class_name.lower() == "import" and local_ordinal_type.lower() == "'+local_ordinal_type+'" and global_ordinal_type.lower() == "'+global_ordinal_type+'" and node_type.lower() == "'+node_type+'":\n')
                        fh.write('\t\treturn Tpetra.'+related_class_name+'\n')
                        related_classes[related_class_name] = True
                if class_name == 'crsmatrix':
                    related_class_name = 'CrsGraph_'+local_ordinal_type_internal+'_'+global_ordinal_type_internal+'_'+node_type_internal+'_t'
                    if not related_class_name in related_classes:
                        # Need to add the Map
                        fh.write('\tif class_name.lower() == "crsgraph" and local_ordinal_type.lower() == "'+local_ordinal_type+'" and global_ordinal_type.lower() == "'+global_ordinal_type+'" and node_type.lower() == "'+node_type+'":\n')
                        fh.write('\t\treturn Tpetra.'+related_class_name+'\n')
                        related_classes[related_class_name] = True
        fh.write('\tprint("Warning: Unknown type, the function returns None.")\n')
        fh.write('\treturn None\n\n')


if __name__ == '__main__':
    CMAKE_CURRENT_SOURCE_DIR = sys.argv[1]
    list_all_ETI_files = sys.argv[2]
    list_all_classes_to_ETI = sys.argv[3]
    output_file = sys.argv[4]

    with open(list_all_ETI_files, 'r') as fh:
        all_ETI_files = fh.read().splitlines()

    with open(list_all_classes_to_ETI, 'r') as fh:
        all_ETI_classes = fh.read().splitlines()

    reduce_list = True

    if reduce_list:
        print('all_ETI_files = '+str(all_ETI_files))
        print('all_ETI_classes = '+str(all_ETI_classes))
        list_ETI_files = get_list_of_ETI_files_to_include(all_ETI_files, all_ETI_classes)
        print('list_ETI_files = '+str(list_ETI_files))
    else:
        list_ETI_files = all_ETI_files
    write_ETI_include_file(CMAKE_CURRENT_SOURCE_DIR, output_file, list_ETI_files)
    write_ETI_getTpetraTypeName_file(CMAKE_CURRENT_SOURCE_DIR+'/python',  'getTpetraTypeName.py', list_ETI_files)
