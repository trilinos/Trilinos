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

def get_without_subfolder(line):
    last_index=-1
    for i in range(len(line)):
        if line[i] == '/':
            last_index = i
        if line[i] == '<':
            first_index = i
    if last_index == -1:
        return line
    return line[:first_index+1]+line[last_index+1:]

def get_angular_include(line, remove_subfolder=False):
    first=True
    i0 = 0
    i1 = 0
    newline = line
    for i in range(len(newline)):
        if newline[i] == '"':
            if first:
                if newline[i+1] == '.':
                    return line
                newline = newline[:i] + '<' + newline[i+1:]
                first = False
                i0 = i+1
            else:
                newline = newline[:i] + '>' + newline[i+1:]
                i1 = i
    if newline[i0:i1] in ['storage_class.h', 'cuda_cc7_asm_atomic_op.inc_predicate', 'cuda_cc7_asm_atomic_fetch_op.inc_predicate']:
        return line
    if remove_subfolder:
        return get_without_subfolder(newline)
    return newline

def make_all_includes(all_include_filename, folders):
    all_includes = []
    for folder in folders:
        for filename in (glob.glob(f'{folder}/**/*.hpp', recursive=True) +
                        glob.glob(f'{folder}/**/*.cpp', recursive=True) +
                        glob.glob(f'{folder}/**/*.h', recursive=True) +
                        glob.glob(f'{folder}/**/*.cc', recursive=True) +
                        glob.glob(f'{folder}/**/*.c', recursive=True)):
            with open(filename, 'r') as fh:
                for line in fh:
                    if line.startswith('#include'):
                        all_includes.append(get_angular_include(line).strip())
    all_includes = list(set(all_includes))
    # This is to ensure that the list is always the same and doesn't
    # depend on the filesystem state.  Not technically necessary, but
    # will cause inconsistent errors without it.
    all_includes.sort()
    with open(all_include_filename, 'w') as fh:
        for include in all_includes:
            fh.write(f'{include}\n')
    return all_include_filename

def make_all_includes_from_filenames(all_include_filename, filenames):
    all_includes = []
    for filename in filenames:
        with open(filename, 'r') as fh:
            for line in fh:
                if line.startswith('#include'):
                    all_includes.append(get_angular_include(line).strip())
    all_includes = list(set(all_includes))
    # This is to ensure that the list is always the same and doesn't
    # depend on the filesystem state.  Not technically necessary, but
    # will cause inconsistent errors without it.
    all_includes.sort()
    with open(all_include_filename, 'w') as fh:
        for include in all_includes:
            fh.write(f'{include}\n')
    return all_include_filename

# https://github.com/RosettaCommons/binder/issues/212

def copy_and_angular_includes(filenames, filenames_witout_dir, to_dir):
    # loops over the files, replace include " " by include < > and write them in the to_dir:
    for i in range(len(filenames)):
        filename = filenames[i]
        filename_witout_dir = filenames_witout_dir[i]
        path=os.path.dirname(filename_witout_dir)
        if not os.path.exists(to_dir+'/'+path):
            os.makedirs(to_dir+'/'+path)
        try:
            with open(filename, 'r') as from_f:
                lines = from_f.readlines()
        except UnicodeDecodeError:
            with open(filename, 'r', encoding='iso-8859-1') as from_f:
                lines = from_f.readlines()
        except PermissionError:
            continue

        with open(to_dir+'/'+filename_witout_dir, 'w') as to_f:
            for line in lines:
                if line.startswith('#include'):
                    line = get_angular_include(line, False)
                to_f.write(f'{line}')

if __name__ == '__main__':
    CMAKE_CURRENT_SOURCE_DIR = sys.argv[1]
    CMAKE_CURRENT_BINARY_DIR = sys.argv[2]
    all_header_list_with_dir = sys.argv[3]
    all_header_list_without_dir = sys.argv[4]
    binder_include_name = sys.argv[5]

    with open(all_header_list_with_dir, 'r') as fh:
        all_include_filenames_with_dir = fh.read().splitlines()

    with open(all_header_list_without_dir, 'r') as fh:
        all_include_filenames_without_dir = fh.read().splitlines()

    copy_and_angular_includes(all_include_filenames_with_dir, all_include_filenames_without_dir, CMAKE_CURRENT_BINARY_DIR+'/include_tmp')
    make_all_includes_from_filenames(binder_include_name, [CMAKE_CURRENT_SOURCE_DIR+'/src/PyTrilinos2_Binder_Input.hpp'])
