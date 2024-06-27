#!/bin/env python

# @HEADER
# *****************************************************************************
#       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
#
# Copyright 2009 NTESS and the Ifpack2 contributors.
# SPDX-License-Identifier: BSD-3-Clause
# *****************************************************************************
# @HEADER

# Script to parse the Hypre headers:
# Usage:
#  From SOURCE_TREE/packages/ifpack2/src call
#  python ../utils/parseHypre.py $HYPRE_INCLUDE_DIR Ifpack2_HypreParameterMap.hpp


import re
from sys import argv

assert len(argv) == 3

hypreInclude = argv[1]
outputFile = argv[2]

functionInt = re.compile(r'\s*HYPRE_Int\s+(HYPRE_[a-zA-Z]+Set[a-zA-Z]+)\s*\(\s*HYPRE_Solver\s*solver\s*,\s*HYPRE_Int\s*[^\s*]+\s*\)')
functionDouble = re.compile(r'\s*HYPRE_Int\s+(HYPRE_[a-zA-Z]+Set[a-zA-Z]+)\s*\(\s*HYPRE_Solver\s*solver\s*,\s*HYPRE_Real\s*[^\s*]+\s*\)')
functionDoubleInt = re.compile(r'\s*HYPRE_Int\s+(HYPRE_[a-zA-Z]+Set[a-zA-Z]+)\s*\(\s*HYPRE_Solver\s*solver\s*,\s*HYPRE_Real\s*[^\s*]+\s*,\s*HYPRE_Int\s*[^\s*]+\s*\)')
functionIntDouble = re.compile(r'\s*HYPRE_Int\s+(HYPRE_[a-zA-Z]+Set[a-zA-Z]+)\s*\(\s*HYPRE_Solver\s*solver\s*,\s*HYPRE_Int\s*[^\s*]+\s*,\s*HYPRE_Real\s*[^\s*]+\s*\)')
functionIntInt = re.compile(r'\s*HYPRE_Int\s+(HYPRE_[a-zA-Z]+Set[a-zA-Z]+)\s*\(\s*HYPRE_Solver\s*solver\s*,\s*HYPRE_Int\s*[^\s*]+\s*,\s*HYPRE_Int\s*[^\s*]+\s*\)')
functionIntStar = re.compile(r'\s*HYPRE_Int\s+(HYPRE_[a-zA-Z]+Set[a-zA-Z]+)\s*\(\s*HYPRE_Solver\s*solver\s*,\s*HYPRE_Int\s*\*[^\s*]+\s*\)')
functionDoubleStar = re.compile(r'\s*HYPRE_Int\s+(HYPRE_[a-zA-Z]+Set[a-zA-Z]+)\s*\(\s*HYPRE_Solver\s*solver\s*,\s*HYPRE_Real\s*\*[^\s*]+\s*\)')
functionIntIntDoubleDouble = re.compile(r'\s*HYPRE_Int\s+(HYPRE_[a-zA-Z]+Set[a-zA-Z]+)\s*\(\s*HYPRE_Solver\s*solver\s*,\s*HYPRE_Int\s*[^\s*]+\s*,\s*HYPRE_Int\s*[^\s*]+\s*,\s*HYPRE_Real\s*[^\s*]+\s*,\s*HYPRE_Real\s*[^\s*]+\s*\)')
functionIntIntIntDoubleIntInt = re.compile(r'\s*HYPRE_Int\s+(HYPRE_[a-zA-Z]+Set[a-zA-Z]+)\s*\(\s*HYPRE_Solver\s*solver\s*,\s*HYPRE_Int\s*[^\s*]+\s*,\s*HYPRE_Int\s*[^\s*]+\s*,\s*HYPRE_Int\s*[^\s*]+\s*,\s*HYPRE_Real\s*[^\s*]+\s*,\s*HYPRE_Int\s*[^\s*]+\s*,\s*HYPRE_Int\s*[^\s*]+\s*\)')
functionIntStarStar = re.compile(r'\s*HYPRE_Int\s+(HYPRE_[a-zA-Z]+Set[a-zA-Z]+)\s*\(\s*HYPRE_Solver\s*solver\s*,\s*HYPRE_Int\s*\*\*[^\s*]+\s*\)')
functionCharStar = re.compile(r'\s*HYPRE_Int\s+(HYPRE_[a-zA-Z]+Set[a-zA-Z]+)\s*\(\s*HYPRE_Solver\s*solver\s*,\s*char\s*\*[^\s*]+\s*\)')

files = ["HYPRE_IJ_mv.h",
         "HYPRE_parcsr_ls.h",
         "krylov.h",
         "_hypre_parcsr_mv.h",
         "_hypre_IJ_mv.h",
         "HYPRE_parcsr_mv.h",
         "HYPRE.h"]


skip = set(["HYPRE_COGMRESSetSkipRealResidualCheck", "HYPRE_BoomerAMGSetDSLUThreshold", "HYPRE_COGMRESSetRelChange"])

s = ''
for fn in files:
    with open(hypreInclude+'/'+fn, 'r') as f:
        d = f.readlines()
    s += ''.join(d)

with open(outputFile, 'w') as out:
    m = functionInt.findall(s)
    m = list(set(m))
    fs = ['  {{\"{func}\", &{func}}}'.format(func=func) if not func in skip else '//  {{\"{func}\", &{func}}}'.format(func=func) for func in m]
    out.write("const std::map<std::string, int_func> FunctionParameter::hypreMapIntFunc_ = {{\n{}\n}};".format(',\n'.join(fs)))

    out.write('\n\n')
    m = functionDouble.findall(s)
    m = list(set(m))
    fs = ['  {{\"{func}\", &{func}}}'.format(func=func) if not func in skip else '//  {{\"{func}\", &{func}}}'.format(func=func) for func in m]
    out.write("const std::map<std::string, double_func> FunctionParameter::hypreMapDoubleFunc_ = {{\n{}\n}};".format(',\n'.join(fs)))

    out.write('\n\n')
    m = functionDoubleInt.findall(s)
    m = list(set(m))
    fs = ['  {{\"{func}\", &{func}}}'.format(func=func) if not func in skip else '//  {{\"{func}\", &{func}}}'.format(func=func) for func in m]
    out.write("const std::map<std::string, double_int_func> FunctionParameter::hypreMapDoubleIntFunc_ = {{\n{}\n}};".format(',\n'.join(fs)))

    out.write('\n\n')
    m = functionIntDouble.findall(s)
    m = list(set(m))
    fs = ['  {{\"{func}\", &{func}}}'.format(func=func) if not func in skip else '//  {{\"{func}\", &{func}}}'.format(func=func) for func in m]
    out.write("const std::map<std::string, int_double_func> FunctionParameter::hypreMapIntDoubleFunc_ = {{\n{}\n}};".format(',\n'.join(fs)))

    out.write('\n\n')
    m = functionIntInt.findall(s)
    m = list(set(m))
    fs = ['  {{\"{func}\", &{func}}}'.format(func=func) if not func in skip else '//  {{\"{func}\", &{func}}}'.format(func=func) for func in m]
    out.write("const std::map<std::string, int_int_func> FunctionParameter::hypreMapIntIntFunc_ = {{\n{}\n}};".format(',\n'.join(fs)))

    out.write('\n\n')
    m = functionIntStar.findall(s)
    m = list(set(m))
    fs = ['  {{\"{func}\", &{func}}}'.format(func=func) if not func in skip else '//  {{\"{func}\", &{func}}}'.format(func=func) for func in m]
    out.write("const std::map<std::string, int_star_func> FunctionParameter::hypreMapIntStarFunc_ = {{\n{}\n}};".format(',\n'.join(fs)))

    out.write('\n\n')
    m = functionDoubleStar.findall(s)
    m = list(set(m))
    fs = ['  {{\"{func}\", &{func}}}'.format(func=func) if not func in skip else '//  {{\"{func}\", &{func}}}'.format(func=func) for func in m]
    out.write("const std::map<std::string, double_star_func> FunctionParameter::hypreMapDoubleStarFunc_ = {{\n{}\n}};".format(',\n'.join(fs)))

    out.write('\n\n')
    m = functionIntIntDoubleDouble.findall(s)
    m = list(set(m))
    fs = ['  {{\"{func}\", &{func}}}'.format(func=func) if not func in skip else '//  {{\"{func}\", &{func}}}'.format(func=func) for func in m]
    out.write("const std::map<std::string, int_int_double_double_func> FunctionParameter::hypreMapIntIntDoubleDoubleFunc_ = {{\n{}\n}};".format(',\n'.join(fs)))

    out.write('\n\n')
    m = functionIntIntIntDoubleIntInt.findall(s)
    m = list(set(m))
    fs = ['  {{\"{func}\", &{func}}}'.format(func=func) if not func in skip else '//  {{\"{func}\", &{func}}}'.format(func=func) for func in m]
    out.write("const std::map<std::string, int_int_int_double_int_int_func> FunctionParameter::hypreMapIntIntIntDoubleIntIntFunc_ = {{\n{}\n}};".format(',\n'.join(fs)))

    out.write('\n\n')
    m = functionIntStarStar.findall(s)
    m = list(set(m))
    fs = ['  {{\"{func}\", &{func}}}'.format(func=func) if not func in skip else '//  {{\"{func}\", &{func}}}'.format(func=func) for func in m]
    out.write("const std::map<std::string, int_star_star_func> FunctionParameter::hypreMapIntStarStarFunc_ = {{\n{}\n}};".format(',\n'.join(fs)))

    out.write('\n\n')
    m = functionCharStar.findall(s)
    m = list(set(m))
    fs = ['  {{\"{func}\", &{func}}}'.format(func=func) if not func in skip else '//  {{\"{func}\", &{func}}}'.format(func=func) for func in m]
    out.write("const std::map<std::string, char_star_func> FunctionParameter::hypreMapCharStarFunc_ = {{\n{}\n}};".format(',\n'.join(fs)))
