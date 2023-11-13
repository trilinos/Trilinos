#!/usr/bin/env python3

# This script attempts to find unneeded or wrong includes in MueLu and
# prints a list of suggested fixes. Not all of them are correct.
#
# Run in $TRILINOS_SRC/packagesmuelu/src or a subfolder.

from pathlib import Path
import re
import subprocess
from pprint import pprint

# Set to true if the script is supposed to fix things.
# This is an 80% solution. Some of the fixes will break the code.

# takeAction = True
takeAction = False

trilinosPath = Path(subprocess.Popen(['git', 'rev-parse', '--show-toplevel'], stdout=subprocess.PIPE).communicate()[0].rstrip().decode('utf-8'))

p = Path('.')

includeRE = re.compile(r'\s*#\s*include\s+["<]([^"<]*)[">]')
xpetraMueLuHeaderRE = re.compile(r'(MueLu_|Xpetra_)')
forwardRE = re.compile(r'.*_fwd\.hpp')
declRE = re.compile(r'.*_decl\.hpp')
xpetraHeaderRE = re.compile(r'Xpetra_.+\.hpp')
mueluHeaderRE = re.compile(r'MueLu_.+\.hpp')

unETId = set([
    "Xpetra_ConfigDefs.hpp",
    "MueLu_UseShortNames.hpp",
    "MueLu_UseShortNamesOrdinal.hpp",
    "MueLu_Level.hpp",
    "MueLu_Monitor.hpp",
    "MueLu_BaseClass.hpp",
    "MueLu_ConfigDefs.hpp",
    "MueLu_SingleLevelFactoryBase.hpp",
    "MueLu_TwoLevelFactoryBase.hpp",
    "MueLu_FactoryManagerBase.hpp",
    "MueLu_FacadeClassBase.hpp",
    "MueLu_HierarchyManager.hpp",
    "MueLu_SmootherPrototype.hpp",
    "MueLu_GraphBase.hpp",
    "MueLu_AggregationAlgorithmBase.hpp",
    "MueLu_AggregationAlgorithmBase_kokkos.hpp",
    "MueLu_Exceptions.hpp",
    "MueLu_Types.hpp",
    "MueLu_Factory.hpp",
    "MueLu_SmootherBase.hpp",
    "MueLu_ParameterListAcceptor.hpp",
    "MueLu_VerbosityLevel.hpp",
    "MueLu_SolverBase.hpp",
    "MueLu_PreDropFunctionBaseClass.hpp",
    "MueLu_NoFactory.hpp",
    "MueLu_KeepType.hpp",
    "MueLu_FacadeClassFactory.hpp"
])

skip = set([
    "Xpetra_ConfigDefs.hpp",
    "MueLu_ConfigDefs.hpp",
    "MueLu_VerbosityLevel.hpp",
    "MueLu_KeepType.hpp",
    "MueLu_Exceptions.hpp",
    "MueLu_Types.hpp",
    "MueLu_BaseClass.hpp",
    "MueLu_HierarchyUtils.hpp",
    "MueLu_ParameterListUtils.hpp",
    "MueLu_FactoryManagerBase.hpp"
])

ETId = set()
with open(trilinosPath/'packages/muelu/src/Utils/ClassList/SC-LO-GO-NO.classList') as f:
    for line in f.readlines():
        ETId.add(line.replace('\n', '').split(' -')[0])
with open(trilinosPath/'packages/muelu/src/Utils/ClassList/LO-GO-NO.classList') as f:
    for line in f.readlines():
        ETId.add(line.replace('\n', '').split(' -')[0])

for classname in ETId:
    unETId.discard('MueLu_'+classname+'.hpp')


def getIncludes(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()
    includes = set()
    duplicates = set()
    for line in lines:
        m = includeRE.match(line)
        if m:
            include = m.group(1)
            if include in includes and include != "MueLu_UseShortNames.hpp":
                print('{:<80} Duplicate include of \'{}\''.format(str(filename), include))
                duplicates.add(include)
            includes.add(m.group(1))
    return includes, duplicates


# def getParentClasses(filename, className):
#     with open(filename, 'r') as f:
#         text = ''.join(f.readlines())
#     reParents = re.compile('class\s+'+className+'\s*:\s*public\s+([^{]*)')
#     # print(filename)
#     print(className, text.replace('\n', ''))
#     parents = set()
#     for m in re.findall(reParents, text.replace('\n', '')):
#         parents.add(m)
#     return parents


def removeInclude(filename, include):
    with open(filename, 'r') as f:
        lines = f.readlines()
    newLines = []
    includeRE = re.compile(r'\s*#\s*include\s+["<]'+include+'[">]')
    for line in lines:
        m = includeRE.match(line)
        if not m:
            newLines.append(line)
    with open(filename, 'w') as f:
        f.writelines(newLines)


def replaceInclude(filename, include, newInclude):
    with open(filename, 'r') as f:
        lines = f.readlines()
    newLines = []
    includeRE = re.compile(r'\s*#\s*include\s+["<]'+include+'[">]')
    for line in lines:
        m = includeRE.match(line)
        if not m:
            newLines.append(line)
        else:
            newLines.append(line.replace(include, newInclude))
    with open(filename, 'w') as f:
        f.writelines(newLines)


def addInclude(filename, newInclude):
    with open(filename, 'r') as f:
        lines = f.readlines()
    newLines = []
    includeAdded = False
    for line in lines:
        m = includeRE.match(line)
        newLines.append(line)
        if not includeAdded and m:
            newLines.append('#include "{}"\n'.format(newInclude))
            includeAdded = True
    with open(filename, 'w') as f:
        f.writelines(newLines)


def removeDuplicateInclude(filename, include):
    with open(filename, 'r') as f:
        lines = f.readlines()
    includeRE = re.compile(r'\s*#\s*include\s+["<]'+include+'[">]')
    includeFound = False
    newLines = []
    for line in lines:
        m = includeRE.match(line)
        if not m:
            newLines.append(line)
        elif not includeFound:
            newLines.append(line)
            includeFound = True
    with open(filename, 'w') as f:
        f.writelines(newLines)


def usesClass(filename, classname):
    with open(filename, 'r') as f:
        lines = f.readlines()
    for line in lines:
        if line.find(classname) >= 0 and not includeRE.match(line):
            return True
    return False


unclear = set()

includedClasses = set()
classesWithDeclDef = set()

######################################################################
# Loop over _decl files

for decl_file in p.rglob('*_decl.hpp'):
    if decl_file.name == "MueLu_FactoryFactory_decl.hpp":
        continue

    def_file = Path(str(decl_file).replace('_decl.hpp', '_def.hpp'))

    # Need a _def to match the _decl
    if not def_file.exists():
        print('{:<80} No def'.format(str(decl_file)))
        if takeAction:
            decl_file.unlink()
        continue
    else:
        classesWithDeclDef.add(decl_file.name.replace('_decl.hpp', '').replace('MueLu_', ''))

    includes_decl, duplicates_decl = getIncludes(decl_file)
    if takeAction:
        # remove duplicate includes
        for duplicate in duplicates_decl:
            removeDuplicateInclude(decl_file, duplicate)

    includes_def, duplicates_def = getIncludes(def_file)
    if takeAction:
        # remove duplicate includes
        for duplicate in duplicates_def:
            removeDuplicateInclude(def_file, duplicate)

    # loop over includes of decl file
    for include_decl in includes_decl:

        isXpetraMueLu = xpetraMueLuHeaderRE.match(include_decl)
        if not isXpetraMueLu:
            continue

        isFwd = forwardRE.match(include_decl)
        isDecl = declRE.match(include_decl)
        isMueLu = mueluHeaderRE.match(decl_file.name)

        includedClasses.add(include_decl.replace('_decl.hpp', '').replace('_def.hpp', '').replace('_fwd.hpp', '').replace('.hpp', '').replace('MueLu_', '').replace('Xpetra_', ''))

        if isFwd:
            # include_decl is a forward header
            if include_decl.replace('_fwd.hpp', '_decl.hpp') in includes_def:
                # fwd header is included in decl
                continue
            nonFwd = include_decl.replace('_fwd.hpp', '.hpp')

            className = nonFwd.replace('.hpp', '').replace('MueLu_', '').replace('Xpetra_', '')

            if (not usesClass(decl_file, className)) and (not usesClass(def_file, className)):
                print('{:<80} Unused fwd header \'{}\'; \'{}\' not used in decl or def'.format(str(decl_file), include_decl, className))
                if takeAction:
                    removeInclude(decl_file, include_decl)
        elif isDecl:
            className = decl_file.name.replace('_decl.hpp', '').replace('MueLu_', '')
            # print(getParentClasses(decl_file, className))
            fwd = include_decl.replace('_decl.hpp', '_fwd.hpp')
            # print('{:<80} Replace {} with {}'.format(str(decl_file), include_decl, fwd))
            # if takeAction:
            #     replaceInclude(decl_file, include_decl, fwd)
            continue
        elif include_decl in unETId:
            pass
        elif isMueLu and xpetraHeaderRE.match(include_decl):
            fwd = include_decl.replace('.hpp', '_fwd.hpp')
            print('{:<80} Replace \'{}\' with \'{}\''.format(str(decl_file), include_decl, fwd))
            if takeAction:
                replaceInclude(decl_file, include_decl, fwd)
        elif isMueLu and mueluHeaderRE.match(include_decl):
            className = include_decl.replace('.hpp', '').replace('MueLu_', '')
            if className in ETId and not usesClass(decl_file, className):
                if usesClass(def_file, className):
                    fwd = include_decl.replace('.hpp', '_fwd.hpp')
                    print('{:<80} Replace \'{}\' with \'{}\''.format(str(decl_file), include_decl, fwd))
                    if takeAction:
                        replaceInclude(decl_file, include_decl, fwd)
                    if include_decl not in includes_def:
                        print('{:<80} Add \'{}\''.format(str(def_file), include_decl))
                        if takeAction:
                            addInclude(def_file, include_decl)
                else:
                    print('{:<80} Unused include \'{}\'; \'{}\' not used in decl or def'.format(str(decl_file), include_decl, className))
                    if takeAction:
                        removeInclude(decl_file, include_decl)
        else:
            print('{:<80} nonFwd {}'.format(decl_file, include_decl))
            unclear.add(include_decl)


######################################################################
# Loop over _def files

for def_file in p.rglob('*_def.hpp'):
    includes_def, duplicates_def = getIncludes(def_file)
    if takeAction:
        for duplicate in duplicates_def:
            removeDuplicateInclude(def_file, duplicate)

    decl_file = Path(str(def_file).replace('_def.hpp', '_decl.hpp'))
    # Need a _decl to match the _def
    if not decl_file.exists():
        print('{:<80} No decl'.format(str(def_file)))
        if takeAction:
            def_file.unlink()
        continue

    if decl_file.name not in includes_def:
        header = def_file.name.replace('_def.hpp', '.hpp')
        if header in includes_def:
            print('{:<80} Replace \'{}\' with \'{}\''.format(str(def_file), header, decl_file.name))
            if takeAction:
                replaceInclude(def_file, header, decl_file.name)

    for include_def in includes_def:

        isXpetraMueLu = xpetraMueLuHeaderRE.match(include_def)
        if not isXpetraMueLu:
            continue

        isFwd = forwardRE.match(include_def)
        className = include_def.replace('.hpp', '').replace('MueLu_', '')
        isDecl = declRE.match(include_def)
        isMueLu = mueluHeaderRE.match(def_file.name)

        if not isDecl:
            includedClasses.add(include_def.replace('_decl.hpp', '').replace('_def.hpp', '').replace('_fwd.hpp', '').replace('.hpp', '').replace('MueLu_', '').replace('Xpetra_', ''))

        if include_def in skip:
            pass
        elif isFwd:
            print('{:<80} Fwd header \'{}\' in def'.format(str(def_file), include_def))
            if takeAction:
                removeInclude(def_file, include_def)
        elif isDecl:
            pass
        elif mueluHeaderRE.match(include_def) and not usesClass(def_file, className):
            print('{:<80} Unused include \'{}\'; \'{}\' not used in def'.format(str(def_file), include_def, className))
            if takeAction:
                removeInclude(def_file, include_def)



# print(unclear)
# print(ETId)
    # print(decl_file, includes_decl, includes_def)
print('\nNot included by other files:')
pprint(classesWithDeclDef-includedClasses)
