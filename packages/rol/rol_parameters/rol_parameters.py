import re
import sys
import pathlib
import networkx as nx
from find_files import find_files
from read_cpp_source import read_cpp_source

def compile_list_of_sublists(source_dir,binary_dir):
    files = find_files(source_dir, "sublist", include=['*.hpp'])
    pattern = re.compile(r'\.sublist\("\s*([^"]*)"\s*\)', re.MULTILINE)
    all_sublists = set()
    for file in files:
        cpp = read_cpp_source(file)
        matches = list(re.finditer(pattern, cpp))
        if len(matches):
            for m in matches:
                all_sublists.add(m.group(1).strip())

    outfile = binary_dir/'all_sublists.txt'

    with open(outfile,'w') as f:
        for key in sorted(all_sublists):
            f.write(f'{key}\n')

    print(f'Created file {outfile}')


def compile_list_of_keys(source_dir,binary_dir):
    parlist_files = find_files(source_dir, "ParameterList", include=['*.hpp'],exclude=['zoo'])
    pattern = re.compile(r'(\.get\s*\(\s*"[^"]*"\s*,(.*)\)\s*;)',re.MULTILINE)
    all_keys = set()
    for file in parlist_files:
        cpp = read_cpp_source(file)
        matches = list(re.finditer(pattern,cpp))
        if len(matches):
            for m in matches:
                all_keys.add(m.group(0).split('"')[1].strip())

    outfile = binary_dir/'all_keys.txt'

    with open(outfile,'w') as f:
        for key in sorted(all_keys):
            f.write(f'{key}\n')

    print(f'Created file {outfile}')



if __name__ == '__main__':

    assert( len(sys.argv)>2 )

    rol_root = pathlib.Path(sys.argv[1])
    rol_src = rol_root/'src'

    binary_dir = pathlib.Path(sys.argv[2])

    assert(rol_src.exists())
    assert(rol_src.is_dir())
    assert(binary_dir.exists())
    assert(binary_dir.is_dir())

    compile_list_of_sublists(rol_src,binary_dir)
    compile_list_of_keys(rol_src,binary_dir)


#    file = pathlib.Path("/Users/gvonwin/Projects/github/ROL-Trilinos/packages/rol/src/step/ROL_AugmentedLagrangianStep.hpp")
#    cpp = read_cpp_source(file)

#    scope = crop_to_scope(cpp, "Penalty Parameter Reciprocal Lower Bound")
#    print(scope)

#    pattern = re.compile(rf'\.get\s*\(\s*"([^"]*)"\s*,.*\)\s*;', re.MULTILINE)
#    pattern = re.compile(r'\{[^{}]*\.get\s*\(\s*"[^"]*"\s*,(.*)\)[^{}]\}')
#    pattern = re.compile(r'\{[^{}]*\.get\s*\(\s*"[^"]*"\s*,([^{}]*)\)[^{}]*\}', re.MULTILINE)
#
#
#    match = re.search(pattern,cpp)
#    if match:
#        print(match.group())
#    pattern = re.compile(r'(=[^{}]*\.get\s*\(\s*"[^"]*"\s*,(.*)\)\s*;)',re.MULTILINE)
#    matches = list(re.finditer(pattern,cpp))
#    if len(matches):
#        print(file)
#        for m in matches:
#            print(m.group(0))

#    G = nx.DiGraph()
