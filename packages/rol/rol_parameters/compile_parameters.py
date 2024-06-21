import re
import pathlib
import subprocess
from typing import Set, Optional, List, Tuple
from read_cpp_source import read_cpp_source

# Compile regex patterns once
SUBLIST_PATTERN = re.compile(r'\bsublist\s*\(\s*"([^"]+)"\s*\)', re.MULTILINE)
GET_KEY_PATTERN = re.compile(r'\bget\s*\(\s*"([^"]*)"\s*,.*\)\s*;', re.MULTILINE)
SET_KEY_PATTERN = re.compile(r'\bset\s*\(\s*"([^"]*)"\s*,.*\)\s*;', re.MULTILINE)

def find_instances(root_path: pathlib.Path,
                   search_token: str,
                   include: Optional[str|Set[str]] = None,
                   exclude: Optional[str|Set[str]] = None,
                   exclude_dir: Optional[str|Set[str]] = None) -> Set[pathlib.Path]:

    # Ensure the root path is an existant directory
    assert( root_path.exists() )
    assert( root_path.is_dir() )

    def join(arg):
        if isinstance(arg,str):
            return [arg]
        else:
            return list(arg)

    cmd = ['grep','-r',search_token]

    if include is not None:
        for inc in join(include):
            cmd.append(f'--include={inc}')

    if exclude is not None:
        for exc in join(exclude):
            cmd.append(f'--exclude={exc}')

    if exclude_dir is not None:
        for exc_dir in join(exclude_dir):
            cmd.append(f'--exclude-dir={exc_dir}')

    cmd.append(str(root_path))

    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

    # Check if the command was successful
    if result.returncode != 0:
        raise Exception(f"Error executing grep: {result.stderr}")

    make_relative = lambda path_str : pathlib.Path(path_str).relative_to(root_path,walk_up=True)

    files = { make_relative(line.split(':')[0]) for line in result.stdout.splitlines() }
    return files



def parse_cpp_file(file_path: pathlib.Path) -> Set[Tuple[str, ...]]:
    cpp = read_cpp_source(file_path)
    cpp = re.sub(';', '\n', cpp)

    def has_token(line: str) -> bool:
        return ('sublist(' in line) or ('get(' in line) or ('set(' in line)

    lines = [re.sub(r'\s+', ' ', line).strip() for line in cpp.splitlines() if has_token(line) and '"' in line]

    names = {}
    code = []
    instances = set()

    for line in lines:
        line = re.sub(r'->', '.', line)
        if '&' in line:
            assignment = line.split('&')[1].strip()
            lhs, rhs = assignment.split('=')
            names[lhs.strip()] = rhs.strip().split('.')
        else:
            code.append(line.strip() if '=' not in line else line.split('=')[1].strip())

    for k, v in names.items():
        if v[0] in names:
            names[k] = names[v[0]] + v[1:]
        for c in code:
            elem = c.split('.')
            if elem[0] in names:
                elem = names[elem[0]] + elem[1:]
            if len(elem) > 1:
                tpl = tuple(filter(has_token, elem))
                if all((e.count('"') in [2, 4]) for e in tpl):
                    instances.add(tuple(e.split('"')[1] for e in tpl))

    return instances





def write_to_csv(instances: List[Tuple[str, ...]], output_file: str):
    with open(output_file, 'w') as f:
        for line in instances:
            f.write(','.join(line) + '\n')

def main():
    rol_src = pathlib.Path('/Users/gvonwin/Projects/github/ROL-Trilinos/packages/rol/src')

    relative_filepaths = find_instances(rol_src, 'ParameterList',
                                        include={'*.hpp', '*.cpp'},
                                        exclude_dir={'compatibility', 'step', 'zoo'})
    all_instances = set()
    for filepath in relative_filepaths:
        all_instances.update(parse_cpp_file(rol_src / filepath))

    write_to_csv(sorted(all_instances), 'all_parameters.csv')

if __name__ == '__main__':
    main()
