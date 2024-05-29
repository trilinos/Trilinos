import re
import sys
import pathlib
from find_files import find_files
from compile_json import compile_json

if __name__ == '__main__':

    assert( len(sys.argv)>2 )

    rol_root = pathlib.Path(sys.argv[1])
#    rol_root = pathlib.Path('/Users/gvonwin/Projects/github/ROL-Trilinos/packages/rol')
    binary_dir = pathlib.Path(sys.argv[2])#rol_root/'rol_parameters'
    rol_src = rol_root/'src'

    # Create list of all (relative path) header files containing the token `ParameterList` in the C++ source 
    relative_pathfiles = find_files(rol_src,'ParameterList','*.hpp')

    # Breakdown of the `sublist` search pattern:
    # \b        : Asserts a word boundary, ensuring that "sublist" is matched as a whole word.
    # sublist   : Matches the literal string "sublist".
    # \s*       : Matches zero or more whitespace characters.
    # \(        : Matches a literal opening parenthesis (().
    # "([^"]+)" : Capturing group that matches one or more characters that are not double quotes ("),
    #             capturing the content between double quotes.
    # \)        : Matches a literal closing parenthesis ()).
    sublist_pattern = re.compile(r'\bsublist\s*\(\s*"([^"]+)"\s*\)', re.MULTILINE)
    sublist_json = compile_json(sublist_pattern,rol_src,relative_pathfiles)

    with open(binary_dir / 'sublist.json', 'w') as f:
        f.write(sublist_json)

    # Breakdown of the `getkey` search pattern:
    # \b        : Asserts a word boundary, ensuring that "get" is matched as a whole word.
    # get       : Matches the literal string "sublist".
    # \s*       : Matches zero or more whitespace characters.
    # \(        : Matches a literal opening parenthesis (().
    # "([^"]+)" : Capturing group that matches one or more characters that are not double quotes ("),
    #             capturing the content between double quotes.
    # ,         : Matches a literal comma.
    # \)        : Matches a literal closing parenthesis ()).
    # ;         : Matches a literal semicolon
    getkey_pattern = re.compile(rf'\bget\s*\(\s*"([^"]*)"\s*,.*\)\s*;', re.MULTILINE)
    getkey_json = compile_json(getkey_pattern,rol_src,relative_pathfiles)

    with open(binary_dir / 'getkey.json', 'w') as f:
        f.write(getkey_json)




