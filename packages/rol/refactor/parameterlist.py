import sys
import os
import re
import subprocess as sp
from pathlib import Path

#  
# regex pattern | meaning
# --------------+----------------------------------------------
# \s*           | arbitrary amount of whitespace including none
# \(            | literal left parenthesis
# \)            | literal right parenthesis 
# \.            | literal period
#
# ParameterList::get("name","value")
#  
# re.search(r'\.\s*get\s*"[a-zA-Z0-9]+"\s*,', source_code)


def get_rol_headers(rol_path : Path, token : str) -> [Path]:
    result = sp.Popen(['grep','-rl',token,'--include=*.hpp',rol_path],stdout=sp.PIPE)
    return [Path(line.decode('utf-8').strip()) for line in result.stdout]


def read_file(pathfile : Path) -> str:
    with open(pathfile,"r") as f:
        text = f.read()
    return text

def strip_cpp_comments( cpp_source : str ) -> str:

    in_string = False
    in_single_line_comment = False
    in_multi_line_comment = False
    result = []
    i = 0

    while i < len(cpp_source):

        # Check for string start/end
        if cpp_source[i] == '"' and not (in_single_line_comment or in_multi_line_comment):
            in_string = not in_string
            result.append(cpp_source[i])
        # Check for single-line comment start

        elif i+1 < len(cpp_source) and cpp_source[i:i+2] == "//" and not (in_string or in_multi_line_comment):
            in_single_line_comment = True
            i += 1  # Skip next character to avoid parsing '/' twice

        # Check for multi-line comment start
        elif i + 1 < len(cpp_source) and cpp_source[i:i+2] == "/*" and not (in_string or in_single_line_comment):
            in_multi_line_comment = True
            i += 1  # Skip next character to avoid parsing '*' twice

        # Check for single-line comment end
        elif in_single_line_comment and cpp_source[i] == "\n":
            in_single_line_comment = False
            result.append(cpp_source[i])  # Include newline in result

        # Check for multi-line comment end
        elif i + 1 < len(cpp_source) and in_multi_line_comment and cpp_source[i:i+2] == "*/":
            in_multi_line_comment = False
            i += 1  # Skip next character to avoid parsing '/' twice

        # Append character if not in a comment
        elif not (in_single_line_comment or in_multi_line_comment):
            result.append(cpp_source[i])

        i += 1

    return ''.join(result)
 

def contains_escaped_quote_advanced(s: str) -> bool:
    i = 0
    while i < len(s):
        if s[i] == '\\':
            backslash_count = 1
            i += 1

            # Count consecutive backslashes
            while i < len(s) and s[i] == '\\':
                backslash_count += 1
                i += 1

            # If there's an odd number of backslashes followed by a quote, then it is escaped
            if i < len(s) and s[i] == '"' and backslash_count % 2 == 1:
                return True
        else:
            i += 1
    return False


if __name__ == '__main__':

    """
    Currently iterates over all ROL header files that contain the token (default: ParameterList)
    then looks for calls to ParameterList::get and prints the entire "line" from the start of the
    line to the end of the statement (semicolon). Ignores code comments
    """

    assert( len(sys.argv) > 1 )
    rol_root_path = sys.argv[1]

    token = 'ParameterList' if len(sys.argv) < 3 else sys.argv[2]

    pattern = re.compile(r'^.*?(\.get\s*\(\s*"[^"]*"\s*,(.*)\)\s*;)',re.MULTILINE)
    headers = get_rol_headers(rol_root_path,token)
    for h in headers:
        cpp = strip_cpp_comments(read_file(h))
        matches = re.finditer(pattern,cpp)
        print(h)
        for m in matches:
            print(m.group(0))

