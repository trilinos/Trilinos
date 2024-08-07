import re
import subprocess
import pathlib
from pprint import pprint
from typing import Set, Optional



def run_grep_command(src_directory):
    grep_command = [
        'grep',
        '-rE',
        '-e',
        r'(\.|\->)\s*(((s|g)et\s*\(\s*"([a-zA-Z0-9]|\s)+"\s*,\s*\S+\s*\))|sublist)',
        '-e',
        r'(\.|\->)\s*sublist\s*\(\s*\"',
        src_directory
    ]

    try:
        result = subprocess.run(grep_command, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        return e.stderr


def split_cpp_code(code_string):
    # Use a regular expression to split on both '.' and '->'
    # The regex looks for either '->' or '.' as delimiters
    split_pattern = r'->|\.'

    # Split the string and discard the delimiters
    tokens = re.split(split_pattern, code_string)

    # Remove any empty strings from the result and strip whitespace
    tokens = [token.strip() for token in tokens if token.strip()]

    return tokens


def extract_quoted_substring(input_string):
    # Regular expression pattern to match content between double quotes
    pattern = r'"([^"]*)"'

    # Search for the pattern in the input string
    match = re.search(pattern, input_string)

    if match:
        # If a match is found, return the content between the quotes
        return match.group(1)
    else:
        # If no match is found, return None or an empty string
        return None  # or return "" if you prefer











def extract_quoted_strings(string_list):
    return tuple((s.strip('"') for s in string_list if s.startswith('"') and s.endswith('"')))

def custom_sort_key(sublist):
    return sublist[:len(sublist)]

def sort_list_of_lists(list_of_lists):
    return sorted(list_of_lists, key=custom_sort_key)




def parse_cpp_strings(input_list):
    parsed_list = []
    
    for item in input_list:
        # Match a word without parentheses, a quoted string inside parentheses,
        # or a quoted string as the first argument of get() or set()
        match = re.search(r'(\w+)$|"([^"]*)"|\b(?:get|set)\s*\(\s*"([^"]*)"', item)
        if match:
            if match.group(1):  # If it's a word without parentheses
                parsed_list.append(match.group(1))
            elif match.group(2):  # If it's a quoted string inside parentheses
                parsed_list.append(f'"{match.group(2)}"')
            elif match.group(3):  # If it's a quoted string in get() or set()
                parsed_list.append(f'"{match.group(3)}"')

    return parsed_list

def build_hierarchy(data):
    def resolve_list(value_list):
        if not value_list:
            return value_list

        first_item = value_list[0]
        if first_item in data and not first_item.startswith('"'):
            return resolve_list(data[first_item]) + value_list[1:]
        else:
            return [first_item] + resolve_list(value_list[1:])

    return {key: resolve_list(value) for key, value in data.items()}

def build_list_hierarchy(data_dict, input_lists):
    def resolve_list(value_list):
        if not value_list:
            return value_list

        first_item = value_list[0]
        if first_item in data_dict and not first_item.startswith('"'):
            return resolve_list(data_dict[first_item]) + value_list[1:]
        else:
            return [first_item] + resolve_list(value_list[1:])

    return [resolve_list(sublist) for sublist in input_lists]

def create_hierarchical_dict(list_of_lists):
    result = {}
    for path in list_of_lists:
        current = result
        for key in path[:-1]:
            if key not in current:
                current[key] = {}
            current = current[key]
        current[path[-1]] = {}
    return result

if __name__ == '__main__':

    # Every line contains an instance calling at least one of the three functions:
    #
    # - ParameterList::sublist
    # - ParameterList::get
    # - ParameterList::set
    #
    # 1) Defining a local sublist variable
    # 2) Getting a parameter
    # 3) Setting a parameter

    rol_src = pathlib.Path('/Users/gvonwin/Projects/github/ROL-Trilinos/packages/rol/src')

    make_relative = lambda path_str : pathlib.Path(path_str).relative_to(rol_src,walk_up=True)

    strip_excess_whitespace = lambda text : re.sub(r'\s+',' ',text).strip()

    sublist_pattern = re.compile(r'\bsublist\s*\(\s*"([^"]+)"\s*\)')

    data = dict()

    exclusions = ['compatibility','step','zoo']

    local_sublist_pattern = re.compile(r'[ParameterList|auto]\s*[&]\s*(\w+)\s*=\s*(\w+)[\.|\->](.*)')

    output = run_grep_command(rol_src)
    for line in output.splitlines():
        splitline = line.split(':')
        file = str(make_relative(splitline[0]))
        code = strip_excess_whitespace(':'.join(splitline[1:]))
        if not any(f'{e}/' in file for e in exclusions):
            if file not in data.keys():
                data[file] = [code]
            else:
                data[file].append(code)

#    with open('list_of_rol_files.txt','w') as f:
#        f.write('\n'.join(sorted(data.keys())))

    paramset = set()

    for file, code in data.items():
#        print(f'{file}')
        sublist = dict()
        parameters = list()
        for line in code:
#             print(line)
             # Look for locally defined sublists
             match = re.search(local_sublist_pattern,line)
             if match:
                 sublist[match.group(1)] = [match.group(2)] + parse_cpp_strings( split_cpp_code(match.group(3)))
             else:
                 if '=' in line:
                     line = line.split('=')[1].strip()
                 parameters.append(parse_cpp_strings(split_cpp_code(line)))
        sublist = build_hierarchy(sublist)
#        print(sublist)
        parameters = build_list_hierarchy(sublist,parameters)
        [ paramset.add(tuple(p)) for p in map(extract_quoted_strings,parameters)]

    parameters = sorted(filter(len,map(list,paramset)))

#        for p in parameters:
#            print(p)

    parameters = create_hierarchical_dict(parameters)
    

#    pprint(parameters)
#        for p in paramset:
#            print(p)
