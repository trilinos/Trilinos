# File rol_parameters.py
import sys
import json
import pathlib
import re
import subprocess
from typing import Dict, List
from collections import defaultdict

def hierarchy_to_json(tuple_set):
    def add_to_hierarchy(hierarchy, path, item_type):
        current = hierarchy
        for part in path[:-1]:
            if part not in current["Sublists"]:
                current["Sublists"][part] = {"Parameters": [], "Sublists": {}}
            current = current["Sublists"][part]
        
        if item_type == "Parameter":
            current["Parameters"].append(path[-1])
        elif item_type == "Sublist":
            if path[-1] not in current["Sublists"]:
                current["Sublists"][path[-1]] = {"Parameters": [], "Sublists": {}}

    root = {"Parameters": [], "Sublists": {}}
    
    for tuple_path in tuple_set:
        path = [item[0] for item in tuple_path]
        item_type = tuple_path[-1][1]
        add_to_hierarchy(root, path, item_type)
    
    return root    

from xml.etree.ElementTree import Element, SubElement, tostring
from xml.dom import minidom

def json_to_xml(json_obj):
    def create_parameter_list(name, data):
        param_list = Element("ParameterList", name=name)
        
        if "Parameters" in data and data["Parameters"]:
            param_str = ",".join(data["Parameters"])
            param = SubElement(param_list, "Parameter", name="Parameters", type="Array(string)", value=f"{{{param_str}}}")
        
        if "Sublists" in data:
            for sublist_name, sublist_data in data["Sublists"].items():
                sub_param_list = create_parameter_list(sublist_name, sublist_data)
                if len(sub_param_list) > 0:
                    param_list.append(sub_param_list)
        
        return param_list

    root = create_parameter_list("ROL Parameters", json_obj)
    
    # Remove empty ParameterLists
    for elem in root.iter("ParameterList"):
        if len(elem) == 0:
            root.remove(elem)
    
    # Convert to string and pretty print
    rough_string = tostring(root, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def pretty_print_xml(xml_str : str):
    def remove_extra_newlines(string : str):
        return '\n'.join(line for line in string.split('\n') if line.strip())
    parsed_xml = minidom.parseString(xml_str)
    return remove_extra_newlines(parsed_xml.toprettyxml(indent="  "))


def replace_whitespace(input_string):
    return re.sub(r'\s{2,}', ' ', input_string)

def dereference(code):
    return re.sub(r'\&','',code)


def extract_quoted_strings(input_string):
    # Use a regular expression to find all substrings in double quotes
    matches = re.findall(r'"(.*?)"', input_string)
    # Convert the list of matches to a tuple
    return tuple(matches)

def get_lines( search_path : pathlib.Path ) -> List[str]:
    try:
        result = subprocess.run(['./find_parameters.sh', search_path], capture_output=True, text=True, check=True)
        return result.stdout.splitlines()
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
    
def prune_tuples(data):
    # Sort the tuples by length in descending order
    sorted_data = sorted(data, key=len, reverse=True)
    
    result = []
    leading_elements_set = set()
    
    for tup in sorted_data:
        # Create a leading element tuple (all but the last element)
        leading_elements = tup[:-1]
        
        # Check if the leading elements are already in the set
        if leading_elements not in leading_elements_set:
            # If not, add the tuple to the result and update the set
            result.append(tup)
            leading_elements_set.add(leading_elements)
    
    return result


if __name__ == '__main__':

    src = sys.argv[1]
    grep_pattern = r'(\.|\->)\s*((s|g)et|sublist)\s*\(\s*\"'
    exclude_dirs = 'compatability,dynamic,interiorpoint,oed,sol,zoo'

    # Get all lines in C++ source in which a ParameterList's sublist, get, or set method is called
    lines = get_lines(src)

    data = defaultdict(list)

    defines_local_sublist = lambda line: line.split('.')[-1].startswith('sublist')

    for line in lines:
        filename,code = line.split(':',1)
        if not code.startswith('//'): # Excluded commented out lines
            data[filename].append(replace_whitespace(code.strip()))

    pattern = re.compile(r'(?:sublist|get|set)\s*\(\s*"([^"]*)"\s*(?:\)|,)')

    expanded_lines = list()

    for key, value in data.items():
        ldefs = dict()
        for line in value:
            if defines_local_sublist(line):
                lhs, rhs = line.split('=')
                if len(dereference(lhs).split()) == 2:
                    var = dereference(lhs).split()[1]
                    ldefs[var] = rhs.strip()[:-1]

        if len(ldefs):
            for cycle in range(len(ldefs)):
                for k,v in ldefs.items():
                    vl,vr = v.split('.',1)
                    if vl in ldefs.keys():
                        ldefs[k] = f'{ldefs[vl]}.{vr}'

             
        for line in value:
            exline = line
            for k,v in ldefs.items():
                exline = re.sub(rf' {k}.',rf' {v}.', exline)
            expanded_lines.append(exline)
                        
    pair_tuples = set()
    for line in expanded_lines:
        if '=' in line:
            line = line.split('=')[1]
        if 'et(' in line and 'SOL' not in line:     
            line = line.split(',')[0]
            token_tuple = extract_quoted_strings(line)
            depth = len(token_tuple)
            if depth > 1:
                type_tuple = ('Sublist',) * (depth-1) + ('Parameter',)
                pair_tuples.add(tuple(zip(token_tuple,type_tuple)))

    rol_json = hierarchy_to_json(set(pair_tuples))
    with open('rol_parameters.json','w') as jsonfile:
        jsonfile.write(json.dumps(rol_json,indent=2))

    xml_output = json_to_xml(rol_json)
    pretty_xml_output = pretty_print_xml(xml_output)
    with open('rol_parameters.xml','w') as xmlfile:
        xmlfile.write(pretty_xml_output)
