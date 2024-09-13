import re
import subprocess
import pathlib
from typing import Dict, List, Tuple
import xml.etree.ElementTree as ET
from xml.dom import minidom

def create_xml_from_dict(dict_data: Dict, root_name: str = "Inputs") -> ET.Element:
    def create_element(name: str, content: Dict) -> ET.Element:
        element = ET.Element(name)
        if isinstance(content, dict):
            for key, value in content.items():
                if key == "Parameters" and value:
                    param_element = ET.SubElement(element, "Parameter")
                    param_element.set("name", "Valid Keys")
                    param_element.set("type", "Array(string)")
                    param_element.set("value", "{" + ",".join(value) + "}")
                elif key == "Sublists":
                    for sublist_name, sublist_content in value.items():
                        sublist_element = create_element("ParameterList", sublist_content)
                        sublist_element.set("name", sublist_name)
                        element.append(sublist_element)
        return element

    root = create_element("ParameterList", dict_data)
    root.set("name", root_name)
    return root

def prettify(elem: ET.Element) -> str:
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ")

def grep_source_files(src_directory: str) -> str:
    grep_command = [
        'grep',
        '-rE',
        '-e', r'(\.|\->)\s*(((s|g)et\s*\(\s*"([a-zA-Z0-9]|\s)+"\s*,\s*\S+\s*\))|sublist)',
        '-e', r'(\.|\->)\s*sublist\s*\(\s*\"',
        src_directory
    ]

    try:
        result = subprocess.run(grep_command, capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error occurred: {e}")
        return e.stderr

def split_cpp_code(code_string: str) -> List[str]:
    tokens = re.split(r'->|\.', code_string)
    return [token.strip() for token in tokens if token.strip()]

def extract_quoted_strings(string_list: List[str]) -> Tuple[str, ...]:
    return tuple(s.strip('"') for s in string_list if s.startswith('"') and s.endswith('"'))

def parse_cpp_strings(input_list: List[str]) -> List[str]:
    parsed_list = []
    for item in input_list:
        match = re.search(r'(\w+)$|"([^"]*)"|\b(?:get|set)\s*\(\s*"([^"]*)"', item)
        if match:
            if match.group(1):
                parsed_list.append(match.group(1))
            elif match.group(2):
                parsed_list.append(f'"{match.group(2)}"')
            elif match.group(3):
                parsed_list.append(f'"{match.group(3)}"')
    return parsed_list

def build_hierarchy(data: Dict[str, List[str]]) -> Dict[str, List[str]]:
    def resolve_list(value_list: List[str]) -> List[str]:
        if not value_list:
            return value_list
        first_item = value_list[0]
        if first_item in data and not first_item.startswith('"'):
            return resolve_list(data[first_item]) + value_list[1:]
        else:
            return [first_item] + resolve_list(value_list[1:])

    return {key: resolve_list(value) for key, value in data.items()}

def create_hierarchical_dict(list_of_lists: List[List[str]]) -> Dict:
    result = {}
    for path in list_of_lists:
        current = result
        for key in path[:-1]:
            current = current.setdefault(key, {})
        current[path[-1]] = {}
    return result

def parse(params: Dict) -> Dict:
    result = {'Parameters': [], 'Sublists': {}}
    for k, v in params.items():
        if isinstance(v, dict):
            if v:
                result['Sublists'][k] = parse(v)
            else:
                result['Parameters'].append(k)
        else:
            result['Parameters'].append(k)
    return result

if __name__ == '__main__':
    rol_src = (pathlib.Path.cwd().parents[0]/'src').resolve()
    make_relative = lambda path_str: pathlib.Path(path_str).relative_to(rol_src, walk_up=True)
    strip_excess_whitespace = lambda text: re.sub(r'\s+', ' ', text).strip()

    local_sublist_pattern = re.compile(r'[ParameterList|auto]\s*[&]\s*(\w+)\s*=\s*(\w+)[\.|\->](.*)')
    exclusions = ['compatibility', 'step', 'zoo', 'sol', 'oed', 'dynamic']

    output = grep_source_files(str(rol_src))
    data = {}
    for line in output.splitlines():
        file, *code_parts = line.split(':')
        file = str(make_relative(file))
        code = strip_excess_whitespace(':'.join(code_parts))
        if not any(f'{e}/' in file for e in exclusions):
            data.setdefault(file, []).append(code)

    paramset = set()
    for file, code in data.items():
        sublist = {}
        parameters = []
        for line in code:
            match = re.search(local_sublist_pattern, line)
            if match:
                sublist[match.group(1)] = [match.group(2)] + parse_cpp_strings(split_cpp_code(match.group(3)))
            else:
                if '=' in line:
                    line = line.split('=')[1].strip()
                parameters.append(parse_cpp_strings(split_cpp_code(line)))
        
        sublist = build_hierarchy(sublist)
        parameters = [build_hierarchy(sublist)[sublist_key] + param[1:] for param in parameters for sublist_key in sublist]
        paramset.update(map(extract_quoted_strings, parameters))

    parameters = create_hierarchical_dict(sorted(filter(len, map(list, paramset))))
    parameters.pop('SOL', None)

    xml_root = create_xml_from_dict(parse(parameters))
    pretty_xml = prettify(xml_root)
    
    with open('rol_parameters.xml', 'w') as f:
        f.write(pretty_xml)
    
    print("XML file 'rol_parameters.xml' has been created.")
