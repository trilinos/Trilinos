#!/usr/bin/env python3
from xml.etree import ElementTree as ET
import yaml
import argparse

parser = argparse.ArgumentParser(description="Convert a xml input file to yaml")
parser.add_argument("-o", "--output", dest="output",
                    default="", help="yaml output file")
parser.add_argument('xml', type=str,
                    help='input xml')

options = parser.parse_args()

def xml_dict(node, path="", dic=None):
    if dic == None:
        dic = {}
    if 'name' in node.attrib:
        name = node.attrib['name']
    else:
        name = 'root'
    if node.tag == 'ParameterList':
        dic[name] = {}
        for childnode in node:
            xml_dict(childnode, name, dic[name])
    else:
        if node.attrib['type'] == 'int':
            dic[name] = int(node.attrib['value'])
        elif node.attrib['type'] == 'double':
            dic[name] = float(node.attrib['value'])
        elif node.attrib['type'] == 'bool':
            dic[name] = node.attrib['value'] == 'true'
        else:
            dic[name] = node.attrib['value']
    return dic


xml = ET.parse(options.xml)
root_element = xml.getroot()
dic = xml_dict(root_element)
if options.output == "":
    print(yaml.dump(dic, default_flow_style=False))
else:
    with open(options.output, 'w') as f:
        yaml.dump(dic, f, default_flow_style=False)
