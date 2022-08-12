#!/bin/env python
#

import sys

string_to_replace = sys.argv[1]
repalcement_string = sys.argv[2]
file_name = sys.argv[3]

with open(file_name, 'r') as file:
  lines = file.readlines()

newLines = []
for line in lines:
  line = line.replace(string_to_replace, repalcement_string)
  newLines.append(line)

with open(file_name, 'w') as file:
  file.writelines(newLines)
