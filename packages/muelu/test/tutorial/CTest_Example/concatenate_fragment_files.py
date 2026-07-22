import os
import sys

output_file = sys.argv[1]
num_fragments = int(sys.argv[2])

print("Attempting to concatenate", num_fragments, "file fragments...")

fragment_list = []
for i in range(num_fragments):
    fragment_list.append(sys.argv[3+i])

print(fragment_list)

with open(output_file, 'w') as outfile:
    for fname in fragment_list:
        with open(fname) as infile:
            for line in infile:
                outfile.write(line)
