import re
import os

def remove_if_conditions(file_path, regex_pattern):
    # Compile the regex pattern
    pattern = re.compile(regex_pattern)

    # Read the content of the CMake file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    output_lines = []
    skip_block = 0  # Use an integer to track the nesting level of if blocks
    skip_next_endif = False

    for line in lines:
        stripped_line = line.strip()
        print(stripped_line)
        # Check for the start of an if block
        if stripped_line.lower().startswith(('if(', 'if ')):
            matches = pattern.findall(stripped_line)
            if matches:
                print(f" ---> MATCH {stripped_line}")
                removeme = False
                for match in matches:
                    actual_match = next(filter(None, match))
                    print(stripped_line.split(actual_match))
                    if not stripped_line.split(actual_match)[0].strip().lower().endswith("not"):
                        removeme = True
                        break
                else:
                    #print(" ---> Conditional(s) only preceded by NOT, leave it alone")
                    pass
                if removeme:
                    print(" ---> Removing this if")
                    skip_block += 1  # Increment nesting level
                    continue  # Skip this line
            elif skip_block > 0:
                #print(f" ---> Found a nested if statement")
                skip_block += 1
                continue

        # Check for the end of an if block
        if stripped_line.lower().startswith('endif'):
            if skip_block > 0:
                #print(" ---> Endif, backing off nesting level")
                skip_block -= 1  # Decrement nesting level
                continue  # Skip this line
            if skip_next_endif:
                #print(" ---> Skipping trailing endif")
                skip_next_endif = False
                continue

        # Check for else clause
        if stripped_line.lower().startswith('else'):
            if skip_block == 1:
                #print(" ---> Root-level else found, will skip trailing endif (also backing off nesting level)")
                skip_block -= 1
                skip_next_endif = True
                continue  # Continue to the next line

        # If not skipping, add the line to output
        if skip_block == 0:
            output_lines.append(line)

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(output_lines)

def remove_assert_defineds(file_path, regex_pattern):
    # Compile the regex pattern
    pattern = re.compile(regex_pattern)

    # Read the content of the CMake file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    output_lines = []
    in_assert_defined = False

    for line in lines:
        stripped_line = line.strip()
        # Check for the start of an if block
        if stripped_line.lower().startswith(('assert_defined(', 'assert_defined ')):
            in_assert_defined = True

        if in_assert_defined:
            matches = pattern.findall(line)
            if matches:
                continue

        if ")" in line:
            in_assert_defined = False

        output_lines.append(line)

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(output_lines)


def remove_regex(file_path, regex_pattern):
    pattern = re.compile(regex_pattern)

    with open(file_path, 'r') as file:
        lines = file.readlines()

    output_lines = []

    for line in lines:
        if pattern.search(line):
            continue
        output_lines.append(line)

    with open(file_path, 'w') as file:
        file.writelines(output_lines)

def process_directory(directory):
    deprecated_packages = [
        "Amesos",
        "AztecOO",
        "Epetra",
        "EpetraExt",
        "Ifpack",
        "Intrepid",
        "Isorropia",
        "ML",
        "NewPackage",
        "Pliris",
        "PyTrilinos",
        "ShyLU_DDCore",
        "ThyraEpetraAdapters",
        "ThyraEpetraExtAdapters",
        "Triutils",
        "Domi",
        "Moertel",
        "FEI",
        "Komplex",
        "Rythmos",
        "Pike",
        "TriKota"
    ]
    regex_pattern = r"|".join([fr"(\S+ENABLE_{x})[^2]" for x in deprecated_packages])
    regex_pattern_defines = r"|".join([fr"ASSERT_DEFINED.*{x}[^2]*[\s)]" for x in deprecated_packages])
    regex_pattern_global_set = r"|".join([fr"GLOBAL_SET.*{x}[^2]*[\s)]" for x in deprecated_packages])

    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.cmake') or file == 'CMakeLists.txt':
                file_path = os.path.join(root, file)
                remove_if_conditions(file_path, regex_pattern)
                remove_regex(file_path, regex_pattern_defines)
                remove_assert_defineds(file_path, regex_pattern)
                remove_regex(file_path, regex_pattern_global_set)

if __name__ == "__main__":
    # Example usage
    directory_to_process = os.getcwd()
    process_directory(directory_to_process)
