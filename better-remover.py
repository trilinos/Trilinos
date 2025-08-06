import os
import subprocess
import re

ignore_these_macros = {
    "XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES",
    "XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES",
    "EPETRA_NO_32BIT_GLOBAL_INDICES",
    "EPETRA_NO_64BIT_GLOBAL_INDICES"
}

def is_cpp_source_file(filename):
    return filename.endswith(('.c', '.cpp', '.h', '.hpp', '.C', '.H', '.hpp.in', '.h.in', '.H.in'))

def filter_rigs(s):
    return not (s.endswith('_H') or s.endswith('_HPP'))

def find_macros_to_remove(directory, regex_pattern):
    # Compile the regex pattern
    pattern = re.compile(regex_pattern)
    macros_to_remove = set()

    # Walk through the directory
    for root, dirs, files in os.walk(directory):
        for file in files:
            # Process only C/C++ source files
            if is_cpp_source_file(file):
                file_path = os.path.join(root, file)
                print(f"Scanning file: {file_path}")

                # Read the file and find macros
                with open(file_path, 'r') as f:
                    content = f.read()
                matches = pattern.findall(content)
                macros_to_remove.update([group for match in matches for group in match if group])

    # Using filter() to apply the function to the set
    macros_to_remove = set(filter(filter_rigs, macros_to_remove))
    return macros_to_remove - ignore_these_macros

def run_unifdef(directory, macros_to_remove):
    # Create a temporary file with the macros to undefine
    with open('undef_macros.tmp', 'w') as temp_file:
        for macro in macros_to_remove:
            temp_file.write(f"#undef {macro}\n")

    # Walk through the directory again to apply unifdef
    for root, _, files in os.walk(directory):
        for file in files:
            # Process only C/C++ source files
            if is_cpp_source_file(file):
                file_path = os.path.join(root, file)
                print(f"Processing file: {file_path}")

                # Run unifdef command
                try:
                    command = ['/fgs/sebrown/trilinos/Trilinos/unifdef-2.12/unifdef', '-B', '-o', file_path, '-f', 'undef_macros.tmp', file_path]
                    result = subprocess.run(command, capture_output=True, text=True)

                    # Check for errors
                    if result.returncode != 0:
                        print(f"Error processing {file_path}: {result.stderr}")
                    else:
                        print(f"Successfully processed {file_path}")

                except Exception as e:
                    print(f"An error occurred while processing {file_path}: {e}")


def remove_cmakedefines(directory, macros_to_remove):
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith(('.hpp.in', '.h.in')):
                    file_path = os.path.join(root, file)
                    with open(file_path, 'r') as file:
                        lines = file.readlines()
                    output_lines = []
                    for line in lines:
                        if not (line.strip().startswith("#cmakedefine") and any([line.strip().endswith(x) for x in macros_to_remove])):
                            output_lines.append(line)

                    with open(file_path, 'w') as file:
                        file.writelines(output_lines)


if __name__ == "__main__":
    # Specify the directory and regex pattern
    target_directory = os.getcwd()
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
        "DOMI",
        "MOERTEL",
        "FEI",
        "Komplex",
        "Rythmos",
        "Pike",
        "TriKota"
    ]
    regex = r"|".join([fr"\W+({x.upper()}_\w*)\W+|\W+(\w*_{x.upper()})\W+|\W+(\w*_{x.upper()}_\w*)\W+" for x in deprecated_packages])
    print(regex)
    macros_to_remove = find_macros_to_remove(target_directory, regex)
    # with open("macros.txt", 'r') as inf:
    #     macros_to_remove = inf.read().splitlines()
    if macros_to_remove:
        print(f"Macros to remove: {macros_to_remove}")
        #run_unifdef(target_directory, macros_to_remove)
        remove_cmakedefines(target_directory, macros_to_remove)
    else:
        print("No matching macros found.")
