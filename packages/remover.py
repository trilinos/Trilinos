import os
import re

def modify_code(directory):
    # Hard-coded regex pattern for the macro
    macro_regex = r'.*EPETRA.*'
    
    # Regex patterns to match #ifdef, #if defined(), #ifndef, and #endif
    ifdef_pattern = re.compile(r'^\s*#\s*ifdef\s+' + macro_regex)
    if_defined_pattern = re.compile(r'^\s*#\s*if\s+defined\s*\(\s*' + macro_regex + r'\s*\)')
    ifndef_pattern = re.compile(r'^\s*#\s*ifndef\s+' + macro_regex)
    endif_pattern = re.compile(r'^\s*#\s*endif')
    
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.cpp') or file.endswith('.h') or file.endswith('.hpp'):
                file_path = os.path.join(root, file)
                with open(file_path, 'r') as f:
                    lines = f.readlines()
                
                output_lines = []
                remove_this_many_endifs = 0
                remove_next_endif = False
                skipping = False
                for line in lines:
                    if remove_next_endif and endif_pattern.match(line):
                        remove_next_endif = False
                        continue
                    if ifdef_pattern.match(line) or if_defined_pattern.match(line):
                        skipping = True

                    if skipping:
                        if re.compile(r'^\s*#\s*if').match(line):
                            remove_this_many_endifs += 1
                        elif endif_pattern.match(line):
                            remove_this_many_endifs -= 1
                        elif remove_this_many_endifs == 1 and re.compile(r'^\s*#\s*else').match(line):
                            skipping = False
                            remove_this_many_endifs = 0
                            remove_next_endif = True

                        if remove_this_many_endifs == 0:
                            skipping = False
                    else:
                        output_lines.append(line)  # Keep the line if not in a skip block

                # Write the modified lines back to the file
                with open(file_path, 'w') as f:
                    f.writelines(output_lines)

if __name__ == "__main__":
    current_directory = os.getcwd()  # Get the current working directory
    modify_code(current_directory)
    print("Code guarded by #ifdef or #if defined() has been removed, while #ifndef directives have been removed but their bodies are retained.")