import os
import re

def remove_words_from_sets(text, words_to_remove):
    # Create a regex pattern to match the SET commands
    pattern = r'SET\(([^)]+)\)'

    def replace_function(match):
        # Extract the current list of packages
        current_packages = match.group(1)
        # Split the packages into a list and remove unwanted words
        packages = current_packages.split()
        filtered_packages = [pkg for pkg in packages if pkg not in words_to_remove]
        # Join the filtered packages back into a string
        return f'SET({" ".join(filtered_packages)})'

    # Use re.sub to replace all SET commands with the filtered lists
    modified_text = re.sub(pattern, replace_function, text)
    return modified_text

def process_file(file_path, words_to_remove):
    # Read the content of the file
    with open(file_path, 'r') as file:
        content = file.read()

    # Remove specified words from the SET commands
    modified_content = remove_words_from_sets(content, words_to_remove)

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.write(modified_content)

def find_and_process_files(directory, words_to_remove):
    # Walk through the directory to find all Dependencies.cmake files
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file == 'Dependencies.cmake':
                file_path = os.path.join(root, file)
                print(f"Processing file: {file_path}")
                process_file(file_path, words_to_remove)


if __name__ == "__main__":
    # Specify the directory to search
    target_directory = os.getcwd()
    # Define the list of words to remove
    words_to_remove = ['Amesos', 'AztecOO', 'Epetra', 'EpetraExt', 'Ifpack', 'Intrepid', 'Isorropia', 'ML', 'NewPackage', 'Pliris', 'PyTrilinos', 'ShyLU_DDCore', 'ThyraEpetraAdapters', 'ThyraEpetraExtAdapters', 'Triutils']

    # Find and process all Dependencies.cmake files
    find_and_process_files(target_directory, words_to_remove)
