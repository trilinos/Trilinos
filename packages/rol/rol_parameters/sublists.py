import re
from find_files import find_files
from read_cpp_source import read_cpp_source


def find_sublist_instances(root_path):
    assert(root_path.exists())
    assert(root_path.is_dir())

    files = [find_files(root_path, "sublist", include=['*.hpp'])
    sublist_pattern = re.compile(r'^.*?(\.\s*sublist\s*\(\s*"[^"]*"\s*\)(.*);)',re.MULTILINE)
#    sublist_pattern = re.compile(r'^.*?(\.\s*sublist\s*\(\s*"[^"]*"\s*\))',re.MULTILINE)

    results = dict()

    for file in files:
        cpp = read_cpp_source(file)
        matches = list(re.finditer(sublist_pattern, cpp))
        if len(matches):
            if file not in results.keys():
                results[file] = []
            [results[file].append(m.groups()) for m in matches]

    return results




