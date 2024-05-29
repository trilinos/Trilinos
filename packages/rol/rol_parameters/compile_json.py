import re
import pathlib
import json
from collections import OrderedDict
from find_files import find_files
from read_cpp_source import read_cpp_source

def compile_json( pattern            : re.Pattern,
                  root_dir           : pathlib.Path,
                  relative_pathfiles : list[pathlib.Path],
                  num_capture_groups : int = 1 ) -> str:

    all_instances = OrderedDict()

    for relative_pathfile in relative_pathfiles:
        cpp = read_cpp_source(root_dir / relative_pathfile)
        matches = list(re.finditer(pattern, cpp))
        file_str = str(relative_pathfile)

        if len(matches):
            for m in matches:
                key_name = m.group(1).strip()
                if key_name not in all_instances.keys():
                    all_instances[key_name] = {file_str}
                else:
                    all_instances[key_name].add(file_str)

    for k,v in all_instances.items():
        all_instances[k] = list(v)

    return json.dumps(all_instances,indent=4)

