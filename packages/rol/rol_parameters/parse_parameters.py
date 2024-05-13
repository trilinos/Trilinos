import re

def crop_to_scope( cpp, token ):
    pattern = re.compile(rf'\{{[^{{}}]*{token}[^{{}}][^{{}}]*\}}', re.MULTILINE)
    match = re.search(pattern, cpp)
    if match:
        return match.group()

def get_sublist_variable_name(cpp,token):
    pattern = re.compile(rf'\{{[^{{}}]*{token}[^{{}}][^{{}}]*\}}', re.MULTILINE)

