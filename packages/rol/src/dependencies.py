import os
import sys
from collections.abc import Iterable
from collections import Counter

def get_pathfile(cwd,filename):
    for root,dirs,files in os.walk(cwd):
        if filename in files:
            pathfile = os.path.join(root,filename)
            return pathfile
    

def includes(pathfile):
    if pathfile is None:
        return list()
    if os.path.exists(pathfile):
        # Read all lines that include a ROL file
        with open(pathfile,'r') as f:
            lines = [line for line in f.readlines() if "#include" in line and "ROL_" in line ]
        quote_includes = [line.split("\"")[1] for line in lines if "\"" in line]
        bracket_includes = [line.split("<")[1].split(">")[0] for line in lines if ">" in line and "<" in line]
        headers = [ h for h in quote_includes if len(h) ] + \
                  [ h for h in bracket_includes if len(h) ]
        return headers
    else:
        return list()

def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el    


def get_headers(cwd,filename,opened_headers=list()):

    # Add the current file name to the list of files that have been opened
    opened_headers.append(filename)

    # Get the name of every file included in the current header provided it has not 
    # previously been opened
    new_headers = [h for h in includes(get_pathfile(cwd,filename)) if h not in opened_headers]

    # If there are any new header files discovered recursively call this function and 
    # flatten the lists 
    if len(new_headers):
        next_level = [get_headers(cwd,h,opened_headers) for h in new_headers]
        return list(flatten(next_level))

    # No new headers to open. Process completed
    else:
        return opened_headers


if __name__ == '__main__':

    filename = sys.argv[1]
    dep = dict()
    
    cwd = os.getcwd()

    headers = get_headers(cwd,filename)
    headers = [ h for h in flatten(headers) ]
    headers.sort()
    counts = Counter(headers)
    msg = "Checking dependencies for file {0}\n\n".format(filename)
    print(msg)
    

    width = max( len(k) for k in counts.keys() )
    row = "{0:<" + "{0}".format(width+2) + "} : {" + "1}"
    print(row.format("File Name","Times Opened"))
    rows = [row.format(k,v) for k,v in counts.items()]
    
    print('-' * max([len(r) for r in rows]))
    for r in rows:
        print(r)
    




