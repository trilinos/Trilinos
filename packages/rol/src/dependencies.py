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

def get_headers(cwd,filename,headers=list()):
   new_headers = includes(get_pathfile(cwd,filename))
   if all([ h in headers for h in new_headers ]):
       return headers+new_headers
   else:
       x = list()
       [x.append(get_headers(cwd,h,headers+new_headers)) for h in new_headers if h not in headers]
       return x

def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
            yield el    

if __name__ == '__main__':

    filename = sys.argv[1]
    dep = dict()
    
    cwd = os.getcwd()

    headers = get_headers(cwd,filename)
    headers = [ h for h in flatten(headers) ]
    counts = Counter(headers)

    msg = "Checking dependencies for file {0}\n\n".format(filename)
    print(msg)
    

    width = max( len(k) for k in counts.keys() )
    row = "{0:<" + "{0}".format(width+2) + "} : {" + "1}"
    row.format("File Name","Times Opened")
    rows = [row.format(k,v) for k,v in counts.items()]
    
    print('-' * max([len(r) for r in rows]))
    for r in rows:
        print(r)
    




