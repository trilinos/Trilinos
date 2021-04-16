import os
import re
import json

def replace_all(text,subsdict):

    rep = dict((re.escape(k), v) for k, v in subsdict.items()) 
    pattern = re.compile("|".join(rep.keys()))
    return pattern.sub(lambda m: rep[re.escape(m.group(0))], text)
    
def get_cppfiles(path):
    return [ os.path.join(root,f) for root,dirs,files in os.walk(path) \
            for f in files if ".cpp" in f or ".hpp" in f ]

def read_file(pathfile):
    with open(pathfile,"r") as f:
        text = f.read()
    return text

def write_textfile(pathfile,text):
    with open(pathfile,"w") as f:
        f.write(text)

def load_substitutions(jsonfile):
    with open(jsonfile) as f:
        subs = json.load(f)
    return subs

def refactor_tree(path,*subs):
    pathfiles = get_cppfiles(path)
    for f in pathfiles:
        print("Refactoring file {0}".format(f))
        text = read_file(f)
        for sub in subs:
            newtext = replace_all(text,sub)
            text = newtext
        write_textfile(f,text)    


if __name__ == '__main__':

    trilinos_home = os.environ["TRILINOS_HOME"]
    rol_home = os.path.join(trilinos_home,"packages","rol")
    rol_src = os.path.join(rol_home,"src")
    rol_test = os.path.join(rol_home,"test")
    rol_example = os.path.join(rol_home,"example")

    rol_testproblems = os.path.join(rol_home,"src","zoo","testproblems")
    

    replace_filenames    = load_substitutions("rol_filenames.json") 
    replace_testproblems = load_substitutions("rol_testproblems.json") 
    replace_types        = load_substitutions("rol_types.json") 
    
    algo_tests = get_cppfiles(os.path.join(rol_test,"algorithm"))

    refactor_tree(rol_testproblems,replace_testproblems)

#    for f in algo_tests:
#        print("Refactoring file {0}".format(f))
#        text1 = read_file(f)
#        text2 = replace_all(text1,replace_filenames)
#        text3 = replace_all(text2,replace_types)
#        write_textfile(f,text3)
        
