##---------------------------------------------------------------------------##
## Make a new file from a template
##---------------------------------------------------------------------------##

import sys, os, re, time
        
##---------------------------------------------------------------------------##
## REGULAR EXPRESSIONS
##---------------------------------------------------------------------------##

pkg_s   = re.compile("\<pkg\>")
spkg_s  = re.compile("\<spkg\>")
tpkg_s  = re.compile("\<tpkg\>")

basename_s  = re.compile("\<basename\>")
namespace_s = re.compile("\<namespace\>")
classname_s = re.compile("\<class\>")

user_s = re.compile("\<user\>")
date_s = re.compile("\<date\>")
year_s = re.compile("\<year\>")

plus_s  = re.compile("\+")
minus_s = re.compile("\-")

start_s = re.compile("\<start\>")    

##---------------------------------------------------------------------------##
## TAGS and Environment
##---------------------------------------------------------------------------##

class Nemesis_Env:
    nemesis_env  = ""
    template_dir = ""
    current      = ""

    def __init__(self):
        # location of nemesis development environment
        self.nemesis_env = "NEMESIS_ENVIRONMENT_DIR"
        
        # nemesis templates
        self.template_dir = "%s/templates" % (self.nemesis_env)

        # current directory
        self.current = os.getcwd()

##---------------------------------------------------------------------------##

class Tags:
    pkg  = ""
    tpkg = ""
    spkg = ""
    
    basename  = ""
    namespace = ""
    classname = ""
    
    user = ""
    date = ""
    year = ""

    start = ""

    def __init__(self, environment):
        # user name
        self.user = os.environ["LOGNAME"]
        
        # time/date
        self.date = time.strftime("%a %b %d %H:%M:%S %Y", time.localtime())
        self.year = time.strftime("%Y", time.localtime())

        # package
        if os.path.basename(environment.current) == "test":
            self.pkg  = os.path.basename(os.path.dirname(environment.current))
            self.tpkg = "%s/test" % (self.pkg) 
        else:
            self.pkg = os.path.basename(environment.current)

        # make safepackage name
        (self.spkg, count) = plus_s.subn("p", self.pkg)
        (self.spkg, count) = minus_s.subn("m", self.spkg)

        # default namespace to safepackage name
        self.namespace = self.spkg

    def make_file(self, environment, template,  output):
        
        # exit if file exists
        if os.path.isfile(output):
            print ">>> file %s exists, skipping" % (output)
            return

        # template file
        template_file = "%s/%s" % (environment.template_dir, template)
        f = open(template_file, 'r')
        lines = f.read()
        f.close()
    
        # make substitutions
        (ns, count) = pkg_s.subn(self.pkg, lines)
        (ns, count) = spkg_s.subn(self.spkg, ns)
        (ns, count) = tpkg_s.subn(self.tpkg, ns)

        (ns, count) = basename_s.subn(self.basename, ns)
        (ns, count) = namespace_s.subn(self.namespace, ns)
        (ns, count) = classname_s.subn(self.classname, ns)

        (ns, count) = user_s.subn(self.user, ns)
        (ns, count) = date_s.subn(self.date, ns)
        (ns, count) = year_s.subn(self.year, ns)

        (ns, count) = start_s.subn(self.start, ns)

        # make output file
        ofile = open(output, 'w')
        ofile.write(ns)
        ofile.close()

##---------------------------------------------------------------------------##
## 
##---------------------------------------------------------------------------##

