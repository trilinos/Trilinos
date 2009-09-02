#!/usr/bin/env python

## Script which compares orderings provided by Scotch and ParMetis by using
## zdrive on MatrixMarket files.

import re, sys, getopt, commands, time, os

# Paths to executables and names
# User can modify this

executable="../../Obj_linux64scotch/zdrive"
zdriveinp={'Scotch': "zdrive.inp.scotch" , 'ParMetis':   "zdrive.inp.parmetis" }
evaltool='../../Obj_linux64scotch/order_eval'
orderfile='./toscotchperm.py'
gunzip='gunzip'
logdir="logs"
resultdir="results"
currentfilename = "current"

# Displays help

def usage():
    """Print the usage information of the command"""
    print "Usage: ./ordering [options] inputfile"
    print "Options:"
    print "\t--force (-f):   force computations even if the number of processes is not a power of two"
    print "\t--np (-n):      number of processes used to produce input"
    print "\t--range (-r):   minproc-maxproc"
    print "\t--verbose (-v): some not usefull information"
    print "\t--debug (-d):   keep temporary files produced by Zoltan"
    print "\t--help (-h):    This page"


# Extracts the "value" in the string. Used to parse order_eval outputs.

def extractvalue(string, value):
    return  float(re.split("=", re.search("%s=.*" % value, string).group())[1])

def opcextract(filename, procnbr, algorithm):
    (status, output) = commands.getstatusoutput("%s %s.mtx %s/%s-%d.%s " % (evaltool, currentfilename, resultdir, filename, procnbr, algorithm))
    if (status == 0):
        return {'%sNNZ' % algorithm: extractvalue(output, 'NNZ'),
        '%sOPC' % algorithm: extractvalue(output, 'OPC')}
    else:
        return {'%sNNZ' % algorithm: 0, '%sOPC' % algorithm: 0}


# Extracts timing results from parmetis output

def extracttime(string):
    val = re.search("Partitioner Library time[\s:]*Max[\s:]*([\d\.e]*)", string)
    if (val == None):
        return 0
    return float(val.group(1))


# Computes the ordering and launchs evaluation of quality.

def computeordering(filename, procnbr, verbose):
    results = {}

    for key in zdriveinp.keys():
        if (verbose == True):
            print "Running ordering with %s ..." % key
        output = commands.getoutput("mpirun -np %d %s %s | tee %s/%s-%s-%d.%s" % (procnbr, executable, zdriveinp[key], logdir, key, filename, procnbr, time.strftime("%Y%m%d%H%m%S")))
        results.update({"%sTIME" % key: extracttime(output)})
        output = commands.getoutput("%s -n %d -o %s/%s-%d.%s %s" % (orderfile, procnbr, resultdir, filename, procnbr, key, currentfilename))

    if (verbose == True):
        print "Computing OPC ..."

    for key in zdriveinp.keys():
        results.update(opcextract(filename, procnbr, key))

    return results

# Prints the results

def displayresults(procmin, procmax, results, values):
    p = procmin
    print " --- %s ----" % values
    for key in zdriveinp.keys():
        print "\t%s\t" % key,
    print '\n',
    while p<= procmax:
        print "%d"%p,
        for key in zdriveinp.keys():
            print "\t%e" %results[p]['%s%s' % (key, values)],
        print "\n",
        p *=2


# Removes temporary files

def docleaning(filename, procnbr):
    digits = len("%d"%procnbr)
    for p in range(0,procnbr):
        file = "%s.out.%d.%s" % (filename, procnbr, ("%d"%p).zfill(digits))
        if os.path.exists(file):
            os.remove(file)


# Says if a number is a power of two

def ispowerof2(num):
    p=1
    while (p < num):
        p *= 2
    return (p == num)


# Verifies that we have all executable files

def doverify(proc, force):
    if (os.path.exists(executable) == False):
            print "You have to compile %s" % executable
            sys.exit(2)
    for key in zdriveinp.keys():
        if (os.path.exists(zdriveinp[key]) == False):
            print "Zdrive.inp file not found : %s" % zdriveinp[key]
            sys.exit(2)
    if (os.path.exists(evaltool) == False):
        print "You have to compile %s" % evaltool
        sys.exit(2)
    os.chmod(orderfile, 0755)
    for dir in (logdir, resultdir):
        if (not os.path.isdir(dir)):
            os.mkdir(dir)
    if ((not force) and (not ispowerof2(proc))):
        print "ParMetis needs a not null power of two number of processors to work."
        print "You can force run by using --force option."
        sys.exit(2)

# The main routine.

def main(argv=None):
    if argv is None:
        argv = sys.argv
    procmin = 2
    procmax=2

    try:
        opts, args = getopt.gnu_getopt(argv[1:], "hi:vn:r:df", ["help", "input=", "verbose", "np=", "range=", "debug", "force"])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    if (len(args) != 1):
        # print help information and exit:
        usage()
        sys.exit(2)

    pathname = args[0]
    verbose = False
    debug = False
    force = False
    outfile = sys.stdout
    offset = -1
    for o, a in opts:
        if o in ("-v", "--verbose"):
            verbose = True
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-n", "--np"):
            procmin = procmax = int(a)
        if o in ("-r", "--range"):
            if (not re.match("[0-9]+-[0-9]+", a)):
                print "Range of processors has to be a range : 'min-max'"
                sys.exit(2)
            (procmin, procmax) = re.split("-", a)
            procmin = int(procmin)
            procmax = int(procmax)
        if o in ("-d", "--debug"):
            debug = True
        if o in ("-f", "--force"):
            force = True
    # ...
    doverify(procmin, force)

    if os.path.exists("%s.mtx"% currentfilename):
        os.remove("%s.mtx" % currentfilename)

    if (re.search("\.gz$", pathname)):
        # Uncompress gzip file
        if (verbose == True):
            print "Decompressing graph file ..."
        commands.getoutput("%s -c %s > %s.mtx"% (gunzip, pathname, currentfilename))
    else:
        # Only do a symlink
        if (verbose == True):
            print "Creating symbolic link to graph file ..."
        os.symlink(pathname, "%s.mtx" % currentfilename)

    p = procmin
    filename = os.path.basename(pathname)
    results = {}
    while p <= procmax:
        if (verbose):
            print "Now computing with %d processors" % p
        results.update({p: computeordering (filename, p, verbose)})
        if (debug == False):
            docleaning(currentfilename, p)
        p *= 2

    if (not debug):
        os.remove("%s.mtx" % currentfilename)

    displayresults(procmin, procmax, results, "NNZ")
    displayresults(procmin, procmax, results, "OPC")
    displayresults(procmin, procmax, results, "TIME")

if __name__ == "__main__":
    main()

