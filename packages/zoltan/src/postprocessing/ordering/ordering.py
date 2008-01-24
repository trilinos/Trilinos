#!/usr/bin/env python

## Script which compares orderings provided by Scotch and ParMetis by using
## zdrive on MatrixMarket files.

import re, sys, getopt, commands, time, os

executable={'Scotch': "../../Obj_linux64scotch/zdrive", 'ParMetis':   "../../Obj_linux64/zdrive" }
evaltool='../../Obj_linux64scotch/order_eval'
orderfile='./toscotchperm.py'

def usage():
    """Print the usage information of the command"""
    print "Usage: ./ordering [options] inputfile"
    print "Options:"
    print "\t--np (-n):      number of processes used to produce input"
    print "\t--range (-r):   minproc-maxproc"
    print "\t--verbose (-v): some not usefull information"
    print "\t--debug (-d):   keep temporary files produced by Zoltan"
    print "\t--help (-h):    This page"

def extractvalue(string, value):
    return  float(re.split("=", re.search("%s=.*" % value, string).group())[1])

def opcextract(filename, procnbr, algorithm):
    (status, output) = commands.getstatusoutput("%s current.mtx results/%s-%d.%s " % (evaltool, filename, procnbr, algorithm))
    if (status == 0):
        return {'%sNNZ' % algorithm: extractvalue(output, 'NNZ'),
        '%sOPC' % algorithm: extractvalue(output, 'OPC')}
    else:
        return {'%sNNZ' % algorithm: 0, '%sOPC' % algorithm: 0}

def computeordering(filename, procnbr, verbose):
    results = {}

    for key in executable.keys():
        if (verbose == True):
            print "Running ordering with %s ..." % key
        output = commands.getoutput("mpirun -np %d %s | tee logs/%s-%s-%d.%s" % (procnbr, executable[key], key, filename, procnbr, time.strftime("%Y%m%d%H%m%S")))
        output = commands.getoutput("%s -n %d -o results/%s-%d.%s current" % (orderfile, procnbr, filename, procnbr, key))

    if (verbose == True):
        print "Computing OPC ..."

    for key in executable.keys():
        results.update(opcextract(filename, procnbr, key))

    return results

def displayresults(procmin, procmax, results, values):
    p = procmin
    print " --- %s ----" % values
    while p<= procmax:
        print "%d"%p,
        for key in executable.keys():
            print "\t%e" %results[p]['%s%s' % (key, values)],
        print "\n",
        p *=2

def docleaning(filename, procnbr):
    digits = len("%d"%procnbr)
    for p in range(0,procnbr):
        file = "%s.out.%d.%s" % (filename, procnbr, ("%d"%p).zfill(digits))
        if os.path.exists(file):
            os.remove(file)


def main(argv=None):
    if argv is None:
        argv = sys.argv

    currentfilename = "current"
    procmin = 2
    procmax=2

    try:
        opts, args = getopt.gnu_getopt(argv[1:], "hi:vn:r:d", ["help", "input=", "verbose", "np=", "range=", "--debug"])
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
            (procmin, procmax) = re.split("-", a)
            procmin = int(procmin)
            procmax = int(procmax)
        if o in ("-d", "--debug"):
            debug = True
    # ...

    if (verbose == True):
        print "Creating symbolic link to graph file ..."

    if os.path.exists("current.mtx"):
        os.remove("current.mtx")
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

    displayresults(procmin, procmax, results, "NNZ")
    displayresults(procmin, procmax, results, "OPC")

if __name__ == "__main__":
    main()

