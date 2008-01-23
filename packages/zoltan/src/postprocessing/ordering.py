#!/usr/bin/env python

## Script which compares orderings provided by Scotch and ParMetis by using
## zdrive on MatrixMarket files.

import re, sys, getopt, commands, time, os


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
    (status, output) = commands.getstatusoutput("./bin/gotst -M current.mtx results/%s-%d.%s " % (filename, procnbr, algorithm))
    if (status == 0):
        return {'%sNNZ' % algorithm: extractvalue(output, 'NNZ'),
        '%sOPC' % algorithm: extractvalue(output, 'OPC')}
    else:
        return {'%sNNZ' % algorithm: 0, '%sOPC' % algorithm: 0}

def computeordering(filename, procnbr, verbose):
    results = {}
    if (verbose == True):
        print "Running ordering with Scotch ..."

    output = commands.getoutput("mpirun -np %d ./bin/zdrive_scotch | tee logs/scotch-%s-%d.%s" % (procnbr, filename, procnbr, time.strftime("%Y%m%d%H%m%S")))

    output = commands.getoutput("./bin/toscotchperm.py -n %d -o results/%s-%d.scotch current" % (procnbr, filename, procnbr))

    if (verbose == True):
        print "Running ordering with ParMetis ..."

    output = commands.getoutput("mpirun -np %d ./bin/zdrive_parmetis | tee logs/parmetis-%s-%d.%s" % (procnbr, filename, procnbr, time.strftime("%Y%m%d%H%m%S")))

    output = commands.getoutput("./bin/toscotchperm.py -n %d -o results/%s-%d.parmetis current" % (procnbr, filename, procnbr))

    if (verbose == True):
        print "Computing OPC ..."

    results.update(opcextract(filename, procnbr, "scotch"))
    results.update(opcextract(filename, procnbr, "parmetis"))

    return results

def displayresults(procmin, procmax, results, values):
    p = procmin
    print " --- %s ----" % values
    while p<= procmax:
        print "%d\t%e\t%e" %(p, results[p]['scotch%s' % values], results[p]['parmetis%s' % values])
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

