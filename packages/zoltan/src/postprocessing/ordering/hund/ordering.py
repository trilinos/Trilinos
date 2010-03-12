#!/usr/bin/env python

## Script which compares orderings provided by Scotch and ParMetis by using
## zdrive on MatrixMarket files.

import re, sys, getopt, commands, time, os

# Paths to executables and names
# User can modify this

basename="./bin/zdrive"
convertscript="./bin/toscotchperm.py"
zdriveinp={'Scotch': "%s zdrive.inp.scotch" % basename , 'Metis':   "%s zdrive.inp.metis" % basename,
           'Hund-1000' : "%s_hund_1000" % (basename) }
#'Hund-1000s' : "%s_hund_1000s" % (basename),
#           'Hund-10000' : "%s_hund_10000" % (basename), 'Hund-10000s' : "%s_hund_10000s" % (basename) }
#           'Hund-100' : "%s_hund_100" % (basename), 'Hund-100s' : "%s_hund_100s" % (basename),

logdir="logs"
currentfilename="current"
resultdir="results"

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

# Computes the ordering and launchs evaluation of quality.

def computeordering(filename, procnum, verbose):
    for key in zdriveinp.keys():
        print "Running ordering with %s ..." % key
        output = commands.getoutput("mpirun -np %d %s | tee %s/%s-%s.%s" % (procnum, zdriveinp[key], logdir, filename, key, time.strftime("%Y%m%d%H%m%S")))
        print output
        output = commands.getoutput("%s  --offset=1 -n %d %s -o %s/%s-%s-%d.ord" % (convertscript, procnum, currentfilename, resultdir, filename, key, procnum))
        print output



# The main routine.

def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        opts, args = getopt.gnu_getopt(argv[1:], "hvdn:", ["help",  "verbose",  "debug", "procs="])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    verbose = True
    procnum = 1
    for o, a in opts:
        if o in ("-v", "--verbose"):
            verbose = True
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-d", "--debug"):
            debug = True
        if o in ("-n", "--procs"):
            procnum=int(a);
    # ...

    for filename in args:
        filename = os.path.basename(filename).replace(".gz", "")
        filename = os.path.basename(filename).replace(".mtx", "")
        print filename
	if os.path.exists("%s.mtx.gz" % filename):
            f_extension="mtx.gz"
        else:
            f_extension="mtx"
        if os.path.exists("%s.%s" % (currentfilename, f_extension)):
            os.remove("%s.%s" % (currentfilename, f_extension))

        print "Creating symbolic link to %s file ..." % filename
        os.symlink("%s.%s" % (filename, f_extension), "%s.%s" % (currentfilename, f_extension))


        computeordering (filename, procnum, verbose)
        os.remove("%s.%s" % (currentfilename, f_extension))


if __name__ == "__main__":
    main()

