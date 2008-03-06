#!/usr/bin/env python

#Script utilise pour transformer la sortie de Zoltan pour l'ordering
#en sortie compatible par scotch.

#Script which translate Zoltan ordering output files into Scotch format


import re, sys, getopt

def computeoffset(infile):
    for line in infile:
        currline= re.split(' |\t',line)
        if (currline[0].isalpha()):
            continue
        return int(currline[0])


def dotranslation(infile, outfile, offset):
    """Translate Zoltan ordering format into Scotch format"""
    counteroffset = 1 - offset
    for line in infile:
        currline= re.split(' |\t',line)
        if (currline[0].isalpha()):
            continue
        buff = buff = "%d\t%d\n" % (int(currline[0])-counteroffset, int(currline[2])+offset)
        outfile.write(buff)

def computevertices(file):
    size=100
    file.seek(-100, 2)
    for line in file:
        currline= re.split(' |\t',line)
    return int(currline[0])

def usage():
    """Print the usage information of the command"""
    print "Usage: ./toscotch [-n procs] [-o outputfile] [--offset offset] [-v] [-h] inputfile"
    print "Options:"
    print "\t--np (-n):      number of processes used to produce input"
    print "\t--ouput (-o):   outputfile in Scotch ordering format"
    print "\t--offset:       displacement to add to Zoltan results"
    print "\t--verbose (-v): some not usefull information"
    print "\t--help (-h):    This page"


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        opts, args = getopt.gnu_getopt(argv[1:], "hi:o:vn:", ["help", "output=", "input=", "verbose", "np=", "offset="])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    if (len(args) != 1):
        # print help information and exit:
        usage()
        sys.exit(2)

    filename = args[0]
    verbose = False
    outfile = sys.stdout
    procnbr = 1
    offset = -1
    for o, a in opts:
        if o in ("-v", "--verbose"):
            verbose = True
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        if o in ("-o", "--output"):
            outfile = open(a, "w+")
        if o in ("-n", "--np"):
            procnbr = int(a)
        if o == "--offset":
            offset = int(a)
    # ...
    digits = len("%d" % (procnbr-1))

    if (offset <0):
        if (verbose == True):
            print "Computing offset ... "
        file = open("%s.out.%d.%s" % (filename, procnbr, "0".zfill(digits)))
        offset = computeoffset(file)
        file.close()

    if (verbose == True):
        print "Computing number of vertices ..."

    file = open("%s.out.%s.%d" % (filename, ("%d" %procnbr).zfill(digits), procnbr - 1))
    vertices = computevertices(file)
    outfile.write("%d\n" % vertices)

    if (verbose == True):
        print "Processing with :"
        print "input = %s.out.%d.[0-%d]" % (filename, procnbr, procnbr-1)
        print "offset = %d, vertices = %d" % (offset, vertices)

    for proc in range(0, procnbr):
        file = "%s.out.%d.%s" % (filename, procnbr, ("%d"%proc).zfill(digits))
        if verbose == True:
            print "Opening %s" % file
        infile=open(file)
        dotranslation(infile, outfile, offset)
        infile.close()
    outfile.close()

if __name__ == "__main__":
    main()

