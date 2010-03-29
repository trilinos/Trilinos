#!/usr/bin/env python

## Script which compares orderings provided by Scotch and ParMetis by using
## zdrive on MatrixMarket files.

import re, sys, getopt, commands, time, os
from operator import itemgetter
import numpy as np
import matplotlib.pyplot as plt


# Paths to executables and names
# User can modify this

# Displays help

def usage():
    """Print the usage information of the command"""
    print "Usage: ./ordering [options] inputfile"
    print "Options:"
    print "\t--help (-h):    This page"



def displayresults(result, how):
    prev = {'name':"",'method':""}
    for r in sorted(result,key=itemgetter('name', 'method', 'threads')):
        if prev['name'] != r['name'] or prev['method'] != r['method']:
            prev['name']=r['name']
            prev['method']=r['method']
            print "\n",prev['name'],"\t", prev['method'], "\t",
        for key in how:
            print r[key],"\t",


def makechart(result, what):
    chartdata={}
    prev = {'name':"",'method':""}
    ticks=[]
    for r in sorted(result,key=itemgetter('name', 'method', 'threads')):
        if prev['name'] != r['name'] or prev['method'] != r['method']:
            prev['name']=r['name']
            prev['method']=r['method']
            chartdata[prev['method']]=[]
        chartdata[prev['method']].append(r[what])
        ticks.append(r['threads'])
    print chartdata
    N=5
    ind = np.arange(N)
    width = 0.2
    ax = plt.figure()
    rects={}
    count=0
    colors=('r','g','b','y','c','m')
    keys=('MMD', 'Hund-1000', 'Metis', 'Scotch')
    for key in keys:
        rects[key] = plt.bar(ind + count*width, map(lambda x:float(x), chartdata[key]), width, color=colors[count])
        count+=1
    plt.ylabel(what)
    plt.title('Results for %s' % prev['name'])
    plt.xticks(ind+width/2*count, ticks)
    plt.savefig('graphs/%s-%s.png' % (prev['name'], what))



def parseregion(filename):
    result=[]
    begining = re.compile("MATRIX")
    collection = (re.compile(r"MATRIX (?P<name>[a-zA-Z0-9_\-]*) with (?P<method>[a-zA-Z0-9_\-]*) on (?P<threads>\d+) threads", re.MULTILINE),
                  re.compile(r"Factor time  =\s*(?P<ftime>[0-9\.]*)", re.MULTILINE),
                  re.compile(r"Factor flops =\s*(?P<fflop>[0-9\.e\+]*)", re.MULTILINE),
                  re.compile(r"\#NZ in factor L = (?P<nzL>\d+)", re.MULTILINE),
                  re.compile(r"\#NZ in factor U = (?P<nzU>\d+)", re.MULTILINE),
                  re.compile(r"total MB needed (?P<memory>[0-9\.]*)", re.MULTILINE)
                  )

    currentblock =""
    infile = open(filename)
    for line in infile:
        if begining.match(line) != None:
            if currentblock == "":
                currentblock += line
                continue
            currentresult = {}
            for matchobj in collection:
                try:
                    currentresult.update(matchobj.search(currentblock).groupdict())
                except AttributeError:
                    print "Cannot find",matchobj.pattern
            if 'name' in currentresult:
                for key in ('ftime', 'fflop', 'nzL', 'nzU', 'memory'):
                    if key not in currentresult:
                        currentresult.update({key:"0"})
                result.append(currentresult)
            currentblock = line
        if currentblock != "":
            currentblock += line
    for key in ('ftime', 'fflop', 'nzL', 'nzU', 'memory'):
        makechart(result, key)
#    displayresults(result, ('fflop', 'ftime', 'nzL', 'nzU', 'memory'))
#    print result


# The main routine.

def main(argv=None):
    if argv is None:
        argv = sys.argv

    try:
        opts, args = getopt.gnu_getopt(argv[1:], "h", ["help"])
    except getopt.GetoptError:
        # print help information and exit:
        usage()
        sys.exit(2)

    if (len(args) != 1):
        # print help information and exit:
        usage()
        sys.exit(2)

    pathname = args[0]
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
    # ...

    parseregion(argv[1])


if __name__ == "__main__":
    main()

