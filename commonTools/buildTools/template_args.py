#! /usr/bin/env python
# -*- python -*-
"""
template_args.py - A python script for parsing C++ code and listing all unique
                   arguments for a given C++ template.

"""
__version__ = "1.0"
__author__  = "Bill Spotz"
__date__    = "Oct 31 2006"

# System imports
from   optparse import *
import os
import re
import sys

################################################################################

def removeComments(text):
    """
    Return a copy of text with its C-style and C++-style comments removed.
    """
    cCommentRE   = re.compile(r"\/\*.*?\*\/", re.DOTALL   )
    cppCommentRE = re.compile(r"\/\/.*$",     re.MULTILINE)

    for commentRE in (cCommentRE, cppCommentRE):
        match = True
        while match:
            match = commentRE.search(text)
            if match:
                if match.end() > len(text):
                    text = text[:match.start()]
                else:
                    text = text[:match.start()] + text[match.end():]
    return text

################################################################################

def findBlock(text, pos=0):
    """
    Given the input text (potentially multiline) and an optional pos marking the
    starting position, find an opening delimeter -- either (, [, {, <, single
    quote, or double quote -- and return a tuple of integers indicating the
    character indexes of the text block -- closed with ), ], }, >, single quote,
    or double quote, respectively -- while correctly handling nested blocks.
    """

    # Define delimeter strings
    quote1Delimeter = "'"
    quote2Delimeter = '"'
    openDelimeters  = "\(\[\{<"
    closeDelimeters = "\)\]\}>"

    # Define delimeter regular expressions
    quote1RE = re.compile("([" + quote1Delimeter + "])", re.M)
    quote2RE = re.compile("([" + quote2Delimeter + "])", re.M)
    openRE   = re.compile("([" + openDelimeters  +
                                 quote1Delimeter +
                                 quote2Delimeter + "])", re.M)
    anyRE    = re.compile("([" + openDelimeters  +
                                 quote1Delimeter +
                                 quote2Delimeter +
                                 closeDelimeters + "])", re.M)

    # Find the first opening delimeter
    matchObject = openRE.search(text, pos)
    if not matchObject: return (None, None)

    # Initialize the loop
    stack = [ matchObject.group() ]
    start = matchObject.start()
    pos   = start + 1

    # Find the end of the block
    while matchObject:

        # Determine the active delimeter regular expression
        if   stack[-1] == quote1Delimeter:
            activeRE = quote1RE
        elif stack[-1] == quote2Delimeter:
            activeRE = quote2RE
        else:
            activeRE = anyRE

        # Search for the next delimeter
        matchObject = activeRE.search(text, pos)
        if matchObject:
            delimeter = matchObject.group()
            pos       = matchObject.end()

            # Check for matched delimeters
            if (((stack[-1] == quote1Delimeter) and
                 (delimeter == quote1Delimeter)) or
                ((stack[-1] == quote2Delimeter) and
                 (delimeter == quote2Delimeter)) or
                ((stack[-1] == "("            ) and
                 (delimeter == ")"            )) or
                ((stack[-1] == "["            ) and
                 (delimeter == "]"            )) or
                ((stack[-1] == "{"            ) and
                 (delimeter == "}"            )) or
                ((stack[-1] == "<"            ) and
                 (delimeter == ">"            ))):
                stack.pop()                  # Remove the last element from the list
                if len(stack) == 0:
                    return (start, pos)

            # Process unmatched delimeter
            else:
                if (delimeter in openDelimeters  or
                    delimeter == quote1Delimeter or
                    delimeter == quote2Delimeter   ):
                    stack.append(delimeter)  # Add the delimeter to the stack
                else:
                    raise RuntimeError, "findBlock: mismatched delimeters: " + \
                          stack[-1] + " " + delimeter

    # We made it through all of text without finding the end of the block
    raise RuntimeError, "findBlock: open block: " + join(stack)

################################################################################

def main():

    # Set up the command-line parser object
    prog    = os.path.split(__file__)[1]
    usage   = __doc__ + "usage: " + prog + " [options] template file1 ..."
    version = prog + " " + __version__ + " " + __date__
    parser  = OptionParser(usage=usage,version=version)
    parser.set_defaults(verbose=False)
    parser.add_option("-v", "--verbose", action="store_true" , dest="verbose",
                      help="run in verbose mode")
    parser.add_option("-q", "--quiet"  , action="store_false", dest="verbose",
                      help="run in quiet mode [default]")

    # Get the options and arguments from the command line
    (options,args) = parser.parse_args()

    # Check the arguments
    if len(args) < 2:
        print usage
        sys.exit(-1)
    template  = args[0]
    filenames = args[1:]

    # Generate the regular expression
    templateRE = re.compile(template + "<")

    # Initialize the list of template arguments
    template_args = [ ]

    # Loop over the filenames
    for filename in filenames:
        if options.verbose:
            print "Processing '%s' ..." % (filename,),
        try:
            text = open(filename,"r").read()
        except IOError:
            print
            print "Could not open '%s'." % filename
            sys.exit(-2)

        # Remove the comments from the text
        text = removeComments(text)

        # Loop throught the text, searching for our template
        pos   = 0
        match = True
        while match:
            match = templateRE.search(text,pos)
            if match:
                (start,end) = findBlock(text, match.start())
                arg = text[start+1:end-1].strip()
                if not arg in template_args:
                    template_args.append(arg)
                pos = end + 1
        if options.verbose:
            print "ok"

    # Process the template arguments
    template_args.sort()
    if options.verbose:
        print "\nTemplate arguments for %s< ... >:" % template
    for arg in template_args:
        print arg

################################################################################

if __name__ == "__main__":
    main()
