#! /usr/bin/env python

"""
MakefileVariables.py - This script can be run as a shell command or imported as
a module into a python script.  Its purpose is to read a Makefile and obtain a
map of make variables and their values.

Shell usage: MakefileVariables.py [options] [filename]

options:
    -h | --help       Print this message and exit
    -m | --make       Output in Makefile format
    -p | --python     Output as a python dictionary (default)
    -v | --version    Print the version number and exit

Python usage: import MakefileVariables

Available functions:
    removeContinuationLines([string, string, ...]) -> None (combine any strings
                                   that end with '\' with its following string)
    isNotBlankLine(string) -> bool (return True if string is not blank)
    parseExportFile(string) -> dict (open filename and obtain make variable
                                     names and their raw values)
    makeSubstitutions(dict) -> None (Interpret dict as varName/value pairs;
                                     wherever '$(...)' appears in values, make
                                     appropriate substitutions)
    specialNormPath(string) -> string (normalize a path, even if it starts with
                                       -I, -L, etc.)
    uniqiufyList([string1, string2, ...]) -> None (remove duplicate strings from
                                                   the list)
    uniquifyString(string) -> string (remove duplicate substrings from the
                                      string)
    uniquifyDict(dict) -> None (Uniquify each value of the given dictionary)
    processFile(string) -> dict (Given a filename, parse the file for make
                                 variable names and values, make substitutions
                                 wherever '$(...)' is found in the values, and
                                 then uniquify the resulting values)
"""

__version__ = "1.0"
__author__  = "Bill Spotz"
__date__    = "Aug 29 2005"

# Import python modules for command-line options, the operating system, regular
# expressions, and system functions
from   getopt import *
import os
import re
import sys

# Define regular expressions for Makefile assignments, blank line, continuation
# lines, include statements and variable references
assignRE   = re.compile(r"([^=]+)=([^=]+)"  )
blankRE    = re.compile(r"^\s*$"            )
continueRE = re.compile(r"(.*)\\\s*$"       )
includeRE  = re.compile(r"\s*include\s+(.+)")
makeVarRE  = re.compile(r"\$\(([^)]+)\)"    )

#############################################################################

def removeContinuationLines(lines):
    """Given lines, a list of strings, check for the continuation character
    ('\') at the end of each line.  If found, combine the appropriate strings
    into one, leaving blank lines in the list to avoid duplication."""
    for i in range(len(lines)-1):
        line = lines[i]
        match = continueRE.search(line)
        if match:
            lines[i+1] = match.group(1) + lines[i+1]
            lines[i  ] = ""

#############################################################################

def isNotBlankLine(s):
    """Return True if a string is not a blank line"""
    if blankRE.match(s):
        return False
    return True

#############################################################################

def parseExportFile(filename):
    """Open filename, read in the text and return a dictionary of make variable
    names and values.  If an include statement is found, this routine will be
    called recursively."""
    lines = open(filename,"r").readlines()     # Read in the lines of the Makefile
    lines = [s.split("#")[0] for s in lines]   # Remove all comments
    removeContinuationLines(lines)             # Remove continuation lines
    lines = filter(isNotBlankLine, lines)      # Remove all blank lines
    dict = { }
    for line in lines:
        # Process include statements
        match = includeRE.match(line)
        if match:
            tempDict = dict
            tempDict["include files"] = match.group(1).strip()
            makeSubstitutions(tempDict)
            for file in tempDict["include files"].split():
                try:
                    dict.update(parseExportFile(file))
                except IOError:
                    pass
            continue
        # Process assignment statements
        match = assignRE.match(line)
        if match:
            dict[match.group(1).strip()] = match.group(2).strip()
    return dict

#############################################################################

def makeSubstitutions(dict):
    """Loop over the items of a dictionary of variable names and string values.
    If the value contains the substring(s) '$(VARNAME)' and VARNAME is a key in
    the dictionary, then substitute the appropriate value.  If VARNAME is not a
    key in the dictionary, then substitute the null string.  For circular
    substitutions, substitute the null string."""
    active    = dict   # Shallow copy
    completed = { }
    i = 0
    # Dictionary active are the varName/value pairs that are still being
    # actively substitutued.  When its length is zero, we are done
    while len(active) > 0:
        varName = active.keys()[i]
        value   = active[varName]
        done    = False
        pos     = 0
        while not done:
            match = makeVarRE.search(value, pos)
            if match:
                subVarName = match.group(1)
                start      = match.start(1)-2
                pos        = match.end(1)  +1
                if subVarName == varName:  # We have encountered a circular reference
                    subValue = ""          #   effectively delete it
                elif subVarName in completed.keys():
                    subValue = completed[subVarName]
                elif subVarName in active.keys():
                    subValue = active[subVarName]
                else:
                    subValue = ""
                active[varName] = "%s%s%s" % (value[:start], subValue, value[pos:])
            else:
                if pos == 0:
                    # We searched from the beginning and found no $(...), so
                    # this varName/value pair is completed
                    completed[varName] = value
                    del active[varName]
                done = True                 # We are done checking this varName/value pair
                i += 1                      # Go to the next varName/value pair,
                if i >= len(active): i = 0  #   wrapping to the beginning, if necessary

    # Update the given dict
    dict.update(completed)

#############################################################################

def specialNormPath(path):
    """Apply os.path.normpath to argument path, but remove a leading option
    such as '-I' or '-L' before the call and add it back before the result is
    returned."""
    if len(path) <= 2: return path
    start = 0
    if path[0] == "-": start = 2
    return path[:start] + os.path.normpath(path[start:])

#############################################################################

def uniquifyList(list):
    """Remove duplicate items from a list, preserving original order."""
    i = 1
    while i < len(list):
        if list[i] in list[:i]:
            del list[i]
        else:
            i += 1

#############################################################################

def uniquifyString(s, delim=' '):
    """Split a string using the specified delimeter, apply specialNormPath to
    each string, uniquify the resulting list, and return a string that joins the
    unique list with the same delimeter."""
    list = s.split(delim)
    list = [specialNormPath(path.strip()) for path in list]
    uniquifyList(list)
    return delim.join(list)

#############################################################################

def uniquifyDict(dict):
    """Loop over each item in a dictionary of variable names and string values
    and uniquify the string."""
    for key in dict:
        dict[key] = uniquifyString(dict[key])

#############################################################################

def processFile(filename):
    """Open filename, read its contents and parse it for Makefile assignments,
    creating a dictionary of variable names and string values.  Substitute
    variable values when '$(...)' appears in a string value."""
    dict = parseExportFile(filename)
    makeSubstitutions(dict)
    uniquifyDict(dict)
    return dict

#############################################################################

def main():
    """This is the routine that gets called if MakefileVariable.py is invoked
    from the shell.  Process any command-line options that may be given, take
    the first argument as a filename, process it to obtain a dictionary of
    variable name/value pairs, and output the results."""

    # Initialization
    (progDir,progName) = os.path.split(sys.argv[0])
    options      = "hmpv"
    long_options = ["help", "make", "python", "version"]
    outStyle     = "python"

    # Get the options and arguemnts from the command line
    (opts,args) = getopt(sys.argv[1:], options, long_options)

    # Loop over options and implement
    for flag in opts:
        if flag[0] in ("-h","--help"):
            print __doc__
            sys.exit()
        elif flag[0] in ("-m", "--make"):
            outStyle = "make"
        elif flag[0] in ("-p", "--python"):
            outStyle = "python"
        elif flag[0] in ("-v", "--version"):
            print progName, __version__, __date__
            sys.exit()
        else:
            print "Unrecognized flag:", flag[0]
            print __doc__
            sys.exit()

    # Process the filename
    dict = processFile(args[0])

    # Output the variable names and values
    if outStyle == "make":
        keys = dict.keys()
        keys.sort()
        for key in keys:
            print key, "=", dict[key]
    elif outStyle == "python":
        print dict

#############################################################################
# If called from the command line, call main()
#############################################################################

if __name__ == "__main__":
    main()
