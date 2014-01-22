#!/usr/bin/env python
import subprocess
import commands
import optparse
import math
import os
import shutil
import glob
import re

# ========================= constants =========================
CPN          = 24                                               # cores per node
NUM_RUNS     = 2                                                # number of same runs (for reproducibility)
BASECASE     = 1                                                # number of nodes in the smallest run
DIR_PREFIX   = 'run_'                                           # directory prefix for runs (run with 2 nodes is stores in ${DIR_PREFIX}_2 (shell notation)
PETRA_PREFIX = 'petra-'                                         # actual PBS script prefix
LABELS       = ['setup', 'solve', 'total']                      # analysis header labels
TIMELINES    = ['MueLu Setup', 'Belos Solve', 'Global Time']    # analysis search strings

# ========================= parser =========================
def controller():
    p = optparse.OptionParser()

    # action arguments
    p.add_option('-a', '--analyze',    dest="action", action="store_const", const='analyze', default="none")
    p.add_option('-b', '--build',      dest="action", action="store_const", const='build')
    p.add_option('-c', '--clean',      dest="action", action="store_const", const='clean')
    p.add_option('-r', '--run',        dest="action", action="store_const", const='run')

    # env arguments
    p.add_option('-e', '--exec',       dest="binary",       default="MueLu_ScalingTestParamList.exe")   # executable
    p.add_option('-o', '--output',     dest="output",       default="screen")                           # output files for analysis
    p.add_option('-p', '--petra',      dest="petra",        default='both')
    p.add_option('-s',                 dest="nscale",       default="8", type='int')                    # number of weak scaling runs
    p.add_option('-t', '--template',   dest="template",     default="petra.pbs.template")               # template pbs file for all runs

    # run arguments
    p.add_option(      "--cmds",       dest="cmds",         default="")
    p.add_option('-m', '--matrixType', dest="matrix",       default="Laplace3D")
    p.add_option('-n', '--nx',         dest="nx",           default="134", type="int")
    p.add_option('-x', '--xml',        dest="xmlfile",      default="")

    # parse
    options, arguments = p.parse_args()

    if   options.petra == 'epetra': petra = 1
    elif options.petra == 'tpetra': petra = 2
    elif options.petra == 'both'  : petra = 3
    else:
        print("Unknown petra type %s" % options.petra)
        return

    if options.action == 'build':
        # validate options
        if   options.matrix == "Laplace3D" or options.matrix == "Elasticity3D":
            dim = 3
        elif options.matrix == "Laplace2D" or options.matrix == "Elasticity2D":
            dim = 2
        else:
            print("Uknown matrix type %s" % options.matrix)
            return

        cmds      = []
        datafiles = []
        if options.cmds == "":
            if options.xmlfile == "":
                print("Please provide at least one of:\n"
                      " - xmlfile           ['-x'/'--xml']\n"
                      " - command arguments ['--cmds']")
            else:
                datafiles.append(options.xmlfile)
                cmds.append("--xml=" + options.xmlfile)

        else: # options.cmds != ""
            cmds = re.split(',', options.cmds)

            for i in range(0, len(cmds)):
                xmlfile = re.findall('(?<=--xml=)[^\s]*', cmds[i])
                if len(xmlfile) >= 1:
                  xmlfile = xmlfile[0]
                else:
                  xmlfile = ""

                if xmlfile == "":
                    if options.xmlfile != "":
                        xmlfile = options.xmlfile
                        cmds[i] = cmds[i] + " --xml=" + xmlfile
                    else:
                        print("Every command requires a '--xml' option. However, command '" + cmds[i] +
                              "' does not have one, and no default was provided through 'x'/'--xml'.")
                        return
                else:
                    if options.xmlfile != "" and xmlfile != options.xmlfile:
                        print("WARNING: command '" + cmds[i] + "' provides an xmlfile '" + xmlfile +
                              "' overriding the provided through '--xml' option ['" + options.xmlfile + "']")

            datafiles.append(xmlfile)

        # main loop [weak scaling study]
        # start with a problem on one node and problem size NxN (or NxNxN), then increase nodes quadratically (cubically)
        # and each dimension linearly (so that the problem also increases quadratically (cubically))
        for i in range(1, options.nscale+1):
            nnodes = i**dim                                 # number of nodes
            nx     = i * options.nx                         # single dimension of the problem
            build(nnodes=nnodes, nx=nx, binary=options.binary, petra=petra, matrix=options.matrix,
                datafiles=datafiles, cmds=cmds, template=options.template, output=options.output)

    elif options.action == 'run':
        run()

    elif options.action == 'clean':
        yesno = raw_input("Are you sure [yes]? ")
        if (yesno == "yes"):
            clean()

    elif options.action == 'analyze':
        r = analyze(petra, analysis_run=options.output)
        if r : print(r)

    else:
        print("You need to specify at least one action option")
        return


# ========================= main functions =========================
def analyze(petra, analysis_run):
    # test which of [et]petra is being run
    has_epetra = (len(glob.glob(DIR_PREFIX + "**/*.epetra")) > 0) and (petra & 1)
    has_tpetra = (len(glob.glob(DIR_PREFIX + "**/*.tpetra")) > 0) and (petra & 2)

    if has_epetra == False and has_tpetra == False:
        return "Cannot find any of *.[et]petra files"

    # construct header
    analysis_run_string = "Analysis is performed for " + analysis_run
    if   (petra == 3):
      analysis_run_string += ".[et]petra"
    elif (petra == 1):
      analysis_run_string += ".epetra"
    elif (petra == 2):
      analysis_run_string += ".tpetra"

    print(analysis_run_string)
    header = "                    :"
    for name in LABELS:
        if has_epetra:
            header = header + "  " + name + "-etime      eff"
        if has_tpetra:
            header = header + "  " + name + "-ttime      eff"
    print(header)

    # initialize lists
    time_epetra     = list2dict(TIMELINES)
    eff_epetra      = list2dict(TIMELINES)
    time_tpetra     = list2dict(TIMELINES)
    eff_tpetra      = list2dict(TIMELINES)
    basetime_epetra = list2dict(TIMELINES)
    basetime_tpetra = list2dict(TIMELINES)

    for dir in sort_nicely(glob.glob(DIR_PREFIX + "*")):
        os.chdir(dir)

        nnodes = int(dir.replace(DIR_PREFIX, ''))

        fullstr = "%20s:" % dir

        # test if there is anything to analyze
        if len(glob.glob("screen.out.*.*")) == 0:
            if (has_epetra and os.path.exists(analysis_run + ".epetra")) or (has_tpetra and os.path.exists(analysis_run + ".tpetra")):
                fullstr += " running now?"
            else:
                fullstr += " not run"
            print(fullstr)
            os.chdir("..")
            continue

        for s in TIMELINES:
            if has_epetra:
                r = commands.getstatusoutput("grep \"" + s + "\" " + analysis_run + ".epetra | cut -f3 -d')' | cut -f1 -d'('")
                if r[0] != 0:
                    return "Error reading \"" + analysis_run + ".epetra"
                time_epetra[s] = float(r[1])
                if nnodes == BASECASE:
                    basetime_epetra[s] = time_epetra[s]
                eff_epetra[s] = 100 * basetime_epetra[s] / time_epetra[s]
                fullstr += "%13.2f %7.2f%%" % (time_epetra[s], eff_epetra[s])
            if has_tpetra:
                r = commands.getstatusoutput("grep \"" + s + "\" " + analysis_run + ".tpetra | cut -f3 -d')' | cut -f1 -d'('")
                if r[0] != 0:
                    return "Error reading \"" + analysis_run + ".tpetra"
                time_tpetra[s] = float(r[1])
                if nnodes == BASECASE:
                    basetime_tpetra[s] = time_tpetra[s]
                eff_tpetra[s] = 100 * basetime_tpetra[s] / time_tpetra[s]
                fullstr += "%13.2f %7.2f%%" % (time_tpetra[s], eff_tpetra[s])

        print(fullstr)

        os.chdir("..")

def build(nnodes, nx, binary, petra, matrix, datafiles, cmds, template, output):
    dir = DIR_PREFIX + str(nnodes)
    print("Building %s..." % dir)

    ensure_dir(dir)                                         # make directory if necessary
    os.symlink("../" + binary, dir + "/binary.exe")         # link binary file
    for afile in datafiles:                                 # copy all xml files
        shutil.copy(afile, dir)

    # construct PBS script from template
    os_cmd = "cat " + template
    os_cmd += (" | sed \"s/_NODES_/"    + str(nnodes)               + "/g\"" +
               " | sed \"s/_CORES_/"    + str(nnodes*CPN)           + "/g\"" +
               " | sed \"s/_MTYPE_/"    + matrix                    + "/g\"" +
               " | sed \"s/_NX_/"       + str(nx)                   + "/g\"" +
               " | sed \"s/_EPETRA_/"   + str(petra & 1)            + "/g\"" +
               " | sed \"s/_TPETRA_/"   + str(petra & 2)            + "/g\"" +
               " | sed \"s/_NUM_RUNS_/" + str(NUM_RUNS)             + "/g\"" +
               " | sed \"s/_NUM_CMDS_/" + str(len(cmds))            + "/g\"")
    for i in range(len(cmds)):
        os_cmd += " | sed \"s/_CMD" + str(i+1) + "_/" + cmds[i] + "/g\""
    os_cmd += " > " + dir + "/" + PETRA_PREFIX + str(nnodes) + ".pbs"
    os.system(os_cmd)

def clean():
    for dir in sort_nicely(glob.glob(DIR_PREFIX + "*")):
        print("Cleaning %s..." % dir)
        shutil.rmtree(dir)

def run():
    for dir in sort_nicely(glob.glob(DIR_PREFIX + "*")):
        print("Running %s..." % dir)

        os.chdir(dir)
        os.system("qsub " + PETRA_PREFIX + dir.replace(DIR_PREFIX, '') + ".pbs")
        os.chdir("..")

# ========================= utility functions =========================
def sort_nicely(l):
    """ Sort the given list in the way that humans expect. """
    convert      = lambda s: int(s) if s.isdigit() else s
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] # turn a string into a list of string and number chunks ("z23a" -> ["z", 23, "a"])

    l.sort(key=alphanum_key)
    return l

def ensure_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)
def list2dict(l):
    return dict(zip(l, [0]*len(l)))
# ========================= main =========================
def main():
    controller()

if __name__ == '__main__':
    main()
