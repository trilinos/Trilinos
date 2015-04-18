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
PLATFORM     = "hopper"
CPN          = 0                                                # number of physical cores per node
SCHEDULER    = ''                                               # HPC platform scheduler
SCHED_HEADER = ''                                               # header for the platform scheduler
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
    p.add_option('-e', '--exec',       dest="binary",       default="MueLu_Driver.exe")                 # executable
    p.add_option('-o', '--output',     dest="output",       default="screen")                           # output files for analysis
    p.add_option('-p', '--petra',      dest="petra",        default="both")                             # petra mode
    p.add_option('-N', '--nnodes',     dest="nnodes",       default="")                                 # custom node numbers
    p.add_option('-s',                 dest="nscale",       default=8, type='int')                      # number of scaling runs
    p.add_option('-t', '--template',   dest="template",     default="sched.template")                   # template file for all runs
    p.add_option('-l', '--labels',     dest="ltmodule",     default="")                                 # labels and timelines
    p.add_option(      '--cpn',        dest="cpn",          default=CPN, type='int')                    # cores per node
    p.add_option(      '--cps',        dest="cps",          default=-1, type='int')                     # cores per socket
    p.add_option('-T', '--type',       dest="type",         default="weak")                             # scaling type (weak | strong)
    # FIXME (29 Sep 2014): unified interface is buggy, disabling
    p.add_option('-u', '--unified',    dest="unified",      action="store_true", default=False)         # by default, try to use unified
    p.add_option('-d', '--default',    dest="unified",      action="store_false")                       #   but sometimes we want to use
                                                                                                        #   the default one, particularly
                                                                                                        #   when we get segfaults and such

    # run arguments
    p.add_option(      '--cmds',       dest="cmds",         default="")                                 # additional args for the command
    p.add_option('-m', '--matrixType', dest="matrix",       default="Laplace3D")                        # matrix
    p.add_option('-n', '--nx',         dest="nx",           default=134, type='int')                    # number of nodes in any direction
    p.add_option('-x', '--xml',        dest="xmlfile",      default="")                                 # xml file with hierarchy configuration

    # parse
    options, arguments = p.parse_args()

    unified = options.unified
    if unified == True:
        raise RuntimeError("Unified interface is buggy, and must be disabled")

    if   options.petra == 'epetra': petra = 1
    elif options.petra == 'tpetra': petra = 2
    elif options.petra == 'both'  : petra = 3
    elif options.petra == 'ml'    : petra = 4
    else:
        print("Unknown petra type %s" % options.petra)
        return

    if options.type != 'weak' and options.type != 'strong':
        raise RuntimeError("Scaling type must be one of (weak|strong)")

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
                return

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

                if unified == True:
                    if cmds[i] != "--xml=" + xmlfile:
                        print("WARNING: command '" + cmds[i] + "' provides extra (to xml) arguments, "
                              "disabling construction of a single unified xml file")
                        unified = False

        nnodes = []         # number of nodes
        nx     = []         # single dimension of the problem
        if options.nnodes == "":
            # main loop [scaling study]
            for i in range(1, options.nscale+1):
                if options.type == "weak":
                    nnodes.append(i**dim)
                else:
                    nnodes.append(i**2)

        else:
            # custom number of nodes
            for i in re.split(',', options.nnodes):
                nnodes.append(int(i))

        cpn = options.cpn
        for i in range(0, len(nnodes)):
            if options.type == "weak":
                # scale problem size with the number of nodes
                nx.append(int(options.nx * pow(nnodes[i] * float(cpn)/CPN, 1./dim)))

            elif options.type == "strong":
                # use same problem size for all nodes
                nx.append(options.nx)

        global NUM_RUNS
        for i in range(0, len(nnodes)):
            if unified == True:
                # Assemble unified xml file
                unified_xml = "unified.xml"
                unified_cmd = ["--xml=" + unified_xml]

                qn = "\\\\n"
                os.system("echo -e '<ParameterList   name=\"Driver\">'" + " > " + unified_xml)
                # Move number of runs into an xml file, and reset the global
                os.system("echo -e '  <Parameter     name=\"number of reruns\" type=\"int\" value=\"" + str(NUM_RUNS) + "\"/>'"+qn + " >> " + unified_xml)
                for k in range(0,len(datafiles)):
                    os.system("echo -e '  <ParameterList name=\"Run" + str(k+1) + "\">'" + " >> " + unified_xml)
                    os.system("echo -e '    <Parameter   name=\"filename\" type=\"string\" value=\"cmd" + str(k+1) + "\"/>'"+ qn + qn + " >> " + unified_xml)
                    os.system("cat " + datafiles[k] + " >> " + unified_xml)
                    os.system("echo -e '  </ParameterList>'" + " >> " + unified_xml)

                os.system("echo -e '</ParameterList>'" + " >> " + unified_xml)

                build(nnodes=nnodes[i], nx=nx[i], binary=options.binary, petra=petra, matrix=options.matrix,
                      datafiles=[unified_xml], cmds=unified_cmd, template=options.template, output=options.output, cpn=cpn, cps=options.cps, unified=True, num_runs=1)

            else:
                build(nnodes=nnodes[i], nx=nx[i], binary=options.binary, petra=petra, matrix=options.matrix,
                      datafiles=datafiles, cmds=cmds, template=options.template, output=options.output, cpn=cpn, cps=options.cps, unified=False, num_runs=NUM_RUNS)

    elif options.action == 'run':
        run()

    elif options.action == 'clean':
        yesno = raw_input("Are you sure [yes]? ")
        if (yesno == "yes"):
            clean()

    elif options.action == 'analyze':
        if (options.ltmodule == ""):
          labels    = LABELS
          timelines = TIMELINES
          parsefunc = ""

        else:
          labels    = import_from(options.ltmodule, "LABELS")
          timelines = import_from(options.ltmodule, "TIMELINES")
          try:
            parsefunc  = import_from(options.ltmodule, "PARSEFUNC")
          except(AttributeError):
            parsefunc = ""

        analysis_runs = re.split(',', options.output)
        r = analyze(petra, analysis_runs=analysis_runs, labels=labels, timelines=timelines, parsefunc=parsefunc, scaling_type=options.type)
        if r : print(r)

    else:
        print("You need to specify at least one action option")
        return


# ========================= main functions =========================
def analyze(petra, analysis_runs, labels, timelines, parsefunc, scaling_type):
    # test which of [et]petra is being run
    has_epetra = (len(glob.glob(DIR_PREFIX + "**/*.epetra")) > 0) and (petra & 1)
    has_tpetra = (len(glob.glob(DIR_PREFIX + "**/*.tpetra")) > 0) and (petra & 2)
    has_ml     = (len(glob.glob(DIR_PREFIX + "**/*.ml")) > 0)     and (petra & 4)

    if has_epetra == False and has_tpetra == False and has_ml == False:
        return "Cannot find any of *.[et]petra files"

    # construct header
    analysis_run_string = "Analysis is performed for " + str(analysis_runs)
    if   (petra == 4):
      analysis_run_string += ".ml"
    elif (petra == 3):
      analysis_run_string += ".[et]petra"
    elif (petra == 1):
      analysis_run_string += ".epetra"
    elif (petra == 2):
      analysis_run_string += ".tpetra"

    print(analysis_run_string)
    header = "                    |"
    for name in labels:
        if has_epetra:
            header = header + "  " + name + "-etime      eff"
        if has_tpetra:
            header = header + "  " + name + "-ttime      eff"
        if has_ml:
            header = header + "  " + name + "-mltime     eff"
    print(header)
    separator = '-' * len(header)

    # initialize lists
    time_epetra     = list2dict(timelines)
    time_ml         = list2dict(timelines)
    time_tpetra     = list2dict(timelines)
    eff_epetra      = list2dict(timelines)
    eff_tpetra      = list2dict(timelines)
    eff_ml          = list2dict(timelines)
    basetime_epetra = list2dict(timelines)
    basetime_tpetra = list2dict(timelines)
    basetime_ml     = list2dict(timelines)

    for analysis_run in analysis_runs:
        print(separator)

        for dir in sort_nicely(glob.glob(DIR_PREFIX + "*")):
            if os.path.isdir(dir):
                os.chdir(dir)
            else:
                continue

            nnodes = float(dir.replace(DIR_PREFIX, ''))

            fullstr = "%19s |" % dir

            # test if there is anything to analyze
            if len(glob.glob("screen.out.*")) == 0:
                if (has_epetra and os.path.exists(analysis_run + ".epetra")) or \
                   (has_tpetra and os.path.exists(analysis_run + ".tpetra")) or \
                   (has_ml and os.path.exists(analysis_run + ".ml")):
                    fullstr += " running now?"
                else:
                    fullstr += " not run"
                print(fullstr)
                os.chdir("..")
                continue

            for s in timelines:
                if has_epetra:
                    try:
                        rtime_epetra = []
                        for epetra_file in sort_nicely(glob.glob(analysis_run + "*.epetra")):
                            if (parsefunc == ""):
                                r = commands.getstatusoutput("grep \"" + s + "\" " + epetra_file)
                                if r[0] != 0:
                                    raise RuntimeError("Error reading \"" + analysis_run + ".epetra")
                                grep_res = r[1]

                                r = commands.getstatusoutput("echo \"" + grep_res + "\" | tail -n 1 | awk '{print $(NF-4)}'")
                                try:
                                    rtime_epetra.append(float(r[1]))
                                except (ValueError):
                                    # check for serial version (it outputs a single timer, compared to multiple in parallel (min, max, ...))
                                    r = commands.getstatusoutput("echo \"" + grep_res + "\" | tail -n 1 | awk '{print $(NF-1)}'")
                                    rtime_epetra.append(float(r[1]))

                            else:
                                r = commands.getstatusoutput(parsefunc(epetra_file, s))
                                rtime_epetra.append(float(r[1]))


                        time_epetra[s] = stat_time(rtime_epetra)

                        if nnodes == BASECASE:
                            basetime_epetra[s] = time_epetra[s]
                        if scaling_type == "weak":
                            eff_epetra[s] = 100 * basetime_epetra[s] / time_epetra[s]
                        else:
                            eff_epetra[s] = 100 * (basetime_epetra[s] * BASECASE) / (time_epetra[s] * nnodes)
                        fullstr += "%13.2f %7.2f%%" % (time_epetra[s], eff_epetra[s])

                    except (RuntimeError, ValueError):
                        # print("Problem converting \"%s\" to float for timeline \"%s\" in \"%s\"" % (r[1], s, epetra_file))
                        fullstr += "           -   -"

                if has_tpetra:
                    try:
                        rtime_tpetra = []
                        for tpetra_file in sort_nicely(glob.glob(analysis_run + "*.tpetra")):
                            if (parsefunc == ""):
                                r = commands.getstatusoutput("grep \"" + s + "\" " + tpetra_file)
                                if r[0] != 0:
                                    raise RuntimeError("Error reading \"" + analysis_run + ".tpetra")
                                grep_res = r[1]

                                r = commands.getstatusoutput("echo \"" + grep_res + "\" | tail -n 1 | awk '{print $(NF-4)}'")
                                try:
                                    rtime_tpetra.append(float(r[1]))
                                except (ValueError):
                                    # check for serial version (it outputs a single timer, compared to multiple in parallel (min, max, ...))
                                    r = commands.getstatusoutput("echo \"" + grep_res + "\" | tail -n 1 | awk '{print $(NF-1)}'")
                                    rtime_tpetra.append(float(r[1]))

                            else:
                                r = commands.getstatusoutput(parsefunc(tpetra_file, s))
                                rtime_tpetra.append(float(r[1]))

                        time_tpetra[s] = stat_time(rtime_tpetra)

                        if nnodes == BASECASE:
                            basetime_tpetra[s] = time_tpetra[s]

                        if scaling_type == "weak":
                            eff_tpetra[s] = 100 * basetime_tpetra[s] / time_tpetra[s]
                        else:
                            eff_tpetra[s] = 100 * (basetime_tpetra[s] * BASECASE) / (time_tpetra[s] * nnodes)
                        fullstr += "%13.2f %7.2f%%" % (time_tpetra[s], eff_tpetra[s])

                    except (RuntimeError, ValueError):
                        #print("Problem converting \"%s\" to float for timeline \"%s\" in \"%s\"" % (r[1], s, tpetra_file))
                        fullstr += "           -   -"

                if has_ml:
                    ml_file = analysis_run + ".ml"
                    if (parsefunc == ""):
                      return "Error: no parsing function provided"
                    else:
                      theCommand = parsefunc(ml_file,s)

                    r = commands.getstatusoutput(theCommand)

                    # Handle multiple timers w/ same name.  This splits last entry in tuple by line breaks into
                    # an array of strings.  The string array is then converted ("mapped") into an array of floats.
                    tt = map(float, r[-1].split())
                    time_ml[s] = sum(tt)

                    if nnodes == BASECASE:
                        basetime_ml[s] = time_ml[s]
                    eff_ml[s] = 100 * basetime_ml[s] / time_ml[s]
                    fullstr += "%13.2f %7.2f%%" % (time_ml[s], eff_ml[s])

            print(fullstr)

            os.chdir("..")

def build(nnodes, nx, binary, petra, matrix, datafiles, cmds, template, output, cpn, cps, unified, num_runs):
    dir = DIR_PREFIX + str(nnodes)
    print("Building %s..." % dir)

    ensure_dir(dir)                                         # make directory if necessary
    os.symlink("../" + binary, dir + "/binary.exe")         # link binary file
    for afile in datafiles:                                 # copy all xml files
        shutil.copy(afile, dir)

    sched_args = ""
    if PLATFORM == "hopper":
        sched_args  = "aprun"
        sched_args += " -ss"                                    # Demands strict memory containment per NUMA node
        sched_args += " -cc numa_node"                          # Controls how tasks are bound to cores and NUMA nodes
        sched_args += " -N " + str(cpn)                         # Number of tasks per node
        if cps != -1:
            sched_args += " -S " + str(cps)                     # Number of tasks per NUMA node (note: hopper has 4 NUMA nodes)

    elif PLATFORM == "shannon":
        # There some issues on Shannon with openmpi and srun, so we use mpirun instead
        sched_args  = "mpirun"
        sched_args += " --npernode " + str(cpn)
        sched_args += " --bind-to-core"
        # sched_args += " --map-by numa"                        # This conflicts with --npernode

    elif PLATFORM == "volta":
        sched_args  = "srun"
        sched_args += " --ntasks-per-node " + str(cpn)
        if cps != -1:
            sched_args += " --ntasks-per-socket " + str(cps)

    script_path = dir + "/" + PETRA_PREFIX + str(nnodes) + "." + SCHEDULER

    full_template = "sched.full.template"
    os.system("echo -e \"" + SCHED_HEADER + "\" > " + full_template)
    os.system("cat " + template + " >> " + full_template)

    # construct batch script from template
    os_cmd = "cat " + full_template
    os_cmd += (" | sed \"s/_SCHED_ARGS_/" + sched_args                + "/\"" +
               " | sed \"s/_WIDTH_/"      + str(nnodes*CPN)           + "/\"" +
               " | sed \"s/_NODES_/"      + str(nnodes)               + "/\"" +
               " | sed \"s/_CORES_/"      + str(nnodes*cpn)           + "/\"" +
               " | sed \"s/_MTYPE_/"      + matrix                    + "/\"" +
               " | sed \"s/_NX_/"         + str(nx)                   + "/g\"" +
               " | sed \"s/_EPETRA_/"     + str(petra & 1)            + "/\"" +
               " | sed \"s/_TPETRA_/"     + str(petra & 2)            + "/\"" +
               " | sed \"s/_UNIFIED_/"    + str(unified)              + "/\"" +
               " | sed \"s/_NUM_RUNS_/"   + str(num_runs)             + "/\"" +
               " | sed \"s/_NUM_CMDS_/"   + str(len(cmds))            + "/\"")

    for i in range(len(cmds)):
        os_cmd += " | sed \"s/_CMD" + str(i+1) + "_/" + cmds[i] + "/g\""
    os_cmd += " >> " + script_path
    os.system(os_cmd)

    os.system("rm " + full_template)

def clean():
    for dir in sort_nicely(glob.glob(DIR_PREFIX + "*")):
        print("Cleaning %s..." % dir)
        shutil.rmtree(dir)

def run():
    for dir in sort_nicely(glob.glob(DIR_PREFIX + "*")):
        print("Running %s..." % dir)

        os.chdir(dir)
        if   SCHEDULER == "pbs":
            os.system("qsub "   + PETRA_PREFIX + dir.replace(DIR_PREFIX, '') + "." + SCHEDULER)
        elif SCHEDULER == "slurm":
            os.system("sbatch " + PETRA_PREFIX + dir.replace(DIR_PREFIX, '') + "." + SCHEDULER)
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

def import_from(module, name):
    module = __import__(module, fromlist=[name])
    return getattr(module, name)

def stat_time(times):
    return min(times)

# ========================= main =========================
def main():
    global PLATFORM
    global CPN
    global SCHEDULER
    global SCHED_HEADER

    # Try to detect platform
    r = commands.getstatusoutput("hostname")
    if r[0] != 0:
        raise RuntimeError("Cannot run \"hostname\"")
    hostname = r[1]
    if   hostname[0:6] == "hopper":
        PLATFORM = "hopper"
    elif hostname[0:6] == "edison":
        PLATFORM = "edison"
    elif hostname[0:7] == "shannon":
        PLATFORM = "shannon"
    elif hostname[0:5] == "volta":
        PLATFORM = "volta"
    else:
        raise RuntimerError("Unable to detect platform")

    # We double escape the \n, as it is later
    # passed to os.system
    if PLATFORM == "hopper":
        CPN          = 24
        SCHEDULER    = "pbs"
        SCHED_HEADER = ("#PBS -S /bin/bash\\n"
                        "#PBS -q regular\\n" +                     # queue
                        "#PBS -l mppwidth=_WIDTH_\\n"
                        "#PBS -l walltime=00:10:00\\n"
                        "#PBS -N _MTYPE___CORES_\\n"
                        "#PBS -e screen.err.\\$PBS_JOBID\\n"
                        "#PBS -o screen.out.\\$PBS_JOBID\\n"
                        "#PBS -V\\n"
                        "#PBS -A m1327\\n"                         # repository
                        "#PBS -m ae\\n\\n"                         # abort/termination
                        "cd \\$PBS_O_WORKDIR\\n")

    elif PLATFORM == "shannon" or PLATFORM == "volta":
        if   PLATFORM == "shannon":
            CPN      = 16
        elif PLATFORM == "volta":
            CPN      = 24
        SCHEDULER    = "slurm"
        SCHED_HEADER = ("#!/bin/bash\\n"
                        "#SBATCH -N _NODES_\\n"
                        "#SBATCH --exclusive\\n"
                        "#SBATCH -t 00:10:00\\n"
                        "#SBATCH -J _MTYPE___CORES_\\n"
                        "#SBATCH -e screen.err.%J\\n"
                        "#SBATCH -o screen.out.%J\\n")

    else:
        print("Unknown platform type \"%s\"" % PLATFORM)
        return

    controller()

if __name__ == '__main__':
    main()
