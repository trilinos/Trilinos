#!/usr/bin/python
import os
import sys
import math
import subprocess

def getstatusoutput(cmd):
    """Return (status, output) of executing cmd in a shell."""
    pipe = os.popen(cmd + ' 2>&1', 'r')
    text = pipe.read()
    sts  = pipe.close()
    if sts is None: sts = 0
    if text[-1:] == '\n': text = text[:-1]
    return sts, text

def pipepager(text, cmd):
    """Page through text by feeding it to another program."""
    pipe = os.popen(cmd, 'w')
    try:
        pipe.write(text)
        pipe.close()
    except IOError:
        pass # Ignore broken pipes caused by quitting the pager program.

def deleteDir(path):
    r = getstatusoutput("rm -rf " + path)
    if r[0] != 0:
        raise RuntimeError(r[1])

def createDir(path):
    r = getstatusoutput("mkdir " + path)
    if r[0] != 0:
        raise RuntimeError(r[1])

def runCommand(cmd):
    r = getstatusoutput(cmd)
    # if r[0] != 0:
    #     raise RuntimeError(r[1])
    return r[1]

def clearWindow():
    os.system('cls' if os.name == 'nt' else 'clear')

def waitForKey():
    os.system('read -s -n 1 -p "Press any key to continue..."')
    print

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# some colors
class bc:
    HEADER      = '\033[95m'
    OKBLUE      = '\033[94m'
    OKGREEN     = '\033[92m'
    OKDARKGREEN = '\033[32m'
    WARNING     = '\033[93m'
    FAIL        = '\033[91m'
    ENDC        = '\033[0m'

def disable(self):
    self.HEADER      = ''
    self.OKBLUE      = ''
    self.OKGREEN     = ''
    self.OKDARKGREEN = ''
    self.WARNING     = ''
    self.FAIL        = ''
    self.ENDC        = ''


class ProblemHandler():
  """Class for handling demonstration problems"""

  def __init__(self):
    self.problem     = "Laplace 2D"
    self.executable  = "MueLu_laplace2d.exe"
    self.meshx       = 50
    self.meshy       = 50
    self.numprocs    = 2
    self.xmlFileName = "xml/s2_easy.xml"

    self.proc = []
    self.proc.append(subprocess.Popen(['gnuplot', '--persist'], shell=True, stdin=subprocess.PIPE, ))
    self.proc.append(subprocess.Popen(['gnuplot', '--persist'], shell=True, stdin=subprocess.PIPE, ))
    self.proc.append(subprocess.Popen(['gnuplot', '--persist'], shell=True, stdin=subprocess.PIPE, ))
    self.proc.append(subprocess.Popen(['gnuplot', '--persist'], shell=True, stdin=subprocess.PIPE, ))

    self.isDirty = True # flag to store whether or not problem must be run again

    self.editor = "kwrite"    # TODO replace me by local editor...

  def main(self):
    self.printMainMenu()

  def runMenu(self,options,callbacks):
    # display all options
    for i,option in enumerate(options):
      print('%s. %s' % (i, option))

    # select an option
    choice = raw_input('your choice? ')
    if is_number(str(choice)) and int(choice) < len(options):
      callbacks[int(choice)]() # call corresponding function
    else:
      print("ups: choice = " + str(choice) + " len(option)=" + str(len(option)))



  def doLaplace2Dn(self):
    self.problem    = "Laplace 2D"
    self.executable = "MueLu_laplace2d.exe"
    self.meshx      = raw_input("Mesh: Elements in x direction = ")
    self.meshy      = raw_input("Mesh: Elements in y direction = ")
    self.runLaplaceProblem()

  def doLaplace2D50(self):
    self.problem    = "Laplace 2D"
    self.executable = "MueLu_laplace2d.exe"
    self.meshx      = 50
    self.meshy      = 50
    self.runLaplaceProblem()

  def doRecirc2Dn(self):
    self.problem    = "Recirc 2D"
    self.executable = "MueLu_recirc2d.exe"
    self.meshx      = raw_input("Mesh: Elements in x direction = ")
    self.meshy      = raw_input("Mesh: Elements in y direction = ")
    self.runLaplaceProblem() # we can use the same routine as for Laplace...

  def doRecirc2D50(self):
    self.problem    = "Recirc 2D"
    self.executable = "MueLu_recirc2d.exe"
    self.meshx      = 50
    self.meshy      = 50
    self.runLaplaceProblem() # we can use the same routine as for Laplace...

  def runLaplaceProblem(self):
    # check whether xml file exists
    while self.xmlFileName == "" or not os.path.isfile(self.xmlFileName) or not os.access(self.xmlFileName, os.R_OK):
      print(bc.FAIL + "Solver xml parameters: " + bc.ENDC + str(self.xmlFileName) + bc.FAIL + " invalid" + bc.ENDC)
      m = MueLu_XMLgenerator()
      m.askForSolver()
      self.xmlFileName = m.xmlFileName # store xml file

    while True:
      self.printActionMenu()

  def printActionMenu(self):
    options   = ['Rerun simulation',   'Show screen output',    'Change input',  'Open xml file', 'Change number of processors',   'Plot solution',   'Postprocess aggregates',             'Exit']
    callbacks = [   self.runExample, self.printScreenOutput, self.changeSolver, self.openXMLfile,             self.changeProcs, self.plotSolution, self.postprocessAggregates, self.doExitProgram]

    while True:
      clearWindow()
      self.printSettings()
      print("")
      if self.isDirty == True:
        print(bc.FAIL + "DO NOT FORGET TO RUN THE EXAMPLE (option 0)" + bc.ENDC)
      else:
        print(bc.OKDARKGREEN + "Results up to date!" + bc.ENDC)
      print("")
      self.runMenu(options, callbacks)

  def runExample(self):
    # runs example
    print("PREPARE SIMULATON")
    runCommand("rm *.vtp *.mat example*.txt output.log aggs*.out nodes*.txt")

    print("RUN EXAMPLE")
    cmd = "mpirun -np %d %s --nx=%d --ny=%d --xml=%s | tee output.log 2>&1" % (self.numprocs, self.executable, self.meshx, self.meshy, self.xmlFileName)
    print(cmd)
    runCommand(cmd)
    runCommand("echo 'Press q to return.' >> output.log")

    print("POSTPROCESSING...")
    runCommand("cat example*.txt > example.txt")

    print("COMPLETE")
    self.isDirty = False
    waitForKey()

  def plotSolution(self):
    for i in range(0,4):
      self.proc[i].stdin.write("set term x11 " + str(i+1))
      self.proc[i].stdin.write("set dgrid3d " + str(self.meshy) + "," + str(self.meshx) + "\n")
      self.proc[i].stdin.write("set style data lines\n")
      self.proc[i].stdin.write("set nolabel\n")
      self.proc[i].stdin.write("set key off\n")
      self.proc[i].stdin.write("set autoscale\n")

      if   i == 0:
        self.proc[0].stdin.write("set title \"Solution\"\n")
        self.proc[0].stdin.write("splot \"example.txt\" using 3:4:5\n")

      elif i == 1:
        self.proc[1].stdin.write("set title \"Multigrid solution\"\n")
        self.proc[1].stdin.write("splot \"example.txt\" using 3:4:7\n")

      elif i == 2:
        self.proc[2].stdin.write("set title \"Error (Exact vs. Multigrid)\"\n")
        self.proc[2].stdin.write("set palette model RGB defined ( 0 'black', 1 'white')\n")
        self.proc[2].stdin.write("set hidden3d\n")
        self.proc[2].stdin.write("set style line 1 lt 4 lw .5\n")
        self.proc[2].stdin.write("set pm3d\n")
        self.proc[2].stdin.write("splot \"example.txt\" using 3:4:($ 5-$ 7) with lines palette\n")

      elif i == 3:
        self.proc[3].stdin.write("set title \"Distribution of processors\"\n")
        self.proc[3].stdin.write("set palette model RGB defined ( 0 'red', 1 'green', 2 'blue', 3 'yellow', 4 'pink')\n")
        self.proc[3].stdin.write("set hidden3d\n")
        self.proc[3].stdin.write("set style line 1 lt 4 lw .5\n")
        self.proc[3].stdin.write("set pm3d\n")
        self.proc[3].stdin.write("splot \"example.txt\" using 3:4:1 with points palette\n")

      # self.proc[i].stdin.write("quit\n") #close the gnuplot window
      self.proc[i].stdin.flush()

  def postprocessAggregates(self):
    # check whether "example.txt" is available
    if os.path.isfile("example.txt") == False:
      print(bc.FAIL + "Simulation data not available. Run the simulation first." + bc.ENDC)
      waitForKey()
      return

    if os.path.isfile("aggs_level0_proc0.out") == False:
      print(bc.FAIL + "No aggregation debug output found. Do not forget to turn on the \"aggregation: visualize\" in your xml file." + bc.ENDC)
      waitForKey()
      return

    if os.path.isfile("MueLu_Agg2VTK.py"):
      os.remove("MueLu_Agg2VTK.py")
    o = open("MueLu_Agg2VTK.py","a")
    for line in open("tmpl/MueLu_Agg2VTK.py_TMPL"):
      line = line.replace("$NUMPROCS", str(self.numprocs))
      o.write(line)
    o.close()

    print("POSTPROCESS AGGREGATION OUTPUT DATA")
    runCommand("chmod 750 ./MueLu_Agg2VTK.py")
    print(runCommand("./MueLu_Agg2VTK.py"))

    if os.path.isfile("aggs0.vtp") == False:
      print(bc.WARNING + "Seems that the postprocessing failed (vtp files could not be created)." + bc.ENDC)
      waitForKey()
      return

    print(bc.OKDARKGREEN + "Use paraview to visualize generated vtk files for aggregates." + bc.ENDC)
    waitForKey()

  def printScreenOutput(self):
    clearWindow()
    if not os.path.isfile("output.log") or not os.access("output.log", os.R_OK):
      print(bc.FAIL + "Screen output not available." + bc.ENDC)
    else:
      pipepager(runCommand("cat output.log"), "less")

  def openXMLfile(self):
    editor = subprocess.Popen([self.editor + " " + self.xmlFileName], shell=True, stdin=subprocess.PIPE, )

  def printProblemSelectionMenu(self):
    options   = ['Laplace 2D (50x50)',      'Laplace 2D', 'Recirc 2D (50x50)',      'Recirc 2D',             'Exit']
    callbacks = [self.doLaplace2D50,   self.doLaplace2Dn,   self.doRecirc2D50, self.doRecirc2Dn, self.doExitProgram]
    while True:
      self.runMenu(options,callbacks)

  def changeSolver(self):
    self.xmlFileName = raw_input("XML file name: ")
    self.isDirty = True
    while self.xmlFileName == "" or not os.path.isfile(self.xmlFileName) or not os.access(self.xmlFileName, os.R_OK):
      print(bc.FAIL + "Solver xml parameters: " + bc.ENDC + str(self.xmlFileName) + bc.FAIL + " invalid" + bc.ENDC)
      m = MueLu_XMLgenerator()
      m.xmlFileName = (self.xmlFileName if self.xmlFileName != "" else "xml/custom.xml")
      m.generateXMLfile()
      m.askForSolver()
      m.generateXMLfile()
      self.xmlFileName = m.xmlFileName # store xml file

  def changeProcs(self):
    self.numprocs = raw_input("Number of processors: ")
    while not is_number(str(self.numprocs)):
      self.numprocs = raw_input("Number of processors: ")

    self.numprocs = int(self.numprocs)
    self.isDirty  = True

  def printMainMenu(self):
    clearWindow()
    self.printSettings()
    print("")
    print("")
    while True:
      self.printProblemSelectionMenu()

  def doExitProgram(self):
    print("CLEAN UP temporary data")
    runCommand("rm *.vtp *.mat example*.txt output.log aggs*.out nodes*.txt")
    print("QUIT")
    sys.exit()

  def printSettings(self):
    ## print out all made settings for xml file
    print(bc.HEADER  + "***************************   PROBLEM   ****************************" + bc.ENDC)
    print(bc.WARNING + "Problem type                       : " + bc.ENDC + str(self.problem))
    print(bc.WARNING + "Mesh size                          : " + bc.ENDC + str(self.meshx) + " x " + str(self.meshy))
    print("")
    if self.xmlFileName == "" or not os.path.isfile(self.xmlFileName) or not os.access(self.xmlFileName, os.R_OK):
      print(bc.FAIL    + "Solver xml parameters              : " + bc.ENDC + str(self.xmlFileName) + bc.FAIL + " invalid" + bc.ENDC)
    else:
      print(bc.WARNING + "Solver xml parameters              : " + bc.ENDC + str(self.xmlFileName))
    print(bc.WARNING + "Number of processors               : " + bc.ENDC + str(self.numprocs))
    print(bc.HEADER  + "***************************   PROBLEM   ****************************" + bc.ENDC)

class MueLu_XMLgenerator():
  """Simple generator for MueLu XML files."""

  def __init__(self):
    # common MG settings
    self.xmlFileName            = ""        # output file name for xml data
    self.maxMultLevels          = 5         # maximum number of levels
    self.maxCoarseSize          = 1000      # max coarse size

    # aggregate settings
    self.dropTolerance          = 0.0
    self.minAggSize             = 4
    self.maxAggSize             = 9

    # smoother settings
    self.levelSmoother          = "Jacobi"
    self.levelSmootherSweeps    = 1
    self.levelSmootherDamp      = 0.7
    self.coarseSolver           = "Direct"

    # transfer operators
    self.mgAlgo                 = "sa"
    self.saDamp                 = 1.33

    # restriction operators
    self.symProblem             = True

    # rebalancing
    self.doRebalancing          = False
    self.minRowsPerProc         = 800
    self.nnzImbalance           = 1.1
    self.rebStartLevel          = 1

    self.isDirty                = True      # flag to store, whether changes have been saved or not
    self.exitLoop               = False     # set to true to exit current loop

    print(bc.FAIL + '====================================================================================' + bc.ENDC)
    print(               '====================================================================================' + bc.ENDC)

  def main(self):
    #self.view, self.exit_view = self.setup_view()
    #self.loop = urwid.MainLoop(self.view, self.palette,
    #    unhandled_input=self.unhandled_input)
    #self.loop.run()
    while True:
      self.printMainMenu()

  def askForSolver(self):
    self.exitLoop = False
    while self.exitLoop == False:
      self.printMainMenu()

  def doFileName(self):
    self.xmlFileName = raw_input("XML file name: ")
    self.isDirty = True

  def doMaxLevels(self):
    self.maxMultLevels = raw_input("Max multigrid levels: ")
    self.isDirty = True

  def doMaxCoarseSize(self):
    self.maxCoarseSize = raw_input("Max coarse size: ")
    self.isDirty = True

  def doRelaxationSweeps(self):
    self.levelSmootherSweeps = raw_input("Smoother sweeps: ")
    self.isDirty  = True

  def doRelaxationJacobi(self):
    self.levelSmoother = "Jacobi"
    self.levelSmootherSweeps = raw_input("Smoother sweeps: ")
    self.levelSmootherDamp   = raw_input("Smoother damping: ")
    self.isDirty = True

  def doRelaxationGS(self):
    self.levelSmoother = "Gauss-Seidel"
    self.levelSmootherSweeps = raw_input("Smoother sweeps: ")
    self.levelSmootherDamp   = raw_input("Smoother damping: ")
    self.isDirty = True

  def doRelaxationSymGS(self):
    self.levelSmoother = "Sym.Gauss-Seidel"
    self.levelSmootherSweeps = raw_input("Smoother sweeps: ")
    self.levelSmootherDamp   = raw_input("Smoother damping: ")
    self.isDirty = True

  def doDropTolerance(self):
    self.dropTolerance = raw_input("Drop tolerance for matrix graph (default = 0.0): ")
    self.isDirty       = True
  def doMinAggSize(self):
    self.minAggSize    = raw_input("Minimum number of nodes per aggregate: ")
    self.isDirty       = True
  def doMaxAggSize(self):
    self.maxAggSize = raw_input("Maximum number of nodes per aggregate: ")
    self.isDirty       = True

  # Transfer operators
  def doPaAMG(self):
    self.mgAlgo  = "unsmoothed"
    self.isDirty = True
  def doSaAMG(self):
    self.mgAlgo  = "sa"
    self.saDamp  = raw_input("Transfer operator damping: ")
    self.isDirty = True
  def doPgAMG(self):
    self.mgAlgo  = "pg"
    self.isDirty = True
  def doEmin(self):
    self.mgAlgo  = "emin"
    self.isDirty = True

  # Restriction operators
  def doSymR(self):
    self.symProblem = True
    self.isDirty    = True
  def doNonsymR(self):
    self.symProblem = False
    self.isDirty    = True

  # Rebalancing
  def doRebalancingOption(self):
    self.doRebalancing  = True
    self.minRowsPerProc = raw_input("Minimum number of DOFs per processor   : ")
    self.nnzImbalance   = raw_input("Max nonzero imbalance [default 1.1]    : ")
    self.rebStartLevel  = raw_input("Start rebalancing on level [default 1] : ")
    self.isDirty        = True

  def doNoRebalancingOption(self):
    self.doRebalancing = False
    self.isDirty       = True

  def runMenu(self, options, callbacks):
    for i,option in enumerate(options):
      print('%s. %s' % (i, option)) # display all options
    choice = raw_input('your choice? ')
    if is_number(str(choice)) and int(choice) < len(options):
      callbacks[int(choice)]() # call correspondending function

  def doCommonMenu(self):
    options   = ['Max multigrid levels',    'Max coarse size',            'Back']
    callbacks = [      self.doMaxLevels, self.doMaxCoarseSize, self.askForSolver]
    while self.exitLoop == False:
      self.runMenu(options, callbacks)
    self.exitLoop = True #False

  def doAggregatesMenu(self):
    options   = [     'Drop tolerance', 'Min aggregate size', 'Max aggregate size',            'Back']
    callbacks = [self.doDropTolerance,     self.doMinAggSize,    self.doMaxAggSize, self.askForSolver]
    while self.exitLoop == False:
      self.runMenu(options, callbacks)
    self.exitLoop = True #False

  def doSmootherMenu(self):
    options   = [      'Change MG sweeps',                'Jacobi',      'Gauss-Seidel',    'Sym. Gauss-Seidel',            'Back']
    callbacks = [ self.doRelaxationSweeps, self.doRelaxationJacobi, self.doRelaxationGS, self.doRelaxationSymGS, self.askForSolver]
    self.runMenu(options, callbacks)

  def doTransferMenu(self):
    options   = ['Smoothed transfer (SA-AMG)', 'Non-smoothed transfer (PA-AMG)', 'Smoothed transfer (PG-AMG)', 'Energy minimization',            'Back']
    callbacks = [                self.doSaAMG,                     self.doPaAMG,                 self.doPgAMG,           self.doEmin, self.askForSolver]
    self.runMenu(options, callbacks)

  def doRestrictorMenu(self):
    options   = ['Symmetric', 'Non-symmetric',            'Back']
    callbacks = [self.doSymR,  self.doNonsymR, self.askForSolver]
    self.runMenu(options, callbacks)

  def doRebalancingMenu(self):
    options   = [          'No rebalancing',   'Activate rebalancing',            'Back']
    callbacks = [self.doNoRebalancingOption, self.doRebalancingOption, self.askForSolver]
    self.runMenu(options, callbacks)

  def doExitProgram(self):
    #sys.exit()
    print("EXIT")
    self.exitLoop = True

  def printMainMenu(self):
    clearWindow()
    self.printSettings()
    print("")
    print("")

    options   = ['General multigrid settings', 'Aggregation settings', 'Smoothing settings', 'Transfer settings', 'Restriction operators', 'Rebalancing settings', 'Save settings (XML)',              'Back']
    callbacks = [           self.doCommonMenu,  self.doAggregatesMenu,  self.doSmootherMenu, self.doTransferMenu,   self.doRestrictorMenu, self.doRebalancingMenu,  self.generateXMLfile, self.doExitProgram]

    self.runMenu(options, callbacks)

  def printSettings(self):
    ## print out all made settings for xml file
    print(bc.HEADER  + "***************************   SETTINGS   ****************************" + bc.ENDC)
    print(bc.WARNING + "XML file name:           : " + bc.ENDC + str(self.xmlFileName))
    print("")
    print(bc.WARNING + "Max multigrid levels     : " + bc.ENDC + str(self.maxMultLevels))
    print(bc.WARNING + "Max coarse size          : " + bc.ENDC + str(self.maxCoarseSize))
    print("")
    print(bc.WARNING + "Level smoother           : " + bc.ENDC + str(self.levelSmoother))
    print(bc.WARNING + "Level smoothing sweeps   : " + bc.ENDC + str(self.levelSmootherSweeps))
    print(bc.WARNING + "Level damping parameter  : " + bc.ENDC + str(self.levelSmootherDamp))
    print("")
    print(bc.WARNING + "Coarse solver            : " + bc.ENDC + str(self.coarseSolver))
    print("")
    print(bc.WARNING + "Graph drop tolerance     : " + bc.ENDC + str(self.dropTolerance))
    print(bc.WARNING + "Aggregate size (min/max) : " + bc.ENDC + str(self.minAggSize) + "/" + str(self.maxAggSize))
    print("")
    print(bc.WARNING + "Transfer operators       : " + bc.ENDC + str(self.mgAlgo))
    print(bc.WARNING + "Transfer smoothing par.  : " + bc.ENDC + str(self.saDamp))
    print("")
    print(bc.WARNING + "Symmetric problem        : " + bc.ENDC + ("true" if self.symProblem else "false"))
    if self.doRebalancing == False:
      print(bc.WARNING + "NO Rebalancing" + bc.ENDC)
    else:
      print(bc.WARNING + "Rebalancing active       : " + bc.ENDC)
      print(bc.WARNING + "Minimum DOFs per proc    : " + bc.ENDC + str(self.minRowsPerProc))
      print(bc.WARNING + "Nonzero imbalance        : " + bc.ENDC + str(self.nnzImbalance))
      print(bc.WARNING + "Start level for rebal.   : " + bc.ENDC + str(self.rebStartLevel))
    print(bc.HEADER  + "***************************   SETTINGS   ****************************" + bc.ENDC)
    print("")

    if self.isDirty == True:
      print(bc.FAIL        + "CHANGES HAVE NOT BEEN SAVED!" + bc.ENDC)
    else:
      print(bc.OKDARKGREEN + "CHANGES HAVE BEEN SAVED!"     + bc.ENDC)

  def generateXMLfile(self):
    # generate HEAD file for pre_exodus
    if os.path.isfile(self.xmlFileName):
      os.remove(self.xmlFileName)

    o = open(self.xmlFileName, "a")

    for line in open("tmpl/muelu_easy.xml_TMPL"):
      if self.levelSmoother == 'Jacobi'         or \
        self.levelSmoother == 'Gauss-Seidel'    or \
        self.levelSmoother == 'Sym.Gauss-Seidel':
        line = line.replace("$IF_CHEBYSHEV"        , "<!--")
        line = line.replace("$ENDIF_CHEBYSHEV"     , "-->")
        line = line.replace("$IF_RELAXATION"       ,   "");
        line = line.replace("$ENDIF_RELAXATION"    , "");
        line = line.replace("$SMOOTHER"            , "RELAXATION")
        line = line.replace("$RELAXATION"          , str(self.levelSmoother))
        line = line.replace("$SMOO_SWEEPS"          , str(self.levelSmootherSweeps))
        line = line.replace("$SMOO_DAMP"            , str(self.levelSmootherDamp))
      else:
        line = line.replace("$IF_CHEBYSHEV"        , "<!--")
        line = line.replace("$ENDIF_CHEBYSHEV"     , "-->")
        line = line.replace("$IF_RELAXATION"       , "");
        line = line.replace("$ENDIF_RELAXATION"    , "");
        line = line.replace("$SMOOTHER"            , "CHEBYSHEV")

      line = line.replace("$MAXLEVELS"              , str(self.maxMultLevels))
      line = line.replace("$MAXCOARSESIZE"          , str(self.maxCoarseSize))

      line = line.replace("$SYMPROBLEM"             , "true" if self.symProblem else "false")
      line = line.replace("$AMG"                    , str(self.mgAlgo))
      line = line.replace("$SADAMPING"              , str(self.saDamp))

      line = line.replace("$DROPTOL"                , str(self.dropTolerance))
      line = line.replace("$MINAGGS"                , str(self.minAggSize))
      line = line.replace("$MAXAGGS"                , str(self.maxAggSize))

      line = line.replace("$REPARTITION"            , "true" if self.doRebalancing else "false")

      o.write(line)
    o.close()
    self.isDirty = False

if __name__ == '__main__':
  #MueLu_XMLgenerator().main()
  ProblemHandler().main()
