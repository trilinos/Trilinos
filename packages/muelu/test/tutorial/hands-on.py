#!/usr/bin/env python
import os
import sys
import math
import subprocess

# Fix for input vs raw_input for Python 2.x compatibility
try: input = raw_input
except NameError: pass

def getstatusoutput(cmd):
    """Return (status, output) of executing cmd in a shell."""
    pipe = os.popen(cmd + ' 2>&1', 'r')
    text = pipe.read()
    sts = pipe.close()
    if sts is None: sts = 0
    if text[-1:] == '\n': text = text[:-1]
    return sts, text


def deleteDir(path):
    """deletes the path entirely"""
    cmd = "rm -rf "+path
    result = getstatusoutput(cmd)
    if(result[0]!=0):
        raise RuntimeError(result[1])

def createDir(path):
    """deletes the path entirely"""
    cmd = "mkdir "+path
    result = getstatusoutput(cmd)
    if(result[0]!=0):
        raise RuntimeError(result[1])

def runCommand(cmd):
    """deletes the path entirely"""
    result = getstatusoutput(cmd)
    #if(result[0]!=0):
    #    raise RuntimeError(result[1])
    return result[1]

def clearWindow():
    #print("CLEAR")
    os.system('cls' if os.name == 'nt' else 'clear')

def waitForKey():
    os.system("""bash -c 'read -s -n 1 -p "Press any key to continue..."'""")
    print()

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# text output colors
class usecolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    OKDARKGREEN = '\033[32m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

# no output colors (primarily for logging output to file)
class nocolors:
    HEADER = ''
    OKBLUE = ''
    OKGREEN = ''
    OKDARKGREEN = ''
    WARNING = ''
    FAIL = ''
    ENDC = ''

def disable(self):
    self.HEADER = ''
    self.OKBLUE = ''
    self.OKGREEN = ''
    self.OKDARKGREEN = ''
    self.WARNING = ''
    self.FAIL = ''
    self.ENDC = ''


class ProblemHandler():
  """Class for handling demonstration problems"""

  def __init__(self):
    self.problem    = "Laplace2D"
    self.solver     = "cg"
    self.executable = "MueLu_TutorialDriver.exe"
    self.bcolors    = usecolors()
    self.meshx      = 50
    self.meshy      = 50
    self.mgsweeps   = 1
    self.numprocs   = 2
    self.xmlFileName = "s2a.xml"

    self.proc1 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, text=True)
    self.proc2 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, text=True)
    self.proc3 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, text=True)
    self.proc4 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, text=True)
    self.proc5 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, text=True)

    self.isDirty = True                   # flag to store, whether problem has to be rerun or not

    self.editor = os.environ.get("H_EDITOR", "nano") # text editor is specified by the env var H_EDITOR, nano by default

  def main(self):
    self.printMainMenu()

  def runMenu(self,options,callbacks):
    for i,option in enumerate(options):
      print('%s. %s' % (i, option)) # display all options
    choice = input('your choice? ')
    if is_number(str(choice)):
      if int(choice) < len(options) and int(choice) > -1:
        callbacks[int(choice)]() # call corresponding function
      else:
        print("Option is out of range!")
    else:
      print("Invalid option! " + str(choice))



  def doLaplace2Dn(self):
    self.problem    = "Laplace2D"
    self.executable = "MueLu_TutorialDriver.exe"
    self.solver     = "cg"
    self.meshx      = input("Mesh: Elements in x direction = ")
    self.meshy      = input("Mesh: Elements in y direction = ")
    self.runLaplaceProblem()

  def doLaplace2D50(self):
    self.problem    = "Laplace2D"
    self.executable = "MueLu_TutorialDriver.exe"
    self.solver     = "cg"
    self.meshx      = 50
    self.meshy      = 50
    self.runLaplaceProblem()

  def doRecirc2Dn(self):
    self.problem    = "Recirc2D"
    self.executable = "MueLu_TutorialDriver.exe"
    self.solver     = "gmres"
    self.meshx      = input("Mesh: Elements in x direction = ")
    self.meshy      = input("Mesh: Elements in y direction = ")
    self.runLaplaceProblem() # we can use the same routine as for Laplace...

  def doRecirc2D50(self):
    self.problem    = "Recirc2D"
    self.executable = "MueLu_TutorialDriver.exe"
    self.solver     = "gmres"
    self.meshx      = 50
    self.meshy      = 50
    self.runLaplaceProblem() # we can use the same routine as for Laplace...

  def doChallenge1(self):
    m = MueLu_XMLChallengeMode()
    m.numProcs      = 1      # number of processors
    #m.globalNumDofs = 16641   # number of DOFs
    #m.nDofsPerNode  = 1      # DOFs per node
    m.solver        = "gmres"        # Belos solver
    m.tol           = 1e-12       # solver tolerance
    m.executable    = "./MueLu_TutorialDriver.exe" # executable
    m.problem       = "Recirc2D"   # string describing problem
    m.meshx = 129
    m.meshy = 129
    m.main()

  def doChallenge2(self):
    m = MueLu_XMLChallengeMode()
    m.numProcs      = 1      # number of processors
    #m.globalNumDofs = 7020   # number of DOFs
    #m.nDofsPerNode  = 2      # DOFs per node
    m.solver        = "cg"        # Belos solver
    m.tol           = 1e-12       # solver tolerance
    m.executable    = "./MueLu_TutorialDriver.exe" # executable
    m.problem       = "Elasticity2D"   # string describing problem
    m.meshx = 78
    m.meshy = 90
    m.main()

  def runLaplaceProblem(self):
    # check whether xml file exists

    while self.xmlFileName == "" or not os.path.isfile(self.xmlFileName) or not os.access(self.xmlFileName, os.R_OK):
      print(self.bcolors.FAIL+"Solver xml parameters: "+self.bcolors.ENDC + str(self.xmlFileName) + self.bcolors.FAIL + " invalid" + self.bcolors.ENDC)
      m = MueLu_XMLgenerator()
      m.askForSolver()
      self.xmlFileName = m.xmlFileName # store xml file

    while True:
      self.printActionMenu()

  def printActionMenu(self):
    #options = ['Rerun example', 'Show screen output', 'Change solver', 'Change processors', 'Exit']
    #callbacks = [self.runExample,self.printScreenOutput,self.changeSolver,self.changeProcs,self.doExitProgram]
    options = ['Rerun simulation', 'Show screen output', 'Change solver', 'Open xml file', 'Change procs', 'Change MG sweeps','Plot solution','Plot residual norm over ' + self.solver + ' solver iterations','Postprocess aggregates', 'Exit']
    callbacks = [self.runExample,self.printScreenOutput,self.changeSolver,self.openXMLfile,self.changeProcs, self.changeMGsweeps,self.plotSolution,self.doPlotResidual, self.postprocessAggregates, self.doExitProgram]
    while True:
      clearWindow()
      self.printSettings()
      print()
      if self.isDirty == True:
        print(self.bcolors.FAIL+ "DO NOT FORGET TO RUN THE EXAMPLE (option 0)" + self.bcolors.ENDC)
      else:
        print(self.bcolors.OKDARKGREEN + "Results up to date!" + self.bcolors.ENDC)
      print()
      self.runMenu(options,callbacks)

  def runExample(self):
    # runs example
    print("PREPARE SIMULATON")
    cmd = "rm *.vtp *.mat example*.txt output.log aggs*.txt nodes*.txt"
    runCommand(cmd)
    print("RUN EXAMPLE")
    cmd = "mpirun -np " + str(self.numprocs) + " " + str(self.executable) + " --matrixType=" + str(self.problem) + " --nx=" + str(self.meshx) + " --ny=" + str(self.meshy) + " --mgridSweeps=" + str(self.mgsweeps) + " --xml=" + str(self.xmlFileName) + " | tee output.log 2>&1"
    print(cmd)
    runCommand(cmd)
    runCommand("echo 'Press q to return.' >> output.log")
    print("POSTPROCESSING...")
    runCommand("cat example*.txt > example.txt")
    print("COMPLETE")
    self.isDirty = False
    waitForKey()

  def plotSolution(self):
    #cmd = "gnuplot -persist << _TTT_"
    #print(cmd)
    #runCommand(cmd)
    #cmd = "set dgrid3d " + str(self.meshy) + "," + str(self.meshx) + "\n set style data lines\n set nolabel \n set key off\n set autoscale\n splot " + "example.txt" + " using 3:4:5\n quit\n_TTT_"
    #runCommand(cmd)

    #proc1 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, text=True)
    self.proc1.stdin.write("set term x11 1\n")
    self.proc1.stdin.write("set title \"Exact solution\"\n")
    self.proc1.stdin.write("set dgrid3d " + str(self.meshy) + "," + str(self.meshx) + "\n")
    self.proc1.stdin.write("set style data lines\n")
    self.proc1.stdin.write("set nolabel\n")
    self.proc1.stdin.write("set key off\n")
    self.proc1.stdin.write("set autoscale\n")
    self.proc1.stdin.write("splot \"example.txt\" using 3:4:5\n")
    #self.proc1.stdin.write("quit\n") #close the gnuplot window
    self.proc1.stdin.flush()

    #proc2 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, text=True)
    self.proc2.stdin.write("set term x11 2\n") #wxt
    if (self.mgsweeps==1):
      self.proc2.stdin.write("set title \"Multigrid solution after " + str(self.mgsweeps) + " multigrid sweep\"\n")
    else:
      self.proc2.stdin.write("set title \"Multigrid solution after " + str(self.mgsweeps) + " multigrid sweeps\"\n")
    self.proc2.stdin.write("set dgrid3d " + str(self.meshy) + "," + str(self.meshx) + "\n")
    self.proc2.stdin.write("set style data lines\n")
    self.proc2.stdin.write("set nolabel\n")
    self.proc2.stdin.write("set key off\n")
    self.proc2.stdin.write("set autoscale\n")
    self.proc2.stdin.write("splot \"example.txt\" using 3:4:7\n")
    #self.proc2.stdin.write("quit\n") #close the gnuplot window
    self.proc2.stdin.flush()

    #proc3 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, text=True)
    self.proc3.stdin.write("set term x11 3\n")
    if (self.mgsweeps==1):
      self.proc3.stdin.write("set title \"Error (Exact vs. " + str(self.mgsweeps) + " multigrid sweep)\"\n")
    else:
      self.proc3.stdin.write("set title \"Error (Exact vs. " + str(self.mgsweeps) + " multigrid sweeps)\"\n")
    self.proc3.stdin.write("set dgrid3d " + str(self.meshy) + "," + str(self.meshx) + "\n")
    self.proc3.stdin.write("set style data lines\n")
    self.proc3.stdin.write("set palette model RGB defined ( 0 'black', 1 'white')\n")
    self.proc3.stdin.write("set nolabel\n")
    self.proc3.stdin.write("set key off\n")
    self.proc3.stdin.write("set autoscale\n")
    self.proc3.stdin.write("set hidden3d\n")
    self.proc3.stdin.write("set style line 1 lt 4 lw .5\n")
    self.proc3.stdin.write("set pm3d\n")
    self.proc3.stdin.write("splot \"example.txt\" using 3:4:($ 5-$ 7) with lines palette\n")
    #self.proc3.stdin.write("quit\n") #close the gnuplot window
    self.proc3.stdin.flush()

    #proc4 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, text=True)
    self.proc4.stdin.write("set term x11 4\n")
    self.proc4.stdin.write("set title \"Distribution of processors\"\n")
    self.proc4.stdin.write("set dgrid3d " + str(self.meshy) + "," + str(self.meshx) + "\n")
    self.proc4.stdin.write("set style data lines\n")
    self.proc4.stdin.write("set palette model RGB defined ( 0 'red', 1 'green', 2 'blue', 3 'yellow', 4 'pink')\n")
    self.proc4.stdin.write("set nolabel\n")
    self.proc4.stdin.write("set key off\n")
    self.proc4.stdin.write("set autoscale\n")
    self.proc4.stdin.write("set hidden3d\n")
    self.proc4.stdin.write("set style line 1 lt 4 lw .5\n")
    self.proc4.stdin.write("set pm3d\n")
    self.proc4.stdin.write("splot \"example.txt\" using 3:4:1:1 with points palette\n")
    #self.proc4.stdin.write("quit\n") #close the gnuplot window
    self.proc4.stdin.flush()

  def postprocessAggregates(self):
    # check whether "example.txt" is available
    if os.path.isfile("example.txt") == False:
      print(self.bcolors.FAIL+"Simulation data not available. Run the simulation first." + self.bcolors.ENDC)
      waitForKey()
      return

    if os.path.isfile("aggs_level0_proc0.out") == False:
      print(self.bcolors.FAIL+"No aggregation debug output found. Do not forget to turn on the AggregationExport factory in your xml file." + self.bcolors.ENDC)
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
    cmd = "chmod 750 ./MueLu_Agg2VTK.py"
    runCommand(cmd)
    cmd = "./MueLu_Agg2VTK.py"
    print(runCommand(cmd))

    if os.path.isfile("aggs0.vtp") == False:
      print(self.bcolors.WARNING+"Seems that the postprocessing failed (vtp files could not be created)." + self.bcolors.ENDC)
      waitForKey()
      return

    print(self.bcolors.OKDARKGREEN+"Use paraview to visualize generated vtk files for aggregates." + self.bcolors.ENDC)
    waitForKey()

  def printScreenOutput(self):
    clearWindow()
    if not os.path.isfile("output.log") or not os.access("output.log", os.R_OK):
      print(self.bcolors.FAIL+"Screen output not available."+self.bcolors.ENDC)
    else:
      print(runCommand("less output.log"))
    waitForKey()

  def openXMLfile(self):
    if self.editor in ["nano", "vim", "vi"]:
        subprocess.run([self.editor, self.xmlFileName])
    else:
        subprocess.Popen([self.editor + " " + self.xmlFileName], shell=True, stdin=subprocess.PIPE, text=True)

  def printProblemSelectionMenu(self):
    options = ['Laplace2D (50x50)', 'Laplace2D', 'Recirc2D (50x50)', 'Recirc2D', 'Challenge: Convection diffusion', 'Challenge: Elasticity problem', 'Exit']
    callbacks = [self.doLaplace2D50,self.doLaplace2Dn,self.doRecirc2D50,self.doRecirc2Dn,self.doChallenge1,self.doChallenge2, self.doExitProgram]
    while True:
      self.runMenu(options,callbacks)

  def changeSolver(self):
    self.xmlFileName = input("XML file name: ")
    self.isDirty = True
    while self.xmlFileName == "" or not os.path.isfile(self.xmlFileName) or not os.access(self.xmlFileName, os.R_OK):
      print(self.bcolors.FAIL+"Solver xml parameters: "+self.bcolors.ENDC + str(self.xmlFileName) + self.bcolors.FAIL + " invalid" + self.bcolors.ENDC)
      m = MueLu_XMLgenerator()
      m.xmlFileName=self.xmlFileName
      m.generateXMLfile()
      m.askForSolver()
      m.generateXMLfile()
      self.xmlFileName = m.xmlFileName # store xml file

  def changeProcs(self):
    self.numprocs = input("Number of processors: ")
    while not is_number(str(self.numprocs)):
      self.numprocs = input("Number of processors: ")
    self.isDirty = True

  def changeMGsweeps(self):
    self.mgsweeps = input("Number of Multigrid sweeps: ")
    while not is_number(str(self.mgsweeps)):
      self.mgsweeps = input("Number of Multigrid sweeps: ")
    self.isDirty = True

  def doPlotResidual(self):

    # prepare residual output file
    cmd = "grep iter: output.log > output.res"
    runCommand(cmd) 

    self.proc5.stdin.write("set term x11 1\n")
    self.proc5.stdin.write("set title \"Residual norm over " + str(self.solver) + " iterations\"\n")
    self.proc5.stdin.write("set style data lines\n")
    self.proc5.stdin.write("set xlabel \"# iterations\"\n")
    self.proc5.stdin.write("set ylabel \"Relative residual\"\n")
    self.proc5.stdin.write("set autoscale\n")
    self.proc5.stdin.write("set logscale y\n")
    self.proc5.stdin.write("plot '< sed \"s/[{}]//g\" output.res' using 2:5 w linespoints title \"" + str(self.xmlFileName) + "\"\n")
    self.proc5.stdin.flush()

  def printMainMenu(self):
    clearWindow()
    self.printSettings()
    print()
    print()
    while True:
      self.printProblemSelectionMenu()

  def doExitProgram(self):
    print("CLEAN UP temporary data")
    cmd = "rm *.vtp *.mat example*.txt output.log output.res aggs*.txt nodes*.txt"
    runCommand(cmd)
    print("QUIT")
    sys.exit()

  def printSettings(self):
    ## print out all made settings for xml file
    print(self.bcolors.HEADER+"***************************   PROBLEM   ****************************"+self.bcolors.ENDC)
    print(self.bcolors.WARNING+"Problem type:          "+self.bcolors.ENDC + str(self.problem))
    print(self.bcolors.WARNING+"Mesh:                  "+self.bcolors.ENDC + str(self.meshx) + "x" + str(self.meshy))
    print()
    if self.xmlFileName == "" or not os.path.isfile(self.xmlFileName) or not os.access(self.xmlFileName, os.R_OK):
      print(self.bcolors.FAIL+"Solver xml parameters: "+self.bcolors.ENDC + str(self.xmlFileName) + self.bcolors.FAIL + " invalid" + self.bcolors.ENDC)
    else:
      print(self.bcolors.WARNING+"Solver xml parameters:              "+self.bcolors.ENDC + str(self.xmlFileName))
    print(self.bcolors.WARNING+"Number of processors:               "+self.bcolors.ENDC + str(self.numprocs))
    print(self.bcolors.WARNING+"Number of Multigrid solving sweeps: "+self.bcolors.ENDC + str(self.mgsweeps))
    print(self.bcolors.HEADER+"***************************   PROBLEM   ****************************"+self.bcolors.ENDC)

class MueLu_XMLChallengeMode():
  """ Menu and options for challenge mode """

  def __init__(self):

    self.numProcs      = 1      # number of processors
    self.globalNumDofs = 7020   # number of DOFs
    self.nDofsPerNode  = 2      # DOFs per node
    self.solver        = "cg"        # Belos solver
    self.tol           = 1e-12       # solver tolerance
    self.executable    = "./MueLu_TutorialDriver.exe" # executable
    self.problem       = "Elasticity2D"   # string describing problem
    self.xmlReferenceFileName = ""
    self.xmlFileName   = ""
    self.bcolors       = usecolors()
    self.isDirty       = True  # dirty flag
    self.editor = os.environ.get("H_EDITOR", "nano") # text editor is specified by the env var H_EDITOR, nano by default
    self.has_coords    = False ### never used?

    self.proc1 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, text=True)

  def main(self):

    # check if tar.gz file with data is in subfolder challenges:
#    if os.path.isfile("challenges/" + self.problem + ".tar.gz") == False:
#      cmd = "rm -Rf challenges"
#      runCommand(cmd)
#      print("Download additional files")
#      print(self.bcolors.WARNING+"https://trilinos.org/wordpress/wp-content/uploads/2015/07/MueLu_tutorial_challenges.tar.gz"+self.bcolors.ENDC)
#      cmd = "wget --no-check-certificate https://trilinos.org/wordpress/wp-content/uploads/2015/07/MueLu_tutorial_challenges.tar.gz"
#      runCommand(cmd)
#      print("Extract files...")
#      cmd = "tar xvf MueLu_tutorial_challenges.tar.gz"
#      runCommand(cmd)
#      print(self.bcolors.OKDARKGREEN + "Success!" + self.bcolors.ENDC)

    # generate results for reference xml files
    self.xmlReferenceFileName = "challenges/" + self.problem + "_reference.xml"

    # copy file with reference parameters for this example
    cmd = "cp challenges/" + self.problem + "_reference.xml " + self.problem + "_parameters.xml"
    runCommand(cmd)

    self.xmlFileName   = self.problem + "_parameters.xml"     # xml parameter file

    if os.path.isfile("challenges/" + self.problem + "_coords.txt"):
      self.has_coords = True

    self.doRunReference()

    while True:
      self.printMainMenu()

  def runMenu(self,options,callbacks):
    for i,option in enumerate(options):
      print('%s. %s' % (i, option)) # display all options
    choice = input('your choice? ')
    if is_number(str(choice)):
      if int(choice) < len(options) and int(choice) > -1:
        callbacks[int(choice)]() # call corresponding function
      else:
        print("Option is out of range!")
    else:
      print("Invalid option! " + str(choice))

  # print main menu for challenge mode
  def printMainMenu(self):
    clearWindow()
    self.printSettings()
    print()
    if self.isDirty == True:
      print(self.bcolors.FAIL+ "DO NOT FORGET TO RUN THE EXAMPLE (option 0)" + self.bcolors.ENDC)
    else:
      print(self.bcolors.OKDARKGREEN + "Results up to date!" + self.bcolors.ENDC)
      print()
      self.printResults()
    print()

    options = ['Run example','Show screen output', 'Change XML parameter file', 'Open xml file', 'Change procs', 'Change linear solver', 'Plot residual', 'Exit']
    callbacks = [self.doRunExample,self.printScreenOutput,self.changeSolver,self.openXMLfile,self.changeProcs, self.doSolverMenu,self.doPlotResidual, self.doExitProgram]

    self.runMenu(options,callbacks)

  def printSettings(self):
    ## print out all made settings for xml file
    print(self.bcolors.HEADER+"***************************   PROBLEM   ****************************"+self.bcolors.ENDC)
    print(self.bcolors.WARNING+"Problem type:          "+self.bcolors.ENDC + str(self.problem))
    print(self.bcolors.WARNING+"Problem size:          "+self.bcolors.ENDC + str(self.globalNumDofs))
    print()
    if self.xmlFileName == "" or not os.path.isfile(self.xmlFileName) or not os.access(self.xmlFileName, os.R_OK):
      print(self.bcolors.FAIL+"Solver xml parameters: "+self.bcolors.ENDC + str(self.xmlFileName) + self.bcolors.FAIL + " invalid" + self.bcolors.ENDC)
    else:
      print(self.bcolors.WARNING+"Solver xml parameters:              "+self.bcolors.ENDC + str(self.xmlFileName))
    print(self.bcolors.WARNING+"Number of processors:               "+self.bcolors.ENDC + str(self.numProcs))
    print(self.bcolors.WARNING+"Solver (Tolerance): "+self.bcolors.ENDC + str(self.solver) + " (" + str(self.tol) + ")")
    print(self.bcolors.HEADER+"***************************   PROBLEM   ****************************"+self.bcolors.ENDC)

  def printResults(self):
    cmd = "grep 'total iterations:' output.log"
    iter = runCommand(cmd)
    cmd = "grep 'Solution time:' output.log"
    time = runCommand(cmd)
    cmd = "grep 'total iterations:' reference.log"
    refiter = runCommand(cmd)
    cmd = "grep 'Solution time:' reference.log"
    reftime = runCommand(cmd)
    print(self.bcolors.HEADER+"***************************   RESULTS   ****************************"+self.bcolors.ENDC)
    print("Reference settings:")
    print(str(refiter))
    print(str(reftime))
    print("Your settings:")
    print(str(iter))
    print(str(time))
    print(self.bcolors.HEADER+"***************************   RESULTS   ****************************"+self.bcolors.ENDC)

  def doRunReference(self):
    print("Please wait...")
    # runs example
    cmd = "rm -f *.vtp *.mat example*.txt output.log output.res reference.log reference.res aggs*.txt nodes*.txt"
    runCommand(cmd)
    #cmd = "mpirun -np " + str(self.numProcs) + " " + str(self.executable) + " --globalNumDofs=" + str(self.globalNumDofs) + " --nDofsPerNode=" + str(self.nDofsPerNode) + " --solver=" + str(self.solver) + " --tol=" + str(self.tol) + " --xml=" + self.xmlReferenceFileName + " --problem=challenges/" + str(self.problem) + " --coordinates=challenges/" + str(self.problem) + "_coords.txt" + " | tee reference.log 2>&1"
    cmd = "mpirun -np " + str(self.numProcs) + " " + str(self.executable) + " --matrixType=" + str(self.problem) + " --nx=" + str(self.meshx) + " --ny=" + str(self.meshy) + " --solver=" + str(self.solver) + " --tol=" + str(self.tol) + " --xml=" + self.xmlReferenceFileName + " | tee reference.log 2>&1"
    runCommand(cmd)
    self.isDirty = False


  def doRunExample(self):
    # runs example
    print("PREPARE SIMULATON")
    cmd = "rm -f *.vtp *.mat example*.txt output.log output.res aggs*.txt nodes*.txt"
    runCommand(cmd)
    print("RUN EXAMPLE")
    #cmd = "mpirun -np " + str(self.numProcs) + " " + str(self.executable) + " --globalNumDofs=" + str(self.globalNumDofs) + " --nDofsPerNode=" + str(self.nDofsPerNode) + " --solver=" + str(self.solver) + " --tol=" + str(self.tol) + " --xml=" + self.xmlFileName + " --problem=challenges/" + str(self.problem) + " --coordinates=challenges/" + str(self.problem) + "_coords.txt" + " | tee output.log 2>&1"
    cmd = "mpirun -np " + str(self.numProcs) + " " + str(self.executable) + " --matrixType=" + str(self.problem) + " --nx=" + str(self.meshx) + " --ny=" + str(self.meshy) + " --solver=" + str(self.solver) + " --tol=" + str(self.tol) + " --xml=" + self.xmlFileName + " | tee output.log 2>&1"
    runCommand(cmd)
    print(cmd)
    runCommand(cmd)
    runCommand("echo 'Press q to return.' >> output.log")
    print("POSTPROCESSING...")
    runCommand("cat example*.txt > example.txt")
    print("COMPLETE")
    self.isDirty = False
    waitForKey()

  def printScreenOutput(self):
    clearWindow()
    if not os.path.isfile("output.log") or not os.access("output.log", os.R_OK):
      print(self.bcolors.FAIL+"Screen output not available."+self.bcolors.ENDC)
    else:
      print(runCommand("less output.log"))
    waitForKey()

  def changeSolver(self):
    self.xmlFileName = input("XML file name: ")
    self.isDirty = True
    while self.xmlFileName == "" or not os.path.isfile(self.xmlFileName) or not os.access(self.xmlFileName, os.R_OK):
      print(self.bcolors.FAIL+"Solver xml parameters: "+self.bcolors.ENDC + str(self.xmlFileName) + self.bcolors.FAIL + " invalid" + self.bcolors.ENDC)
      m = MueLu_XMLgenerator()
      m.xmlFileName=self.xmlFileName
      m.generateXMLfile()
      m.askForSolver()
      m.generateXMLfile()
      self.xmlFileName = m.xmlFileName # store xml file

  def openXMLfile(self):
    if self.editor in ["nano", "vim", "vi"]:
        subprocess.run([self.editor, self.xmlFileName])
    else:
        subprocess.Popen([self.editor + " " + self.xmlFileName], shell=True, stdin=subprocess.PIPE, text=True)

  def changeProcs(self):
    self.numProcs = input("Number of processors: ")
    while not is_number(str(self.numProcs)):
      self.numProcs = input("Number of processors: ")
    self.isDirty = True
    self.doRunReference()

  def doCGIteration(self):
    self.solver = "cg"
    self.isDirty = True
    self.doRunReference()
  def doGMRESIteration(self):
    self.solver = "gmres"
    self.isDirty = True
    self.doRunReference()

  def doSolverMenu(self):
    options = ['CG method', 'GMRES method']
    callbacks = [self.doCGIteration,self.doGMRESIteration]
    #while self.exitLoop == False:
    self.runMenu(options,callbacks)
    #self.exitLoop=True #False

  def doPlotResidual(self):

    # prepare residual output file
    cmd = "grep iter: output.log > output.res"
    runCommand(cmd)

    # prepare reference data
    cmd = "grep iter: reference.log > reference.res"
    runCommand(cmd)
    self.proc1.stdin.write("set term x11 1\n")
    self.proc1.stdin.write("set title \"Residual norm over " + str(self.solver) + " iterations\"\n")
    self.proc1.stdin.write("set style data lines\n")
    self.proc1.stdin.write("set xlabel \"# iterations\"\n")
    self.proc1.stdin.write("set ylabel \"Relative residual\"\n")
    self.proc1.stdin.write("set autoscale\n")
    self.proc1.stdin.write("set logscale y\n")
    self.proc1.stdin.write("plot '< sed \"s/[{}]//g\" output.res' using 2:5 w linespoints title \"" + str(self.xmlFileName) + "\"\n")

    self.proc1.stdin.flush()

  def doExitProgram(self):
    runCommand("rm output.log reference.log output.res reference.res")
    sys.exit() # terminate full program


    # gnuplot commands
    # set logscale y
    # plot "temp.log" using 5 w linespoints

class MueLu_XMLgenerator():
  """Simple generator for MueLu xml files."""

  def __init__(self):
    # common MG settings
    self.xmlFileName = ""    # output file name for xml data
    self.maxMultLevels = 5   # maximum number of levels
    self.maxCoarseSize = 1000 # max. coarse size

    # aggregate settings
    self.dropTolerance = 0.0
    self.minAggSize = 4
    self.maxAggSize = 9
    self.maxNeighCount = 0

    # smoother settings
    self.levelSmoother = "Jacobi"
    self.levelSmootherSweeps = 1
    self.levelSmootherDamp   = 0.7
    self.coarseSolver = "Direct"

    # transfer operators
    self.transferOps = "PA-AMG"
    self.transferOpDamp = 1.33

    # restriction operators
    self.restrictionOp = "TransPFactory"

    # rebalancing
    self.doRebalancing = False
    self.minRowsPerProc = 800
    self.nnzImbalance = 1.1
    self.rebStartLevel = 1

    self.bcolors = usecolors()
    self.isDirty = True                   # flag to store, whether changes have been saved or not
    self.exitLoop = False                 # set to true to exit current loop

    print(self.bcolors.FAIL+'===================================================================================='+self.bcolors.ENDC)
    print('===================================================================================='+self.bcolors.ENDC)




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
    self.xmlFileName = input("XML file name: ")
    self.isDirty = True

  def doRelaxationMaxLevels(self):
    self.maxMultLevels = input("Max. multigrid levels: ")
    self.isDirty = True

  def doRelaxationMaxCoarseSize(self):
    self.maxCoarseSize = input("Max. coarse size: ")
    self.isDirty = True

  def doRelaxationJacobi(self):
    self.levelSmoother = "Jacobi"
    self.levelSmootherSweeps = input("Smoother sweeps: ")
    self.levelSmootherDamp   = input("Smoother damping: ")
    self.isDirty = True

  def doRelaxationGS(self):
    self.levelSmoother = "Gauss-Seidel"
    self.levelSmootherSweeps = input("Smoother sweeps: ")
    self.levelSmootherDamp   = input("Smoother damping: ")
    self.isDirty = True

  def doRelaxationSymGS(self):
    self.levelSmoother = "Sym.Gauss-Seidel"
    self.levelSmootherSweeps = input("Smoother sweeps: ")
    self.levelSmootherDamp   = input("Smoother damping: ")
    self.isDirty = True

  def doDropTolerance(self):
    self.dropTolerance = input("Drop tolerance for matrix graph (default = 0.0): ")
    self.isDirty = True
  def doMinAggSize(self):
    self.minAggSize = input("Minimum number of nodes per aggregate: ")
    self.isDirty
  def doMaxAggSize(self):
    self.maxAggSize = input("Maximum number of nodes per aggregate: ")
    self.isDirty
  def doMaxNeigh(self):
    self.maxNeighCount = input("Maximum number of already aggregated neighbor nodes (default = 0): ")
    self.isDirty

  # Transfer operators
  def doPaAMG(self):
    self.transferOps = "PA-AMG"
    self.transferOpDamp = 0.0
    if self.restrictionOp == "GenericRFactory":
      self.restrictionOp = "TransPFactory"
      print(self.bcolors.WARNING + "GenericRFactory cannot be used with non-smoothed PA-AMG prolongation operators. We change it back to TransPFactory."+self.bcolors.ENDC)
      print()
      print("Press any key to proceed")
      waitForKey()

    self.isDirty = True
  def doSaAMG(self):
    self.transferOps = "SA-AMG"
    self.transferOpDamp = input("Transfer operator damping: ")
    self.isDirty = True
  def doPgAMG(self):
    self.transferOps = "PG-AMG"
    self.transferOpDamp = 0.0
    self.isDirty = True

  # Restriction operators
  def doSymR(self):
    self.restrictionOp = "TransPFactory"
    self.isDirty = True
  def doNonsymR(self):
    self.restrictionOp = "GenericRFactory"
    if self.transferOps == "PA-AMG":
      self.restrictionOp = "TransPFactory"
      print(self.bcolors.WARNING+"GenericRFactory cannot be used with non-smoothed PA-AMG prolongation operators. We change it back to TransPFactory.")
      print("To use GenericRFactory you have to select either SaPFactory or PgPFactory for prolongation."+self.bcolors.ENDC)
      print()
      print("Press any key to proceed")
      waitForKey()
    self.isDirty = True

  # Rebalancing
  def doRebalancingOption(self):
    self.doRebalancing = True
    self.minRowsPerProc = input("Minimum number of DOFs per processor: ")
    self.nnzImbalance = input("Max. nonzero imbalance (default 1.1): ")
    self.rebStartLevel = input("Start rebalancing on level (default 1): ")
    self.isDirty = True

  def doNoRebalancingOption(self):
    self.doRebalancing = False
    self.isDirty = True

  def SetExitLoop(self):
    self.exitLoop = True

  def runMenu(self,options,callbacks):
    for i,option in enumerate(options):
      print('%s. %s' % (i, option)) # display all options
    choice = input('your choice? ')
    if is_number(str(choice)):
      if int(choice) < len(options) and int(choice) > -1:
        callbacks[int(choice)]() # call corresponding function
      else:
        print("Option is out of range!")
    else:
      print("Invalid option! " + str(choice))

  def doCommonMenu(self):
    options = ['Max. multigrid levels', 'Max. coarse size', 'Back']
    callbacks = [self.doRelaxationMaxLevels,self.doRelaxationMaxCoarseSize, self.askForSolver]
    while self.exitLoop == False:
      self.runMenu(options,callbacks)
    self.exitLoop=True #False

  def doAggregatesMenu(self):
    options = ['Drop tolerance', 'Min. aggregate size', 'Max. aggregate size', 'Max. Neighbor Count', 'Back']
    callbacks = [self.doDropTolerance,self.doMinAggSize, self.doMaxAggSize, self.doMaxNeigh, self.askForSolver]
    while self.exitLoop == False:
      self.runMenu(options,callbacks)
    self.exitLoop=True #False

  def doSmootherMenu(self):
    options = ['Jacobi', 'Gauss-Seidel', 'Sym. Gauss-Seidel', 'Back']
    callbacks = [self.doRelaxationJacobi,self.doRelaxationGS, self.doRelaxationSymGS, self.askForSolver]
    self.runMenu(options,callbacks)

  def doTransferMenu(self):
    options = ['Non-smoothed transfer (PA-AMG)', 'Smoothed transfer (SA-AMG)', 'Smoothed transfer (PG-AMG)', 'Back']
    callbacks = [self.doPaAMG,self.doSaAMG, self.doPgAMG, self.askForSolver]
    self.runMenu(options,callbacks)

  def doRestrictorMenu(self):
    options = ['Symmetric', 'Non-symmetric', 'Back']
    callbacks = [self.doSymR,self.doNonsymR, self.askForSolver]
    self.runMenu(options,callbacks)

  def doRebalancingMenu(self):
    options = ['No rebalancing', 'Activate rebalancing', 'Back']
    callbacks = [self.doNoRebalancingOption,self.doRebalancingOption, self.askForSolver]
    self.runMenu(options,callbacks)

  def doExitProgram(self):
    #sys.exit()
    print("doEXIT")
    self.exitLoop = True

  def printMainMenu(self):
    clearWindow()
    self.printSettings()
    print()
    print()

    #options = ['Set Output file name','Common Multigrid settings', 'Level smoother settings', 'Transfer operators', 'Restriction operators', 'Save XML file', 'Exit']
    #callbacks = [self.doFileName, self.doCommonMenu, self.doSmootherMenu, self.doTransferMenu, self.doRestrictorMenu, self.generateXMLfile, self.doExitProgram]
    options = ['Common Multigrid settings', 'Aggregate settings', 'Level smoother settings', 'Transfer operators', 'Restriction operators', 'Rebalancing options', 'Save XML file', 'Back']
    callbacks = [self.doCommonMenu, self.doAggregatesMenu, self.doSmootherMenu, self.doTransferMenu, self.doRestrictorMenu, self.doRebalancingMenu, self.generateXMLfile, self.doExitProgram]

    self.runMenu(options,callbacks)

  def printSettings(self):
    ## print out all made settings for xml file
    print(self.bcolors.HEADER+"***************************   SETTINGS   ****************************"+self.bcolors.ENDC)
    print(self.bcolors.WARNING+"XML file name:           "+self.bcolors.ENDC + str(self.xmlFileName))
    print()
    print(self.bcolors.WARNING+"Max. MultiGrid levels:   "+self.bcolors.ENDC + str(self.maxMultLevels))
    print(self.bcolors.WARNING+"Max. CoarseSize:         "+self.bcolors.ENDC + str(self.maxCoarseSize))
    print()
    print(self.bcolors.WARNING+"Level smoother:          "+self.bcolors.ENDC + str(self.levelSmoother))
    print(self.bcolors.WARNING+"Level smoothing sweeps:  "+self.bcolors.ENDC + str(self.levelSmootherSweeps))
    print(self.bcolors.WARNING+"Level damping parameter: "+self.bcolors.ENDC + str(self.levelSmootherDamp))
    print()
    print(self.bcolors.WARNING+"Coarse solver:           "+self.bcolors.ENDC + str(self.coarseSolver))
    print()
    print(self.bcolors.WARNING+"Graph drop tolerance:    "+self.bcolors.ENDC + str(self.dropTolerance))
    print(self.bcolors.WARNING+"Aggregate size (min/max):"+self.bcolors.ENDC + str(self.minAggSize) + "/" + str(self.maxAggSize))
    print(self.bcolors.WARNING+"Max. neighbor count:     "+self.bcolors.ENDC + str(self.maxNeighCount))
    print()
    print(self.bcolors.WARNING+"Transfer operators:     "+self.bcolors.ENDC + str(self.transferOps))
    print(self.bcolors.WARNING+"Transfer smoothing par.:"+self.bcolors.ENDC + str(self.transferOpDamp))
    print()
    print(self.bcolors.WARNING+"Restriction operator:   "+self.bcolors.ENDC + str(self.restrictionOp))
    print()
    if self.doRebalancing == False:
      print(self.bcolors.WARNING+"NO Rebalancing"+self.bcolors.ENDC)
    else:
      print(self.bcolors.WARNING+"Rebalancing active:"+ self.bcolors.ENDC)
      print(self.bcolors.WARNING+"Minimum DOFs per proc:  "+ self.bcolors.ENDC + str(self.minRowsPerProc))
      print(self.bcolors.WARNING+"Nonzero imbalance:      "+ self.bcolors.ENDC + str(self.nnzImbalance))
      print(self.bcolors.WARNING+"Start level for rebal.: "+ self.bcolors.ENDC + str(self.rebStartLevel))
    print(self.bcolors.HEADER+"***************************   SETTINGS   ****************************"+self.bcolors.ENDC)

    print()
    if self.isDirty == True:
      print(self.bcolors.FAIL+ "CHANGES HAVE NOT BEEN SAVED!" + self.bcolors.ENDC)
    else:
      print(self.bcolors.OKDARKGREEN + "CHANGES HAVE BEEN SAVED!" + self.bcolors.ENDC)

  def generateXMLfile(self):
    # generate HEAD file for pre_exodus
    if os.path.isfile(self.xmlFileName):
      os.remove(self.xmlFileName)
    o = open(self.xmlFileName,"a")
    for line in open("tmpl/muelu.xml_TMPL"):
    #for line in open(headfile_tmpl):
      line = line.replace("$SMOO_SWEEPS", str(self.levelSmootherSweeps))
      line = line.replace("$SMOO_DAMP"  , str(self.levelSmootherDamp))
      line = line.replace("$SMOOTHER"   , str(self.levelSmoother))
      line = line.replace("$MAXLEVELS"  , str(self.maxMultLevels))
      line = line.replace("$MAXCOARSESIZE", str(self.maxCoarseSize))
      line = line.replace("$RESTRICTOR",  str(self.restrictionOp))
      line = line.replace("$PROLONGATOR", str(self.transferOps))
      line = line.replace("$SADAMPING"  , str(self.transferOpDamp))
      line = line.replace("$DROPTOL"    , str(self.dropTolerance))
      line = line.replace("$MAXNEIGH"    , str(self.maxNeighCount))
      line = line.replace("$MINAGGS"    , str(self.minAggSize))
      line = line.replace("$MAXAGGS"    , str(self.maxAggSize))

      if self.doRebalancing == False:
        line = line.replace("$MANAGER_PROLONGATOR", str(self.transferOps))
        line = line.replace("$MANAGER_RESTRICTOR",  "myRestrictorFact")
        line = line.replace("$MANAGER_RAP", "myRAPFact")
        line = line.replace("$MANAGER_NULLSPACE", "PA-AMG")
      else:
        line = line.replace("$MANAGER_PROLONGATOR", "myRebalanceProlongatorFact")
        line = line.replace("$MANAGER_RESTRICTOR",  "myRebalanceRestrictionFact")
        line = line.replace("$MANAGER_RAP", "myRebalanceAFact")
        line = line.replace("$MANAGER_NULLSPACE", "myRebalanceProlongatorFact")
      o.write(line)
    o.close()
    self.isDirty = False

if __name__ == '__main__':
  #MueLu_XMLgenerator().main()
  ProblemHandler().main()