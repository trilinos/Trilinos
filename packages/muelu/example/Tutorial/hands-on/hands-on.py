#!/usr/bin/python
import os
import sys
import math
import subprocess

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
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    OKDARKGREEN = '\033[32m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

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
    self.problem    = "Laplace 2D"
    self.executable = "MueLu_laplace2d.exe"
    self.meshx      = 50
    self.meshy      = 50
    self.mgsweeps   = 1
    self.numprocs   = 2
    self.xmlFileName = "xml/s2a.xml"
    
    self.proc1 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, )
    self.proc2 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, )
    self.proc3 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, )
    self.proc4 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, )
    
    self.isDirty = True                   # flag to store, whether problem has to be rerun or not

    self.editor = "kwrite"    # TODO replace me by local editor...
    
  def main(self):
    self.printMainMenu()
    
  def runMenu(self,options,callbacks):
    for i,option in enumerate(options):
      print('%s. %s' % (i, option)) # display all options
    choice = raw_input('your choice? ')
    if is_number(str(choice)) and int(choice) < len(options):
      callbacks[int(choice)]() # call correspondending function
    else:
      print "ups: choice = " + str(choice) + " len(option)=" + str(len(option))



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
      print bcolors.FAIL+"Solver xml parameters: "+bcolors.ENDC + str(self.xmlFileName) + bcolors.FAIL + " invalid" + bcolors.ENDC
      m = MueLu_XMLgenerator()
      m.askForSolver()
      self.xmlFileName = m.xmlFileName # store xml file
      
    while True:
      self.printActionMenu()
      
  def printActionMenu(self):
    #options = ['Rerun example', 'Show screen output', 'Change solver', 'Change processors', 'Exit']
    #callbacks = [self.runExample,self.printScreenOutput,self.changeSolver,self.changeProcs,self.doExitProgram]
    options = ['Rerun simulation', 'Show screen output', 'Change solver', 'Open xml file', 'Change procs', 'Change MG sweeps','Plot solution','Postprocess aggregates', 'Exit']
    callbacks = [self.runExample,self.printScreenOutput,self.changeSolver,self.openXMLfile,self.changeProcs, self.changeMGsweeps,self.plotSolution, self.postprocessAggregates, self.doExitProgram]
    while True:
      clearWindow()    
      self.printSettings()
      print ""
      if self.isDirty == True:
        print bcolors.FAIL+ "DO NOT FORGET TO RUN THE EXAMPLE (option 0)" + bcolors.ENDC
      else:
        print bcolors.OKDARKGREEN + "Results up to date!" + bcolors.ENDC   
      print ""
      self.runMenu(options,callbacks)
      
  def runExample(self):
    # runs example
    print "PREPARE SIMULATON"
    cmd = "rm *.vtp *.mat example*.txt output.log aggs*.txt nodes*.txt"
    runCommand(cmd)
    print "RUN EXAMPLE"
    cmd = "mpirun -np " + str(self.numprocs) + " " + str(self.executable) + " --nx=" + str(self.meshx) + " --ny=" + str(self.meshy) + " --mgridSweeps=" + str(self.mgsweeps) + " --xml=" + str(self.xmlFileName) + " | tee output.log 2>&1"
    print cmd
    runCommand(cmd)
    runCommand("echo 'Press q to return.' >> output.log")
    print "POSTPROCESSING..."
    runCommand("cat example*.txt > example.txt")
    print "COMPLETE"
    self.isDirty = False
    waitForKey()
  
  def plotSolution(self):
    #cmd = "gnuplot -persist << _TTT_"
    #print cmd
    #runCommand(cmd)
    #cmd = "set dgrid3d " + str(self.meshy) + "," + str(self.meshx) + "\n set style data lines\n set nolabel \n set key off\n set autoscale\n splot " + "example.txt" + " using 3:4:5\n quit\n_TTT_"
    #runCommand(cmd)
    
    #proc1 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, )
    self.proc1.stdin.write("set term x11 1\n")
    self.proc1.stdin.write("set title \"Solution\"\n")
    self.proc1.stdin.write("set dgrid3d " + str(self.meshy) + "," + str(self.meshx) + "\n")
    self.proc1.stdin.write("set style data lines\n")
    self.proc1.stdin.write("set nolabel\n")
    self.proc1.stdin.write("set key off\n")
    self.proc1.stdin.write("set autoscale\n")
    self.proc1.stdin.write("splot \"example.txt\" using 3:4:5\n")
    #self.proc1.stdin.write("quit\n") #close the gnuplot window
    self.proc1.stdin.flush()

    #proc2 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, )
    self.proc2.stdin.write("set term x11 2\n") #wxt
    self.proc2.stdin.write("set title \"Multigrid solution\"\n")    
    self.proc2.stdin.write("set dgrid3d " + str(self.meshy) + "," + str(self.meshx) + "\n")
    self.proc2.stdin.write("set style data lines\n")
    self.proc2.stdin.write("set nolabel\n")
    self.proc2.stdin.write("set key off\n")
    self.proc2.stdin.write("set autoscale\n")
    self.proc2.stdin.write("splot \"example.txt\" using 3:4:7\n")
    #self.proc2.stdin.write("quit\n") #close the gnuplot window
    self.proc2.stdin.flush()
    
    #proc3 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, )
    self.proc3.stdin.write("set term x11 3\n")
    self.proc3.stdin.write("set title \"Error (Exact vs. Multigrid)\"\n")    
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
    
    #proc4 = subprocess.Popen(['gnuplot','-p'], shell=True, stdin=subprocess.PIPE, )
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
    self.proc4.stdin.write("splot \"example.txt\" using 3:4:1 with points palette\n")
    #self.proc3.stdin.write("quit\n") #close the gnuplot window
    self.proc3.stdin.flush()
    
  def postprocessAggregates(self): 
    # check whether "example.txt" is available
    if os.path.isfile("example.txt") == False:
      print bcolors.FAIL+"Simulation data not available. Run the simulation first." + bcolors.ENDC
      waitForKey()
      return
      
    if os.path.isfile("aggs_level0_proc0.txt") == False:
      print bcolors.FAIL+"No aggregation debug output found. Do not forget to turn on the AggregationExport factory in your xml file." + bcolors.ENDC
      waitForKey()
      return
    
    if os.path.isfile("MueLu_Agg2VTK.py"):
      os.remove("MueLu_Agg2VTK.py")
    o = open("MueLu_Agg2VTK.py","a")
    for line in open("tmpl/MueLu_Agg2VTK.py_TMPL"):
      line = line.replace("$NUMPROCS", str(self.numprocs))     
      o.write(line)
    o.close()    
    
    print "POSTPROCESS AGGREGATION OUTPUT DATA"
    cmd = "chmod 750 ./MueLu_Agg2VTK.py"
    runCommand(cmd)
    cmd = "./MueLu_Agg2VTK.py"
    print runCommand(cmd)
    
    if os.path.isfile("aggs0.vtp") == False:
      print bcolors.WARNING+"Seems that the postprocessing failed (vtp files could not be created)." + bcolors.ENDC
      waitForKey()
      return
    
    print bcolors.OKDARKGREEN+"Use paraview to visualize generated vtk files for aggregates." + bcolors.ENDC
    waitForKey()
  
  def printScreenOutput(self):
    clearWindow()
    if not os.path.isfile("output.log") or not os.access("output.log", os.R_OK):
      print bcolors.FAIL+"Screen output not available."+bcolors.ENDC
    else: 
      print runCommand("less output.log")
    waitForKey()
    
  def openXMLfile(self):
    editor = subprocess.Popen([self.editor + " " + self.xmlFileName], shell=True, stdin=subprocess.PIPE, )

  def printProblemSelectionMenu(self):
    options = ['Laplace 2D (50x50)', 'Laplace 2D', 'Recirc 2D (50x50)', 'Recirc 2D', 'Exit']
    callbacks = [self.doLaplace2D50,self.doLaplace2Dn,self.doRecirc2D50,self.doRecirc2Dn, self.doExitProgram]
    while True:
      self.runMenu(options,callbacks)
  
  def changeSolver(self):
    self.xmlFileName = raw_input("XML file name: ")
    self.isDirty = True
    while self.xmlFileName == "" or not os.path.isfile(self.xmlFileName) or not os.access(self.xmlFileName, os.R_OK):
      print bcolors.FAIL+"Solver xml parameters: "+bcolors.ENDC + str(self.xmlFileName) + bcolors.FAIL + " invalid" + bcolors.ENDC
      m = MueLu_XMLgenerator()
      m.xmlFileName=self.xmlFileName
      m.generateXMLfile()
      m.askForSolver()
      m.generateXMLfile()
      self.xmlFileName = m.xmlFileName # store xml file  
    
  def changeProcs(self):
    self.numprocs = raw_input("Number of processors: ")
    while not is_number(str(self.numprocs)):
      self.numprocs = raw_input("Number of processors: ")
    self.isDirty = True
  
  def changeMGsweeps(self):
    self.mgsweeps = raw_input("Number of Multigrid sweeps: ")
    while not is_number(str(self.mgsweeps)):
      self.mgsweeps = raw_input("Number of Multigrid sweeps: ")
    self.isDirty = True
		       
  def printMainMenu(self):
    clearWindow()
    self.printSettings()
    print ""
    print ""    
    while True:
      self.printProblemSelectionMenu()
      
  def doExitProgram(self):
    print "CLEAN UP temporary data"
    cmd = "rm *.vtp *.mat example*.txt output.log aggs*.txt nodes*.txt"
    runCommand(cmd)
    print "QUIT"
    sys.exit()  
 
  def printSettings(self):
    ## print out all made settings for xml file
    print bcolors.HEADER+"***************************   PROBLEM   ****************************"+bcolors.ENDC
    print bcolors.WARNING+"Problem type:          "+bcolors.ENDC + str(self.problem)
    print bcolors.WARNING+"Mesh:                  "+bcolors.ENDC + str(self.meshx) + "x" + str(self.meshy)
    print ""
    if self.xmlFileName == "" or not os.path.isfile(self.xmlFileName) or not os.access(self.xmlFileName, os.R_OK):
      print bcolors.FAIL+"Solver xml parameters: "+bcolors.ENDC + str(self.xmlFileName) + bcolors.FAIL + " invalid" + bcolors.ENDC
    else:
      print bcolors.WARNING+"Solver xml parameters:              "+bcolors.ENDC + str(self.xmlFileName)
    print bcolors.WARNING+"Number of processors:               "+bcolors.ENDC + str(self.numprocs)
    print bcolors.WARNING+"Number of Multigrid solving sweeps: "+bcolors.ENDC + str(self.mgsweeps)
    print bcolors.HEADER+"***************************   PROBLEM   ****************************"+bcolors.ENDC

      
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
    
    self.isDirty = True                   # flag to store, whether changes have been saved or not
    self.exitLoop = False                 # set to true to exit current loop
    
    print bcolors.FAIL+'===================================================================================='+bcolors.ENDC
    print '===================================================================================='+bcolors.ENDC

    
    
    
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
    
  def doRelaxationMaxLevels(self):
    self.maxMultLevels = raw_input("Max. multigrid levels: ")
    self.isDirty = True

  def doRelaxationMaxCoarseSize(self):
    self.maxCoarseSize = raw_input("Max. coarse size: ")
    self.isDirty = True
    
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
    self.isDirty = True
  def doMinAggSize(self):
    self.minAggSize = raw_input("Minimum number of nodes per aggregate: ")
    self.isDirty
  def doMaxAggSize(self):
    self.maxAggSize = raw_input("Maximum number of nodes per aggregate: ")
    self.isDirty
  def doMaxNeigh(self):
    self.maxNeighCount = raw_input("Maximum number of already aggregated neighbor nodes (default = 0): ")
    self.isDirty
    
  # Transfer operators
  def doPaAMG(self):
    self.transferOps = "PA-AMG"
    self.transferOpDamp = 0.0
    if self.restrictionOp == "GenericRFactory":
      self.restrictionOp = "TransPFactory"
      print bcolors.WARNING + "GenericRFactory cannot be used with non-smoothed PA-AMG prolongation operators. We change it back to TransPFactory."+bcolors.ENDC
      print ""
      print "Press any key to proceed"
      waitForKey()

    self.isDirty = True
  def doSaAMG(self):
    self.transferOps = "SA-AMG"
    self.transferOpDamp = raw_input("Transfer operator damping: ")
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
      print bcolors.WARNING+"GenericRFactory cannot be used with non-smoothed PA-AMG prolongation operators. We change it back to TransPFactory."
      print "To use GenericRFactory you have to select either SaPFactory or PgPFactory for prolongation."+bcolors.ENDC
      print ""
      print "Press any key to proceed"
      waitForKey()
    self.isDirty = True
  
  # Rebalancing
  def doRebalancingOption(self):
    self.doRebalancing = True
    self.minRowsPerProc = raw_input("Minimum number of DOFs per processor: ")
    self.nnzImbalance = raw_input("Max. nonzero imbalance (default 1.1): ")
    self.rebStartLevel = raw_input("Start rebalancing on level (default 1): ")
    self.isDirty = True
    
  def doNoRebalancingOption(self):
    self.doRebalancing = False
    self.isDirty = True
  
  def runMenu(self,options,callbacks):
    for i,option in enumerate(options):
      print('%s. %s' % (i, option)) # display all options
    choice = raw_input('your choice? ')
    if is_number(str(choice)) and int(choice) < len(options):
      callbacks[int(choice)]() # call correspondending function
  
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
    print "doEXIT"
    self.exitLoop = True
    
  def printMainMenu(self):
    clearWindow()    
    self.printSettings()
    print ""
    print ""
    
    #options = ['Set Output file name','Common Multigrid settings', 'Level smoother settings', 'Transfer operators', 'Restriction operators', 'Save XML file', 'Exit']
    #callbacks = [self.doFileName, self.doCommonMenu, self.doSmootherMenu, self.doTransferMenu, self.doRestrictorMenu, self.generateXMLfile, self.doExitProgram]
    options = ['Common Multigrid settings', 'Aggregate settings', 'Level smoother settings', 'Transfer operators', 'Restriction operators', 'Rebalancing options', 'Save XML file', 'Back']
    callbacks = [self.doCommonMenu, self.doAggregatesMenu, self.doSmootherMenu, self.doTransferMenu, self.doRestrictorMenu, self.doRebalancingMenu, self.generateXMLfile, self.doExitProgram]
    
    self.runMenu(options,callbacks)      
        
  def printSettings(self):
    ## print out all made settings for xml file
    print bcolors.HEADER+"***************************   SETTINGS   ****************************"+bcolors.ENDC
    print bcolors.WARNING+"XML file name:           "+bcolors.ENDC + str(self.xmlFileName)
    print ""
    print bcolors.WARNING+"Max. MultiGrid levels:   "+bcolors.ENDC + str(self.maxMultLevels)
    print bcolors.WARNING+"Max. CoarseSize:         "+bcolors.ENDC + str(self.maxCoarseSize)
    print ""
    print bcolors.WARNING+"Level smoother:          "+bcolors.ENDC + str(self.levelSmoother)
    print bcolors.WARNING+"Level smoothing sweeps:  "+bcolors.ENDC + str(self.levelSmootherSweeps)
    print bcolors.WARNING+"Level damping parameter: "+bcolors.ENDC + str(self.levelSmootherDamp)
    print ""
    print bcolors.WARNING+"Coarse solver:           "+bcolors.ENDC + str(self.coarseSolver)
    print ""
    print bcolors.WARNING+"Graph drop tolerance:    "+bcolors.ENDC + str(self.dropTolerance)
    print bcolors.WARNING+"Aggregate size (min/max):"+bcolors.ENDC + str(self.minAggSize) + "/" + str(self.maxAggSize)
    print bcolors.WARNING+"Max. neighbor count:     "+bcolors.ENDC + str(self.maxNeighCount)
    print ""
    print bcolors.WARNING+"Transfer operators:     "+bcolors.ENDC + str(self.transferOps)
    print bcolors.WARNING+"Transfer smoothing par.:"+bcolors.ENDC + str(self.transferOpDamp)
    print ""
    print bcolors.WARNING+"Restriction operator:   "+bcolors.ENDC + str(self.restrictionOp)
    print ""
    if self.doRebalancing == False:
      print bcolors.WARNING+"NO Rebalancing"+bcolors.ENDC
    else:
      print bcolors.WARNING+"Rebalancing active:"+ bcolors.ENDC
      print bcolors.WARNING+"Minimum DOFs per proc:  "+ bcolors.ENDC + str(self.minRowsPerProc)
      print bcolors.WARNING+"Nonzero imbalance:      "+ bcolors.ENDC + str(self.nnzImbalance)      
      print bcolors.WARNING+"Start level for rebal.: "+ bcolors.ENDC + str(self.rebStartLevel)
    print bcolors.HEADER+"***************************   SETTINGS   ****************************"+bcolors.ENDC
    
    print ""
    if self.isDirty == True:
      print bcolors.FAIL+ "CHANGES HAVE NOT BEEN SAVED!" + bcolors.ENDC
    else:
      print bcolors.OKDARKGREEN + "CHANGES HAVE BEEN SAVED!" + bcolors.ENDC
      
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