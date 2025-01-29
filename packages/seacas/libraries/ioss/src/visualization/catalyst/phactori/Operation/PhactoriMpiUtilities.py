# Copyright(C) 1999-2020, 2024 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     * Neither the name of NTESS nor the names of its
#       contributors may be used to endorse or promote products derived
#       from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from phactori import *
from paraview.simple import *

#phactori_combine_to_single_python_file_subpiece_begin_1

global BarrierLockCount
BarrierLockCount = 0
def BarrierLock(tagString = ""):
  """called to sync up all processes"""
  global BarrierLockCount
  BarrierLockCount += 1
  numProcs = SmartGetNumberOfProcesses()
  myPid = SmartGetLocalProcessId()
  localIntList = [myPid]
  if PhactoriDbg(100):
    myDebugPrint3("BarrierLock A " + str(BarrierLockCount) + " " + tagString + "\n")
  globalPidSumList = UseReduceOnIntegerList(localIntList, 2)
  if PhactoriDbg(100):
    myDebugPrint3("BarrierLock B " + str(BarrierLockCount) + " " + tagString + "\n")

def UseReduceOnFloatList(inFloatList, inAllReduceType):
  """inAllReduceType 0-max 1-min 2-sum"""
  #if PhactoriDbg(100):
  #  myDebugPrint3("UseReduceOnFloatList entered " + str(inAllReduceType) + "\n", 100)

  #if PhactoriDbg(100):
  #  myDebugPrint3("num values: " + str(len(inFloatList)) + "\n")
  #  if len(inFloatList) > 0:
  #    myDebugPrint3("local first: " + str(inFloatList[0]) + "   last: " + str(inFloatList[-1]) + "\n")

  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()

  numVals = len(inFloatList)
  localarray = vtk.vtkDoubleArray()
  localarray.SetNumberOfTuples(numVals)
  globalarray = vtk.vtkDoubleArray()
  globalarray.SetNumberOfTuples(numVals)

  for idx, oneVal in enumerate(inFloatList):
    localarray.SetValue(idx, oneVal)

  globalController.AllReduce(localarray, globalarray, inAllReduceType)

  returnArray = []
  for idx in range(0, numVals):
    returnArray.append(globalarray.GetTuple1(idx))

  #if PhactoriDbg(100):
  #  if len(returnArray) > 0:
  #    myDebugPrint3("global first: " + str(returnArray[0]) + "   last: " + str(returnArray[-1]) + "\n")

  #if PhactoriDbg(100):
  #  myDebugPrint3("UseReduceOnFloatList returning\n", 100)
  return returnArray

def UseReduceOnIntegerList(inIntegerList, inAllReduceType):
  #if PhactoriDbg(100):
  #  myDebugPrint3("UseReduceOnIntegerList entered\n", 100)

  #if PhactoriDbg(100):
  #  myDebugPrint3("num values: " + str(len(inIntegerList)) + "\n")
  #  if len(inIntegerList) > 0:
  #    myDebugPrint3("local first: " + str(inIntegerList[0]) + "   last: " + str(inIntegerList[-1]) + "\n")

  pm = paraview.servermanager.vtkProcessModule.GetProcessModule()
  globalController = pm.GetGlobalController()

  numVals = len(inIntegerList)
  localarray = vtk.vtkIntArray()
  localarray.SetNumberOfTuples(numVals)
  globalarray = vtk.vtkIntArray()
  globalarray.SetNumberOfTuples(numVals)

  for idx, oneVal in enumerate(inIntegerList):
    localarray.SetValue(idx, oneVal)

  globalController.AllReduce(localarray, globalarray, inAllReduceType)

  returnArray = []
  for idx in range(0, numVals):
    returnArray.append(int(globalarray.GetTuple1(idx)))

  #if PhactoriDbg(100):
  #  if len(returnArray) > 0:
  #    myDebugPrint3("global first: " + str(returnArray[0]) + "   last: " + str(returnArray[-1]) + "\n")

  #if PhactoriDbg(100):
  #  myDebugPrint3("UseReduceOnIntegerList returning\n", 100)
  return returnArray

def UseMpiToShareAllProcessArraySizes(inArray):
  if PhactoriDbg(100):
    myDebugPrint3("UseMpiToShareAllProcessArraySizes entered\n", 100)
  numProcs = SmartGetNumberOfProcesses()
  myPid = SmartGetLocalProcessId()

  localArraySize = []
  for ii in range(0, numProcs):
    localArraySize.append(int(0))

  localArraySize[myPid] = int(len(inArray))
  globalArraySize = UseReduceOnIntegerList(localArraySize, 0)
  if PhactoriDbg(100):
    myDebugPrint3(str(localArraySize) + "\n")
    myDebugPrint3(str(globalArraySize) + "\n")
    myDebugPrint3("UseMpiToShareAllProcessArraySizes returning\n", 100)
  return globalArraySize

def UseMpiToSendAllProcessesFloatArrayToOneProcess(thisProcessFloatArray, pidToCollectArray):
  """given a float array on this process, send the float arrays from all
     processes to the process pidToCollectArray.  Gets sizes from all
     processes first and then uses reduce to send it all"""

  if PhactoriDbg(100):
    myDebugPrint3("UseMpiToSendAllProcessesFloatArrayToOneProcess entered\n", 100)

  numProcs = SmartGetNumberOfProcesses()
  myPid = SmartGetLocalProcessId()

  globalArraySize = UseMpiToShareAllProcessArraySizes(thisProcessFloatArray)

  localFloatArray = []
  for pidIndex, oneSize in enumerate(globalArraySize):
    if pidIndex == myPid:
      for vv in thisProcessFloatArray:
        localFloatArray.append(vv)
    else:
      for cnt1 in range(0, oneSize):
        localFloatArray.append(0.0)

  #we really don't need to use AllReduce here, just Reduce to one process
  #fix later
  globalFloatArray = UseReduceOnFloatList(localFloatArray, 2)

  #we really should only have value for pidToCollectArray
  if myPid == pidToCollectArray:
    retVal = globalFloatArray
  else:
    retVal = []

  if PhactoriDbg(100):
    if len(globalFloatArray) == 0:
      myDebugPrint3(str(len(globalFloatArray)) + "\n")
    else:
      myDebugPrint3(str(len(globalFloatArray)) + "  " + str(globalFloatArray[0]) + "  " + str(globalFloatArray[-1]) + "\n")
    myDebugPrint3("UseMpiToSendAllProcessesFloatArrayToOneProcess returning\n", 100)
  return retVal

def UseMpiToSendAllProcessesIntArrayToOneProcess(thisProcessIntArray, pidToCollectArray):
  """given a int array on this process, send the int arrays from all
     processes to the process pidToCollectArray.  Gets sizes from all
     processes first and then uses reduce to send it all"""
  if PhactoriDbg(100):
    myDebugPrint3("UseMpiToSendAllProcessesintArrayToOneProcess entered\n", 100)
  numProcs = SmartGetNumberOfProcesses()
  myPid = SmartGetLocalProcessId()

  globalArraySize = UseMpiToShareAllProcessArraySizes(thisProcessIntArray)

  localIntArray = []
  for pidIndex, oneSize in enumerate(globalArraySize):
    if pidIndex == myPid:
      for vv in thisProcessIntArray:
        localIntArray.append(vv)
    else:
      for cnt1 in range(0, oneSize):
        localIntArray.append(int(0))

  #we really don't need to use AllReduce here, just Reduce to one process
  #fix later
  globalIntArray = UseReduceOnIntegerList(localIntArray, 2)

  #we really should only have value for pidToCollectArray
  if myPid == pidToCollectArray:
    retVal = globalIntArray
  else:
    retVal = []

  if PhactoriDbg(100):
    if len(globalIntArray) == 0:
      myDebugPrint3(str(len(globalIntArray)) + "\n")
    else:
      myDebugPrint3(str(len(globalIntArray)) + "  " + str(globalIntArray[0]) + "  " + str(globalIntArray[-1]) + "\n")
    myDebugPrint3("UseMpiToSendAllProcessesintArrayToOneProcess returning\n", 100)
  return retVal

def ReadAndMpiBroadcastJsonFile(inJsonFileName):
  """given a json format file, read in the contents using the python json
     library on process zero, use MPI to broadcast the contents to all
     processes (via converting to a long char string)"""
  if SmartGetLocalProcessId() == 0:
    #process zero loads json script and broadcasts
    returnJson = HandleJsonScriptLoadProcessZero(inJsonFileName)
  else:
    #other processes receive json script
    returnJson = HandleJsonScriptLoadProcessNotZero()
  return returnJson

#phactori_combine_to_single_python_file_subpiece_end_1

