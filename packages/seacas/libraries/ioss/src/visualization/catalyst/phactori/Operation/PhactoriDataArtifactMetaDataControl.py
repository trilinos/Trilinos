# Copyright(C) 1999-2020 National Technology & Engineering Solutions
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
import datetime

#phactori_combine_to_single_python_file_subpiece_begin_1

class PhactoriDataArtifactMetaDataControl:
  def __init__(self, myProcessId, numProcesses):
    self.outputEnabled = True
    self.artifactListCsvFileName = None
    self.artifactListCsvFileHandle = None
    self.artifactListTryOpenCount = 0
    self.thisProcessIsWritingProcess = -1
    if myProcessId == int(numProcesses/2):
      self.thisProcessIsWritingProcess = 1
    else:
      self.thisProcessIsWritingProcess = 0
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriDataArtifactMetaDataControl:__init__ " + \
        str(myProcessId) + ", " + str(numProcesses) + ", " + \
        str(self.thisProcessIsWritingProcess) + "\n", 100)

  def EnableOutput(self, onOff):
    self.outputEnabled = onOff

  def CloseFiles(self):
    if self.artifactListCsvFileHandle != None:
      self.artifactListCsvFileHandle.close()

  def ThisProcessIsWritingProcess(self):
    return self.thisProcessIsWritingProcess > 0

  def DataArtifactOutputListOpen(self):
    if self.artifactListCsvFileHandle == None:
      if self.artifactListTryOpenCount > 0:
        return False
      self.artifactListTryOpenCount += 1
      myNow = datetime.datetime.now()
      myFormat = "%Y-%m-%d-%H:%M:%S"
      self.artifactListCsvFileName = "DataArtifactList_" + myNow.strftime(myFormat) + ".csv"
      self.artifactListCsvFileHandle = open(self.artifactListCsvFileName, "w")
      if self.artifactListCsvFileHandle == None:
        print("AddImageToDataArtifactOutputList: unable to open file:\n",
          self.artifactListCsvFileName)
        return False
    return True

  def CheckIfDoOutput(self):
    if self.outputEnabled == False:
      return False

    if self.ThisProcessIsWritingProcess() == False:
      return False

    if self.DataArtifactOutputListOpen() == False:
      return False

    return True

  def AddDataArtifactToDataArtifactOutputList(self, fileName):

    if self.CheckIfDoOutput() == False:
      return

    myNow = datetime.datetime.now()
    outDate = str(myNow)
    outStr = fileName + "," + outDate + "\n"
    self.artifactListCsvFileHandle.write(outStr)
    self.artifactListCsvFileHandle.flush()

  def AddDataExportToDataArtifactOutputList(self, fileName):
    self.AddDataArtifactToDataArtifactOutputList(fileName)

  def AddImageToDataArtifactOutputList(self, fileName):
    self.AddDataArtifactToDataArtifactOutputList(fileName)

if __name__ == '__main__':
  pdamdc = PhactoriDataArtifactMetaDataControl(1,2)
  pdamdc.AddImageToDataArtifactOutputList("CatalystOutput/test1.png")
  pdamdc.AddDataExportToDataArtifactOutputList("CatalystVtkDataOutput/test2.vtm")

#phactori_combine_to_single_python_file_subpiece_end_1

    

      

