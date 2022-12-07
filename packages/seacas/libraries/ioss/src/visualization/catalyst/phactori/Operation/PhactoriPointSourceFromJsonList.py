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
from paraview.simple import *
from .PhactoriMpiUtilities import *

#phactori_combine_to_single_python_file_subpiece_begin_1

class PhactoriPointSourceFromJsonList(PhactoriOperationSpecifics):
  """Filter/operation which reads in a .json file which is a list of 3d
     coordinates and creates a new point source from that list. Resulting
     point source will have 1 element an N points. This source is intented
     to work correctly in parallel for Catalyst or pvbatch symmetric mode.
     The json file which is read in will only be read on one process and
     mpi broadcast is used to distribute the list, rather than having each
     process read the json file."""
  def __init__(self):
    self.JsonListFileName = "PhactoriPointSourceFromJsonList.json"
    self.JsonList = None
    self.ParaviewPointSource = None
    self.SourceVtkPolyData = None
    self.myVtkPoints = vtk.vtkPoints()
    self.myVtkPolyVertex = vtk.vtkPolyVertex()
    self.myVtkCellArray = vtk.vtkCellArray()
    self.myVtkPolyData = vtk.vtkPolyData()

  def ValidateJsonPointList(self):
    numPoints = len(self.JsonList)
    if numPoints < 1:
      myDebugPrint3AndException(
          "PhactoriPointSourceFromJsonList::ValidateJsonPointList\n"
          "list must have at least one element\n")

    for ptNdx in range(0,numPoints):
      jsonPt = self.JsonList[ptNdx]
      if len(jsonPt) != 3:
        errStr = "PhactoriPointSourceFromJsonList::ValidateJsonPointList\n" \
          "point with index " + str(ptNdx) + "does not have three elements\n"
        myDebugPrint3AndException(errStr)

    return True

  def ParseParametersFromJson(self, inJson):
    if 'filename' in inJson:
      self.JsonListFileName = inJson['filename']
    else:
      myDebugPrint3AndException(
          "PhactoriPointSourceFromJsonList::ParseParametersFromJson\n"
          "Error:  must have 'filename' key\n")

  def CreateVtkPolyDataFromJsonList(self):
    numPoints = len(self.JsonList)
    self.myVtkPoints.SetNumberOfPoints(numPoints)
    self.myVtkPolyVertex.GetPointIds().SetNumberOfIds(numPoints)

    vtkPolyVertPtIds = self.myVtkPolyVertex.GetPointIds()
    usePoint = [0.0, 0.0, 0.0]
    for ptNdx in range(0,numPoints):
      #self.myVtkPoints.SetPoint(ptNdx, self.JsonList[ptNdx])
      jsonPt = self.JsonList[ptNdx]
      usePoint[0] = jsonPt[0]
      usePoint[1] = jsonPt[1]
      usePoint[2] = jsonPt[2]
      self.myVtkPoints.SetPoint(ptNdx, usePoint)
      vtkPolyVertPtIds.SetId(ptNdx, ptNdx)

    self.myVtkPolyData.SetPoints(self.myVtkPoints)
    self.myVtkCellArray.InsertNextCell(self.myVtkPolyVertex)
    self.myVtkPolyData.SetVerts(self.myVtkCellArray)

  def CreateParaViewFilter(self, inInputFilter):
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPointSourceFromJsonList.CreateParaViewFilter entered\n", 100)

    savedActiveSource = GetActiveSource()

    self.JsonList = ReadAndMpiBroadcastJsonFile(self.JsonListFileName)

    self.ValidateJsonPointList()

    self.CreateVtkPolyDataFromJsonList()

    self.ParaviewPointSource = PVTrivialProducer()
    self.ParaviewPointSource.GetClientSideObject().SetOutput(self.myVtkPolyData)

    SetActiveSource(self.ParaviewPointSource)
    SetActiveSource(savedActiveSource)

    UpdatePipelineWithCurrentTimeArgument(self.ParaviewPointSource)

    if PhactoriDbg(100):
      pssfjlcso = self.ParaviewPointSource.GetClientSideObject()
      pssfjlodo = pssfjlcso.GetOutputDataObject(0)
      myDebugPrint3("pssfjlodo: " + str(pssfjlodo.GetClassName()) + " numcells " + str(pssfjlodo.GetNumberOfCells()) + " numpoints " + str(pssfjlodo.GetNumberOfPoints()) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPointSourceFromJsonList.CreateParaViewFilter returning\n", 100)
    return self.ParaviewPointSource

#phactori_combine_to_single_python_file_subpiece_end_1
