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

#phactori_combine_to_single_python_file_subpiece_begin_1
class PhactoriPartitionedToMultiBlockOperation(PhactoriOperationSpecifics):
  """manages PartitionedToMultiBlock filter"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.InputFilter = None
    self.OutVtkMultiBlockDataSet = None
    self.OutputPVTrivialProducer = None

  def ParseParametersFromJson(self, inJson):
    dummy = 0

  def GetNodeType(self, inNode):
    if inNode == None:
      myDebugPrint3AndException("PhactoriVtkDataExportOperation::GetNodeType\n"
          "fix to handle None inNode\n")
    icsdClassname = inNode.GetClassName()
    if icsdClassname == "vtkMultiBlockDataSet":
      nodeType = 0
    elif icsdClassname == "vtkExodusIIMultiBlockDataSet":
      nodeType = 1
    elif icsdClassname == "vtkMultiPieceDataSet":
      nodeType = 2
    elif icsdClassname == "vtkPartitionedDataSet":
      nodeType = 3
    elif icsdClassname == "vtkStructuredGrid":
      nodeType = 4
    elif icsdClassname == "vtkUntructuredGrid":
      nodeType = 5
    else:
      nodeType = 6
    return nodeType

  def AddItemToMultiblockOrMultiPiece(self, nonLeafType, inCurrentAlteredNode, oneChild, childName):
    if nonLeafType == 2:
      newChildIndex = inCurrentAlteredNode.GetNumberOfPieces()
      inCurrentAlteredNode.SetPiece(newChildIndex, oneChild)
    elif nonLeafType == 3:
      #newChildIndex = inCurrentAlteredNode.GetNumberOfPartitions()
      #inCurrentAlteredNode.SetPartition(newChildIndex, oneChild)
      myDebugPrint3AndException(
        "PhactoriVtkDataExportOperation::ParseParametersFromJson\n"
        "unexpected nonLeafType (vtkPartitionedDataSet) at this place\n")
    else:
      newChildIndex = inCurrentAlteredNode.GetNumberOfBlocks()
      inCurrentAlteredNode.SetBlock(newChildIndex, oneChild)
    if childName != None:
      newAlteredNodeMetaData = inCurrentAlteredNode.GetMetaData(newChildIndex)
      newAlteredNodeMetaData.Set(vtk.vtkCompositeDataSet.NAME(), childName)

  def UpdateAlteredVtkMultiBlockDataSetRecurse1(self, inCurrentAlteredNode, inCurrentSourceNode):
    #at this point, both inputs are not leafs, and inCurrentSource is not
    #VtkPartitionedDataSet
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriVtkDataExportOperation." \
          "UpdateAlteredVtkMultiBlockDataSetRecurse1 entered\n"
          "inCurrentAlteredNode type: " + str(inCurrentAlteredNode.GetClassName()) + "\n"
          "inCurrentSourceNode type:  " + str(inCurrentSourceNode.GetClassName()) + "\n", 100)

    nonLeafType = self.GetNodeType(inCurrentSourceNode)
    if PhactoriDbg(100):
      myDebugPrint3("inCurrentAlteredNode nonLeafType " + str(nonLeafType) + "\n")

    if nonLeafType > 3:
      myDebugPrint3AndException(
        "PhactoriVtkDataExportOperation::UpdateAlteredVtkMultiBlockDataSetRecurse1\n"
        "this is a leaf node or unexpected class here.\n"
        "GetClassName(): " + str(inCurrentSourceNode.GetClassName()) + "\n")

    if nonLeafType == 2:
      numChildren = inCurrentSourceNode.GetNumberOfPieces()
    elif nonLeafType == 3:
      numChildren = inCurrentSourceNode.GetNumberOfPartitions()
    else:
      numChildren = inCurrentSourceNode.GetNumberOfBlocks()

    if PhactoriDbg(100):
      myDebugPrint3("nonLeafType: " + str(nonLeafType) + " numChildren: " + str(numChildren) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("begin quick children check\n")
      for childIndex in range(0, numChildren):
        if nonLeafType == 2:
          oneChild = inCurrentSourceNode.GetPiece(childIndex)
        elif nonLeafType == 3:
          oneChild = inCurrentSourceNode.GetPartition(childIndex)
        else:
          oneChild = inCurrentSourceNode.GetBlock(childIndex)
        if oneChild != None:
          myDebugPrint3("childIndex " + str(childIndex) + " is " + str(oneChild.GetClassName()) + "\n")
        else:
          myDebugPrint3("childIndex " + str(childIndex) + " is None\n")
      myDebugPrint3("end quick children check\n")

    for childIndex in range(0, numChildren):
      if PhactoriDbg(100):
        myDebugPrint3("doing childIndex: " + str(childIndex) + "\n")

      if nonLeafType == 2:
        oneChild = inCurrentSourceNode.GetPiece(childIndex)
      elif nonLeafType == 3:
        oneChild = inCurrentSourceNode.GetPartition(childIndex)
      else:
        oneChild = inCurrentSourceNode.GetBlock(childIndex)

      childMetaData = inCurrentSourceNode.GetMetaData(childIndex)
      if childMetaData != None:
        childName = childMetaData.Get(vtk.vtkCompositeDataSet.NAME())
      else:
        childName = None

      if oneChild == None:
        #for now, just add a None block to match
        self.AddItemToMultiblockOrMultiPiece(nonLeafType, inCurrentAlteredNode, oneChild, childName)
        if PhactoriDbg(100):
          myDebugPrint3("done childIndex: " + str(childIndex) + " (oneChild was None)\n")
        continue

      childNodeType = self.GetNodeType(oneChild)

      if PhactoriDbg(100):
        myDebugPrint3("childIndex " + str(childIndex) + " childNodeType " + str(childNodeType) + " childName " + str(childName) + "\n")
      if (childNodeType == 4) or (childNodeType == 5):
        self.AddItemToMultiblockOrMultiPiece(nonLeafType, inCurrentAlteredNode, oneChild, childName)
      elif (childNodeType == 0) or (childNodeType == 1):
        #newAlteredNode = vtk.VtkMultiBlockDataSet()
        #newAlteredNode = oneChild.NewInstance()
        #we are just using self.OutVtkMultiBlockDataSet here to get
        #access to the desired type of NewInstance regardless if
        #childNodeType is 0 or 1
        newAlteredNode = self.OutVtkMultiBlockDataSet.NewInstance()
        newChildIndex = inCurrentAlteredNode.GetNumberOfBlocks()
        self.AddItemToMultiblockOrMultiPiece(nonLeafType, inCurrentAlteredNode, newAlteredNode, childName)
        self.UpdateAlteredVtkMultiBlockDataSetRecurse1(newAlteredNode, oneChild)
      elif childNodeType == 2:
        #newAlteredNode = vtk.VtkMultiPieceDataSet()
        newAlteredNode = oneChild.NewInstance()
        newChildIndex = inCurrentAlteredNode.GetNumberOfPieces()
        self.AddItemToMultiblockOrMultiPiece(nonLeafType, inCurrentAlteredNode, newAlteredNode, childName)
        self.UpdateAlteredVtkMultiBlockDataSetRecurse1(newAlteredNode, oneChild)
      elif childNodeType == 3:
        #for now we're making assumption all the partitions under here are leaf grids
        numSubChildren = oneChild.GetNumberOfPartitions()
        if PhactoriDbg(100):
          myDebugPrint3("numSubChildren: " + str(numSubChildren) + "\n")
        for subChildIndex in range(0, numSubChildren):
          if PhactoriDbg(100):
            myDebugPrint3("doing subChildIndex: " + str(subChildIndex) + "\n")
          oneSubChild = oneChild.GetPartition(subChildIndex)

          if oneSubChild == None:
            self.AddItemToMultiblockOrMultiPiece(nonLeafType, inCurrentAlteredNode, oneSubChild, childName)
            if PhactoriDbg(100):
              myDebugPrint3("done subChildIndex: " + str(subChildIndex) + " (which was a None subchild)\n")
            continue

          subChildNodeType = self.GetNodeType(oneSubChild)
          if PhactoriDbg(100):
            myDebugPrint3("subChildIndex " + str(subChildIndex) + " subChildNodeType " + str(subChildNodeType) + "\n")
          if (subChildNodeType == 4) or (subChildNodeType == 5):
            if (nonLeafType == 0) or (nonLeafType == 1) or (nonLeafType == 2):
              self.AddItemToMultiblockOrMultiPiece(nonLeafType, inCurrentAlteredNode, oneSubChild, childName)
            else:
              myDebugPrint3AndException(
                "PhactoriVtkDataExportOperation::ParseParametersFromJson\n"
                "unexpected nonLeafType at this place: " + str(nonLeafType) + "\n")
          else:
            myDebugPrint3AndException(
              "PhactoriVtkDataExportOperation::ParseParametersFromJson\n"
              "unexpected subChildNodeType at this place: " + str(subChildNodeType) + "\n")
          if PhactoriDbg(100):
            myDebugPrint3("done subChildIndex: " + str(subChildIndex) + "\n")
      else:
        myDebugPrint3AndException(
          "PhactoriVtkDataExportOperation::ParseParametersFromJson\n"
          "unexpected childNodeType\n")
      if PhactoriDbg(100):
        myDebugPrint3("done childIndex: " + str(childIndex) + "\n")

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriVtkDataExportOperation.ExportOperationData returning\n")


  def MakeOrUpdateVtkMultiBlockDataSetWithtoutVtkPartitionedDataSets(self, inVtkMultiBlockDataSet):
    if self.OutVtkMultiBlockDataSet == None:
      self.OutVtkMultiBlockDataSet = inVtkMultiBlockDataSet.NewInstance()
    else:
      self.OutVtkMultiBlockDataSet.SetNumberOfBlocks(0)
    self.UpdateAlteredVtkMultiBlockDataSetRecurse1(self.OutVtkMultiBlockDataSet, inVtkMultiBlockDataSet)

  def DoUpdateDueToChangeInData(self, inIncomingPvFilter,
      outOutgoingPvFilter):
    #default behavior is to do nothing on new data
    if PhactoriDbg():
      myDebugPrint3("PhactoriPartitionedToMultiBlockOperation.DoUpdateDueToChangeInData entered\n")
    if self.InputFilter != None:
      if PhactoriDbg():
        myDebugPrint3("self.InputFilter != None; updating\n")
      UpdatePipelineWithCurrentTimeArgument(self.InputFilter)
      theVtkMultiBlockDataSetObject = self.InputFilter.GetClientSideObject().GetOutputDataObject(0)
      self.MakeOrUpdateVtkMultiBlockDataSetWithtoutVtkPartitionedDataSets(theVtkMultiBlockDataSetObject)
      UpdatePipelineWithCurrentTimeArgument(self.OutputPVTrivialProducer)
    if PhactoriDbg():
      myDebugPrint3("PhactoriPartitionedToMultiBlockOperation.DoUpdateDueToChangeInData returning\n")
    return True

  def CreateParaViewFilter(self, inInputFilter):
    """create the PartitionedToMultiBlock filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPartitionedToMultiBlockOperation.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    self.InputFilter = inInputFilter
    theVtkMultiBlockDataSetObject = self.InputFilter.GetClientSideObject().GetOutputDataObject(0)

    self.MakeOrUpdateVtkMultiBlockDataSetWithtoutVtkPartitionedDataSets(theVtkMultiBlockDataSetObject)

    #make trivial producer with this block
    newParaViewFilter = PVTrivialProducer()
    newParaViewFilter.GetClientSideObject().SetOutput(self.OutVtkMultiBlockDataSet)

    self.OutputPVTrivialProducer = newParaViewFilter

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriPartitionedToMultiBlockOperation.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter

#phactori_combine_to_single_python_file_subpiece_end_1
