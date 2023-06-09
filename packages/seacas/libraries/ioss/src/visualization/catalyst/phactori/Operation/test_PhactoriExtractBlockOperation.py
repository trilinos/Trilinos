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
import unittest
#from vtk import *
import vtk
from Operation.PhactoriExtractBlockOperation import *
import paraview.simple

class TestPhactoriExtractBlockOperation(unittest.TestCase):

  def test_ExtractOneBlockFromAGroupWithThreeBlocks(self):
    cone1 = Cone(registrationName = "cone1")
    sphere1 = Sphere(registrationName = "sphere1")
    cylinder1 = Sphere(registrationName = "cylinder1")
    group1 = GroupDatasets(registrationName = "group1", Input = [cone1, sphere1, cylinder1])
    group1.UpdatePipeline()

    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "testoperation"
    operationParams = {"include blocks": ["/Root/cone1"]}
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'extractblock',
              PhactoriExtractBlockOperation,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, group1)
    newOperationBlock.GetPvFilter().UpdatePipeline()
    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
    if PhactoriDbg(100):
      myDebugPrint3("group number of blocks:" + str(csdata.GetNumberOfBlocks()))
    self.assertEqual(csdata.GetNumberOfBlocks(), 1)

  def test_ExtractTwoBlocksFromAGroupWithThgeeBlocks(self):
    cone1 = Cone(registrationName = "cone1")
    sphere1 = Sphere(registrationName = "sphere1")
    cylinder1 = Sphere(registrationName = "cylinder1")
    group1 = GroupDatasets(registrationName = "group1", Input = [cone1, sphere1, cylinder1])
    group1.UpdatePipeline()

    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "testoperation"
    operationParams = {"include blocks": ["/Root/cone1", "/Root/cylinder1"]}
    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'extractblock',
              PhactoriExtractBlockOperation,
              operationParams)
    ConstructPipelineOperationFromParsedOperationBlockC_ForTest(newOperationBlock, group1)
    newOperationBlock.GetPvFilter().UpdatePipeline()
    csdata = newOperationBlock.GetPvFilter().GetClientSideObject().GetOutputDataObject(0)
    if PhactoriDbg(100):
      myDebugPrint3("group number of blocks:" + str(csdata.GetNumberOfBlocks()))
    self.assertEqual(csdata.GetNumberOfBlocks(), 2)

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()


