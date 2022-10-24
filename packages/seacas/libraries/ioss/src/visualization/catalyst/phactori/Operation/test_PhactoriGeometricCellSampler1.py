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
from Operation.PhactoriSampledCellInfo import *
from Operation.PhactoriGeometricCellSampler1 import *
from paraview.simple import *
import os

class TestPhactoriGeometricCellSampler1(unittest.TestCase):

  def test_PointIsInsideBoundingBox(self):
    pgcs1_inst = PhactoriGeometricCellSampler1()
    testbb = [-1.25,1.25,-1.5,2.0,1.75,4.25]
    self.assertFalse(pgcs1_inst.PointIsInsideBoundingBox([-1.26, 0.0, 2.0], testbb))
    self.assertFalse(pgcs1_inst.PointIsInsideBoundingBox([1.26, 0.0, 2.0], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([-1.25, 0.0, 2.0], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([1.25, 0.0, 2.0], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([-1.24, 0.0, 2.0], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([1.24, 0.0, 2.0], testbb))

    self.assertFalse(pgcs1_inst.PointIsInsideBoundingBox([0.0, -1.6, 2.0], testbb))
    self.assertFalse(pgcs1_inst.PointIsInsideBoundingBox([0.0, 2.1, 2.0], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([0.0, -1.5, 2.0], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([0.0, 1.9, 2.0], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([0.0, -1.4, 2.0], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([0.0, 1.9, 2.0], testbb))

    self.assertFalse(pgcs1_inst.PointIsInsideBoundingBox([0.0, 0.0, 1.74], testbb))
    self.assertFalse(pgcs1_inst.PointIsInsideBoundingBox([0.0, 0.0, 4.26], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([0.0, 0.0, 1.75], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([0.0, 0.0, 4.25], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([0.0, 0.0, 1.76], testbb))
    self.assertTrue(pgcs1_inst.PointIsInsideBoundingBox([0.0, 0.0, 4.24], testbb))

  def test_ParseParametersFromJson(self):
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "phactorigeometriccellsampler1"
    badOperationParams = {
    "type":"geometriccellsampler1",
    "cell data array names":["RTData"],
    "cell data array tuple size": 1,
    "do programmable filter": False,
    "data controlled sampling method":"cells with ratio of min/max data value",
    "data controlled ratio of min/max":0.75,
    "sampling geometry bounding box":[-1.5, -4.5, -2.5, 3.5, 3.25, 7.25]
    }
    with self.assertRaises(Exception):
      ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'geometriccellsampler1',
              PhactoriGeometricCellSampler1,
              badOperationParams)

  def test_CreateInternalListOfGeometricallySampledCellsOnThisProcess(self):
    testWavelet2 = Wavelet()
    testWavelet2.UpdatePipeline()
    testWavelet = PointDatatoCellData(Input=testWavelet2)
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "phactorigeometriccellsampler1"

    testOutFileBasename = "test_WriteAllDataFromOneProcessUsingMPI_output_cells_"

    operationParams = {
    "type":"geometriccellsampler1",
    "cell data array names":["RTData"],
    "cell data array tuple size": 1,
    "do programmable filter": False,
    "data controlled sampling method":"cells with ratio of min/max data value",
    "data controlled ratio of min/max":0.75,
    "sampling geometry bounding box":[-4.5, -1.5, -2.5, 3.5, 3.25, 7.25]
    }

    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'geometriccellsampler1',
              PhactoriGeometricCellSampler1,
              operationParams)

    PhactoriGeometricCellSampler1Instance = newOperationBlock.mOperationSpecifics
    PhactoriGeometricCellSampler1Instance.myCopyOfInputFilter = testWavelet
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfGeometricallySampledCellsOnThisProcess()
    PhactoriGeometricCellSampler1Instance.WriteCellListToFile("test_PhactoriGeometricCellSampler1_1.json",
      PhactoriGeometricCellSampler1Instance.GeometricallySampledCellsForThisProcess)

    goldTestStringJson = """
    {
    "format for cells in lists": {"PhactoriSampledCellInfo output format 1 info":[
    " [cellTestPoint, ijk, dataTuple, pid, index, segmentIndex, collectionAxis]",
    " cellTestPoint is [X, Y, Z], ijk is [i, j, k], dataTuple is [c1, c2, ... cN]"]},
    "sampling geometry bounding box": [-4.5, -1.5, -2.5, 3.5, 3.25, 7.25],
    "number of cells": 112,
    "cell variable names": ["RTData"],
    "data tuple size": 1,
    "sum variable values": [23509.975616455078],
    "average variable values": [209.91049657549178],
    "list of sampled cells": [
    [[-4.5, -2.5, 3.5], [-1, -1, -1], [204.47789001464844], 0, 1, 5345, -1, -1],
    [[-3.5, -2.5, 3.5], [-1, -1, -1], [211.81105041503906], 0, 1, 5346, -1, -1],
    [[-2.5, -2.5, 3.5], [-1, -1, -1], [219.71676635742188], 0, 1, 5347, -1, -1],
    [[-1.5, -2.5, 3.5], [-1, -1, -1], [222.8864288330078], 0, 1, 5348, -1, -1],
    [[-4.5, -1.5, 3.5], [-1, -1, -1], [226.4726104736328], 0, 1, 5365, -1, -1],
    [[-3.5, -1.5, 3.5], [-1, -1, -1], [233.97752380371094], 0, 1, 5366, -1, -1],
    [[-2.5, -1.5, 3.5], [-1, -1, -1], [242.0166473388672], 0, 1, 5367, -1, -1],
    [[-1.5, -1.5, 3.5], [-1, -1, -1], [245.2774658203125], 0, 1, 5368, -1, -1],
    [[-4.5, -0.5, 3.5], [-1, -1, -1], [227.3440399169922], 0, 1, 5385, -1, -1],
    [[-3.5, -0.5, 3.5], [-1, -1, -1], [234.93612670898438], 0, 1, 5386, -1, -1],
    [[-2.5, -0.5, 3.5], [-1, -1, -1], [243.04296875], 0, 1, 5387, -1, -1],
    [[-1.5, -0.5, 3.5], [-1, -1, -1], [246.3500518798828], 0, 1, 5388, -1, -1],
    [[-4.5, 0.5, 3.5], [-1, -1, -1], [209.38912963867188], 0, 1, 5405, -1, -1],
    [[-3.5, 0.5, 3.5], [-1, -1, -1], [216.98121643066406], 0, 1, 5406, -1, -1],
    [[-2.5, 0.5, 3.5], [-1, -1, -1], [225.08804321289062], 0, 1, 5407, -1, -1],
    [[-1.5, 0.5, 3.5], [-1, -1, -1], [228.3951416015625], 0, 1, 5408, -1, -1],
    [[-4.5, 1.5, 3.5], [-1, -1, -1], [205.9775390625], 0, 1, 5425, -1, -1],
    [[-3.5, 1.5, 3.5], [-1, -1, -1], [213.48245239257812], 0, 1, 5426, -1, -1],
    [[-2.5, 1.5, 3.5], [-1, -1, -1], [221.52157592773438], 0, 1, 5427, -1, -1],
    [[-1.5, 1.5, 3.5], [-1, -1, -1], [224.7823944091797], 0, 1, 5428, -1, -1],
    [[-4.5, 2.5, 3.5], [-1, -1, -1], [219.5332794189453], 0, 1, 5445, -1, -1],
    [[-3.5, 2.5, 3.5], [-1, -1, -1], [226.86642456054688], 0, 1, 5446, -1, -1],
    [[-2.5, 2.5, 3.5], [-1, -1, -1], [234.77215576171875], 0, 1, 5447, -1, -1],
    [[-1.5, 2.5, 3.5], [-1, -1, -1], [237.94180297851562], 0, 1, 5448, -1, -1],
    [[-4.5, 3.5, 3.5], [-1, -1, -1], [217.14462280273438], 0, 1, 5465, -1, -1],
    [[-3.5, 3.5, 3.5], [-1, -1, -1], [224.22647094726562], 0, 1, 5466, -1, -1],
    [[-2.5, 3.5, 3.5], [-1, -1, -1], [231.93704223632812], 0, 1, 5467, -1, -1],
    [[-1.5, 3.5, 3.5], [-1, -1, -1], [234.97329711914062], 0, 1, 5468, -1, -1],
    [[-4.5, -2.5, 4.5], [-1, -1, -1], [191.78919982910156], 0, 1, 5745, -1, -1],
    [[-3.5, -2.5, 4.5], [-1, -1, -1], [198.7889404296875], 0, 1, 5746, -1, -1],
    [[-2.5, -2.5, 4.5], [-1, -1, -1], [206.43572998046875], 0, 1, 5747, -1, -1],
    [[-1.5, -2.5, 4.5], [-1, -1, -1], [209.42840576171875], 0, 1, 5748, -1, -1],
    [[-4.5, -1.5, 4.5], [-1, -1, -1], [213.61886596679688], 0, 1, 5765, -1, -1],
    [[-3.5, -1.5, 4.5], [-1, -1, -1], [220.7836456298828], 0, 1, 5766, -1, -1],
    [[-2.5, -1.5, 4.5], [-1, -1, -1], [228.55862426757812], 0, 1, 5767, -1, -1],
    [[-1.5, -1.5, 4.5], [-1, -1, -1], [231.63890075683594], 0, 1, 5768, -1, -1],
    [[-4.5, -0.5, 4.5], [-1, -1, -1], [214.4065399169922], 0, 1, 5785, -1, -1],
    [[-3.5, -0.5, 4.5], [-1, -1, -1], [221.65509033203125], 0, 1, 5786, -1, -1],
    [[-2.5, -0.5, 4.5], [-1, -1, -1], [229.4951171875], 0, 1, 5787, -1, -1],
    [[-1.5, -0.5, 4.5], [-1, -1, -1], [232.61985778808594], 0, 1, 5788, -1, -1],
    [[-4.5, 0.5, 4.5], [-1, -1, -1], [196.45162963867188], 0, 1, 5805, -1, -1],
    [[-3.5, 0.5, 4.5], [-1, -1, -1], [203.70018005371094], 0, 1, 5806, -1, -1],
    [[-2.5, 0.5, 4.5], [-1, -1, -1], [211.5402069091797], 0, 1, 5807, -1, -1],
    [[-1.5, 0.5, 4.5], [-1, -1, -1], [214.66494750976562], 0, 1, 5808, -1, -1],
    [[-4.5, 1.5, 4.5], [-1, -1, -1], [193.12379455566406], 0, 1, 5825, -1, -1],
    [[-3.5, 1.5, 4.5], [-1, -1, -1], [200.28857421875], 0, 1, 5826, -1, -1],
    [[-2.5, 1.5, 4.5], [-1, -1, -1], [208.0635528564453], 0, 1, 5827, -1, -1],
    [[-1.5, 1.5, 4.5], [-1, -1, -1], [211.14382934570312], 0, 1, 5828, -1, -1],
    [[-4.5, 2.5, 4.5], [-1, -1, -1], [206.84458923339844], 0, 1, 5845, -1, -1],
    [[-3.5, 2.5, 4.5], [-1, -1, -1], [213.84432983398438], 0, 1, 5846, -1, -1],
    [[-2.5, 2.5, 4.5], [-1, -1, -1], [221.49111938476562], 0, 1, 5847, -1, -1],
    [[-1.5, 2.5, 4.5], [-1, -1, -1], [224.48379516601562], 0, 1, 5848, -1, -1],
    [[-4.5, 3.5, 4.5], [-1, -1, -1], [204.69740295410156], 0, 1, 5865, -1, -1],
    [[-3.5, 3.5, 4.5], [-1, -1, -1], [211.45567321777344], 0, 1, 5866, -1, -1],
    [[-2.5, 3.5, 4.5], [-1, -1, -1], [218.91493225097656], 0, 1, 5867, -1, -1],
    [[-1.5, 3.5, 4.5], [-1, -1, -1], [221.7794189453125], 0, 1, 5868, -1, -1],
    [[-4.5, -2.5, 5.5], [-1, -1, -1], [184.473388671875], 0, 1, 6145, -1, -1],
    [[-3.5, -2.5, 5.5], [-1, -1, -1], [191.07464599609375], 0, 1, 6146, -1, -1],
    [[-2.5, -2.5, 5.5], [-1, -1, -1], [198.41197204589844], 0, 1, 6147, -1, -1],
    [[-1.5, -2.5, 5.5], [-1, -1, -1], [201.193115234375], 0, 1, 6148, -1, -1],
    [[-4.5, -1.5, 5.5], [-1, -1, -1], [206.10580444335938], 0, 1, 6165, -1, -1],
    [[-3.5, -1.5, 5.5], [-1, -1, -1], [212.86407470703125], 0, 1, 6166, -1, -1],
    [[-2.5, -1.5, 5.5], [-1, -1, -1], [220.32333374023438], 0, 1, 6167, -1, -1],
    [[-1.5, -1.5, 5.5], [-1, -1, -1], [223.18783569335938], 0, 1, 6168, -1, -1],
    [[-4.5, -0.5, 5.5], [-1, -1, -1], [206.79336547851562], 0, 1, 6185, -1, -1],
    [[-3.5, -0.5, 5.5], [-1, -1, -1], [213.63131713867188], 0, 1, 6186, -1, -1],
    [[-2.5, -0.5, 5.5], [-1, -1, -1], [221.15248107910156], 0, 1, 6187, -1, -1],
    [[-1.5, -0.5, 5.5], [-1, -1, -1], [224.05926513671875], 0, 1, 6188, -1, -1],
    [[-4.5, 0.5, 5.5], [-1, -1, -1], [188.8384552001953], 0, 1, 6205, -1, -1],
    [[-3.5, 0.5, 5.5], [-1, -1, -1], [195.67642211914062], 0, 1, 6206, -1, -1],
    [[-2.5, 0.5, 5.5], [-1, -1, -1], [203.19757080078125], 0, 1, 6207, -1, -1],
    [[-1.5, 0.5, 5.5], [-1, -1, -1], [206.10435485839844], 0, 1, 6208, -1, -1],
    [[-4.5, 1.5, 5.5], [-1, -1, -1], [185.61073303222656], 0, 1, 6225, -1, -1],
    [[-3.5, 1.5, 5.5], [-1, -1, -1], [192.36900329589844], 0, 1, 6226, -1, -1],
    [[-2.5, 1.5, 5.5], [-1, -1, -1], [199.82826232910156], 0, 1, 6227, -1, -1],
    [[-1.5, 1.5, 5.5], [-1, -1, -1], [202.69276428222656], 0, 1, 6228, -1, -1],
    [[-4.5, 2.5, 5.5], [-1, -1, -1], [199.52877807617188], 0, 1, 6245, -1, -1],
    [[-3.5, 2.5, 5.5], [-1, -1, -1], [206.13003540039062], 0, 1, 6246, -1, -1],
    [[-2.5, 2.5, 5.5], [-1, -1, -1], [213.46734619140625], 0, 1, 6247, -1, -1],
    [[-1.5, 2.5, 5.5], [-1, -1, -1], [216.24850463867188], 0, 1, 6248, -1, -1],
    [[-4.5, 3.5, 5.5], [-1, -1, -1], [197.67019653320312], 0, 1, 6265, -1, -1],
    [[-3.5, 3.5, 5.5], [-1, -1, -1], [204.04171752929688], 0, 1, 6266, -1, -1],
    [[-2.5, 3.5, 5.5], [-1, -1, -1], [211.2006378173828], 0, 1, 6267, -1, -1],
    [[-1.5, 3.5, 5.5], [-1, -1, -1], [213.85984802246094], 0, 1, 6268, -1, -1],
    [[-4.5, -2.5, 6.5], [-1, -1, -1], [175.79248046875], 0, 1, 6545, -1, -1],
    [[-3.5, -2.5, 6.5], [-1, -1, -1], [181.94105529785156], 0, 1, 6546, -1, -1],
    [[-2.5, -2.5, 6.5], [-1, -1, -1], [188.92681884765625], 0, 1, 6547, -1, -1],
    [[-1.5, -2.5, 6.5], [-1, -1, -1], [191.46768188476562], 0, 1, 6548, -1, -1],
    [[-4.5, -1.5, 6.5], [-1, -1, -1], [197.20082092285156], 0, 1, 6565, -1, -1],
    [[-3.5, -1.5, 6.5], [-1, -1, -1], [203.49728393554688], 0, 1, 6566, -1, -1],
    [[-2.5, -1.5, 6.5], [-1, -1, -1], [210.597900390625], 0, 1, 6567, -1, -1],
    [[-1.5, -1.5, 6.5], [-1, -1, -1], [213.21726989746094], 0, 1, 6568, -1, -1],
    [[-4.5, -0.5, 6.5], [-1, -1, -1], [197.774658203125], 0, 1, 6585, -1, -1],
    [[-3.5, -0.5, 6.5], [-1, -1, -1], [204.14617919921875], 0, 1, 6586, -1, -1],
    [[-2.5, -0.5, 6.5], [-1, -1, -1], [211.3050994873047], 0, 1, 6587, -1, -1],
    [[-1.5, -0.5, 6.5], [-1, -1, -1], [213.9643096923828], 0, 1, 6588, -1, -1],
    [[-4.5, 0.5, 6.5], [-1, -1, -1], [179.8197479248047], 0, 1, 6605, -1, -1],
    [[-3.5, 0.5, 6.5], [-1, -1, -1], [186.19126892089844], 0, 1, 6606, -1, -1],
    [[-2.5, 0.5, 6.5], [-1, -1, -1], [193.35018920898438], 0, 1, 6607, -1, -1],
    [[-1.5, 0.5, 6.5], [-1, -1, -1], [196.0093994140625], 0, 1, 6608, -1, -1],
    [[-4.5, 1.5, 6.5], [-1, -1, -1], [176.70574951171875], 0, 1, 6625, -1, -1],
    [[-3.5, 1.5, 6.5], [-1, -1, -1], [183.00221252441406], 0, 1, 6626, -1, -1],
    [[-2.5, 1.5, 6.5], [-1, -1, -1], [190.10284423828125], 0, 1, 6627, -1, -1],
    [[-1.5, 1.5, 6.5], [-1, -1, -1], [192.72219848632812], 0, 1, 6628, -1, -1],
    [[-4.5, 2.5, 6.5], [-1, -1, -1], [190.84786987304688], 0, 1, 6645, -1, -1],
    [[-3.5, 2.5, 6.5], [-1, -1, -1], [196.99644470214844], 0, 1, 6646, -1, -1],
    [[-2.5, 2.5, 6.5], [-1, -1, -1], [203.98220825195312], 0, 1, 6647, -1, -1],
    [[-1.5, 2.5, 6.5], [-1, -1, -1], [206.5230712890625], 0, 1, 6648, -1, -1],
    [[-4.5, 3.5, 6.5], [-1, -1, -1], [189.317138671875], 0, 1, 6665, -1, -1],
    [[-3.5, 3.5, 6.5], [-1, -1, -1], [195.24932861328125], 0, 1, 6666, -1, -1],
    [[-2.5, 3.5, 6.5], [-1, -1, -1], [202.06704711914062], 0, 1, 6667, -1, -1],
    [[-1.5, 3.5, 6.5], [-1, -1, -1], [204.49305725097656], 0, 1, 6668, -1, -1]]
    }
    """

    goldJson = json.loads(goldTestStringJson)
    ff = open("test_PhactoriGeometricCellSampler1_1.json", "r")
    testJson = json.load(ff)
    ff.close()
    self.assertEqual(goldJson, testJson)

    #remove file that got made during test
    os.remove("test_PhactoriGeometricCellSampler1_1.json")

  def test_CreateInternalListOfDataControlledSampledCellsOnThisProcess_ratio_1(self):
    testWavelet2 = Wavelet()
    testWavelet2.UpdatePipeline()
    testWavelet = PointDatatoCellData(Input=testWavelet2)
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "phactorigeometriccellsampler1"

    testOutFileBasename = "test_WriteAllDataFromOneProcessUsingMPI_output_cells_"

    operationParams = {
    "type":"geometriccellsampler1",
    "cell data array names":["RTData"],
    "cell data array tuple size": 1,
    "do programmable filter": False,
    "data controlled sampling method":"cells with ratio of min/max data value",
    "data controlled ratio of min/max":0.95,
    "sampling geometry bounding box":[-7.75, 7.75, -8.25, 9.75, -9.25, 8.215]
    }

    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'geometriccellsampler1',
              PhactoriGeometricCellSampler1,
              operationParams)

    PhactoriGeometricCellSampler1Instance = newOperationBlock.mOperationSpecifics
    PhactoriGeometricCellSampler1Instance.myCopyOfInputFilter = testWavelet
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfGeometricallySampledCellsOnThisProcess()
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfDataControlledSampledCellsOnThisProcess()

    #PhactoriGeometricCellSampler1Instance.WriteCellListToFile("test_PhactoriGeometricCellSampler1_2.json",
    #  PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess)

    self.assertEqual(len(PhactoriGeometricCellSampler1Instance.GeometricallySampledCellsForThisProcess), 4896)

    self.assertEqual(len(PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess), 66)
    firstCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[0]
    firstCellString = firstCell.ToStrTerseOneLineList()
    lastCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[-1]
    lastCellString = lastCell.ToStrTerseOneLineList()
    goldFirstCellString = "[[-1.5, -1.5, -2.5], [-1, -1, -1], [251.10580444335938], 0, 1, 2968, -1, -1]"
    goldLastCellString = "[[1.5, -0.5, 2.5], [-1, -1, -1], [253.63250732421875], 0, 1, 4991, -1, -1]"
    self.assertEqual(firstCellString, goldFirstCellString)
    self.assertEqual(lastCellString, goldLastCellString)

  def test_CreateInternalListOfDataControlledSampledCellsOnThisProcess_ratio_2(self):
    testWavelet2 = Wavelet()
    testWavelet2.UpdatePipeline()
    testWavelet = PointDatatoCellData(Input=testWavelet2)
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "phactorigeometriccellsampler1"

    testOutFileBasename = "test_WriteAllDataFromOneProcessUsingMPI_output_cells_"

    operationParams = {
    "type":"geometriccellsampler1",
    "cell data array names": ["RTData"],
    "cell data array tuple size": 1,
    "do programmable filter": False,
    "data controlled sampling method":"cells with ratio of min/max data value",
    "data controlled ratio basis":"ratio is from data minimum to data maximum",
    "data controlled ratio of min/max": 0.95,
    "sampling geometry bounding box": [-7.75, 7.75, -8.25, 9.75, -9.25, 8.215]
    }

    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'geometriccellsampler1',
              PhactoriGeometricCellSampler1,
              operationParams)

    PhactoriGeometricCellSampler1Instance = newOperationBlock.mOperationSpecifics
    PhactoriGeometricCellSampler1Instance.myCopyOfInputFilter = testWavelet
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfGeometricallySampledCellsOnThisProcess()
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfDataControlledSampledCellsOnThisProcess()

    #PhactoriGeometricCellSampler1Instance.WriteCellListToFile("test_PhactoriGeometricCellSampler1_7.json",
    #  PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess)

    self.assertEqual(len(PhactoriGeometricCellSampler1Instance.GeometricallySampledCellsForThisProcess), 4896)

    self.assertEqual(len(PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess), 40)
    firstCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[0]
    firstCellString = firstCell.ToStrTerseOneLineList()
    lastCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[-1]
    lastCellString = lastCell.ToStrTerseOneLineList()
    goldFirstCellString = "[[-0.5, -1.5, -2.5], [-1, -1, -1], [254.91671752929688], 0, 1, 2969, -1, -1]"
    goldLastCellString = "[[0.5, -0.5, 2.5], [-1, -1, -1], [254.67347717285156], 0, 1, 4990, -1, -1]"
    self.assertEqual(firstCellString, goldFirstCellString)
    self.assertEqual(lastCellString, goldLastCellString)

  def test_CreateInternalListOfDataControlledSampledCellsOnThisProcess_ratio_3(self):
    testWavelet2 = Wavelet()
    testWavelet2.UpdatePipeline()
    testWavelet = PointDatatoCellData(Input=testWavelet2)
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "phactorigeometriccellsampler1"

    testOutFileBasename = "test_WriteAllDataFromOneProcessUsingMPI_output_cells_"

    operationParams = {
    "type":"geometriccellsampler1",
    "cell data array names":["RTData"],
    "cell data array tuple size": 1,
    "do programmable filter": False,
    "data controlled sampling method":"cells with ratio of min/max data value",
    "data controlled ratio basis":"ratio is from data minimum to data maximum",
    "data controlled ratio of min/max": 0.05,
    "collect cells relative to ratio": "cells less/equal",
    "sampling geometry bounding box":[-7.75, 7.75, -8.25, 9.75, -9.25, 8.215]
    }

    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'geometriccellsampler1',
              PhactoriGeometricCellSampler1,
              operationParams)

    PhactoriGeometricCellSampler1Instance = newOperationBlock.mOperationSpecifics
    PhactoriGeometricCellSampler1Instance.myCopyOfInputFilter = testWavelet
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfGeometricallySampledCellsOnThisProcess()
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfDataControlledSampledCellsOnThisProcess()

    #PhactoriGeometricCellSampler1Instance.WriteCellListToFile("test_PhactoriGeometricCellSampler1_8.json",
    #  PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess)

    self.assertEqual(len(PhactoriGeometricCellSampler1Instance.GeometricallySampledCellsForThisProcess), 4896)

    self.assertEqual(len(PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess), 8)
    firstCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[0]
    firstCellString = firstCell.ToStrTerseOneLineList()
    lastCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[-1]
    lastCellString = lastCell.ToStrTerseOneLineList()
    goldFirstCellString = "[[-7.5, 9.5, -8.5], [-1, -1, -1], [71.15290069580078], 0, 1, 782, -1, -1]"
    goldLastCellString = "[[7.5, 9.5, 7.5], [-1, -1, -1], [77.62138366699219], 0, 1, 7197, -1, -1]"
    self.assertEqual(firstCellString, goldFirstCellString)
    self.assertEqual(lastCellString, goldLastCellString)

  def test_CreateInternalListOfDataControlledSampledCellsOnThisProcess_distance(self):
    testWavelet2 = Wavelet()
    testWavelet2.UpdatePipeline()
    testWavelet = PointDatatoCellData(Input=testWavelet2)
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "phactorigeometriccellsampler1"

    testOutFileBasename = "test_WriteAllDataFromOneProcessUsingMPI_output_cells_"

    operationParams = {
    "type":"geometriccellsampler1",
    "cell data array names":["RTData"],
    "cell data array tuple size": 1,
    "do programmable filter": False,
    "data controlled sampling method":"cells within distance of min/max highest data value cell",
    "data controlled distance":1.25,
    "data controlled sampling use min or max": "max",
    "sampling geometry bounding box":[-7.75, 7.75, -8.25, 9.75, -9.25, 8.215]
    }

    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'geometriccellsampler1',
              PhactoriGeometricCellSampler1,
              operationParams)

    PhactoriGeometricCellSampler1Instance = newOperationBlock.mOperationSpecifics
    PhactoriGeometricCellSampler1Instance.myCopyOfInputFilter = testWavelet
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfGeometricallySampledCellsOnThisProcess()
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfDataControlledSampledCellsOnThisProcess()

    #PhactoriGeometricCellSampler1Instance.WriteCellListToFile("test_PhactoriGeometricCellSampler1_3.json",
    #  PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess)

    self.assertEqual(len(PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess), 7)
    firstCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[0]
    firstCellString = firstCell.ToStrTerseOneLineList()
    lastCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[-1]
    lastCellString = lastCell.ToStrTerseOneLineList()
    goldFirstCellString = "[[-0.5, -0.5, -1.5], [-1, -1, -1], [257.593505859375], 0, 1, 3389, -1, -1]"
    goldLastCellString = "[[-0.5, -0.5, 0.5], [-1, -1, -1], [264.2397155761719], 0, 1, 4189, -1, -1]"
    self.assertEqual(firstCellString, goldFirstCellString)
    self.assertEqual(lastCellString, goldLastCellString)

    operationParams = {
    "type":"geometriccellsampler1",
    "cell data array names":["RTData"],
    "cell data array tuple size": 1,
    "do programmable filter": False,
    "data controlled sampling method":"cells within distance of min/max highest data value cell",
    "data controlled distance":2.25,
    "data controlled sampling use min or max": "min",
    "sampling geometry bounding box":[-7.75, 7.75, -8.25, 9.75, -9.25, 8.215]
    }

    PhactoriGeometricCellSampler1Instance.ParseParametersFromJson(operationParams)
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfDataControlledSampledCellsOnThisProcess()

    #PhactoriGeometricCellSampler1Instance.WriteCellListToFile("test_PhactoriGeometricCellSampler1_4.json",
    #  PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess)

    self.assertEqual(len(PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess), 17)
    firstCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[0]
    firstCellString = firstCell.ToStrTerseOneLineList()
    lastCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[-1]
    lastCellString = lastCell.ToStrTerseOneLineList()
    goldFirstCellString = "[[-7.5, 7.5, -8.5], [-1, -1, -1], [112.75462341308594], 0, 1, 742, -1, -1]"
    goldLastCellString = "[[-6.5, 9.5, -6.5], [-1, -1, -1], [96.05937194824219], 0, 1, 1583, -1, -1]"
    self.assertEqual(firstCellString, goldFirstCellString)
    self.assertEqual(lastCellString, goldLastCellString)

  def test_CreateInternalListOfDataControlledSampledCellsOnThisProcess_boundingbox(self):
    testWavelet2 = Wavelet()
    testWavelet2.UpdatePipeline()
    testWavelet = PointDatatoCellData(Input=testWavelet2)
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "phactorigeometriccellsampler1"

    testOutFileBasename = "test_WriteAllDataFromOneProcessUsingMPI_output_cells_"

    operationParams = {
    "type":"geometriccellsampler1",
    "cell data array names":["RTData"],
    "cell data array tuple size": 1,
    "do programmable filter": False,
    "data controlled sampling method":"cells within bounding box around min/max highest data value cell",
    "data controlled bounding box": [-1.25, 1.25, -2.25, 1.25, -2.25, 2.25],
    "data controlled sampling use min or max": "max",
    "sampling geometry bounding box":[-7.75, 7.75, -8.25, 9.75, -9.25, 8.215]
    }

    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'geometriccellsampler1',
              PhactoriGeometricCellSampler1,
              operationParams)

    PhactoriGeometricCellSampler1Instance = newOperationBlock.mOperationSpecifics
    PhactoriGeometricCellSampler1Instance.myCopyOfInputFilter = testWavelet
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfGeometricallySampledCellsOnThisProcess()
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfDataControlledSampledCellsOnThisProcess()

    #PhactoriGeometricCellSampler1Instance.WriteCellListToFile("test_PhactoriGeometricCellSampler1_5.json",
    #  PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess)

    operationParams = {
    "type":"geometriccellsampler1",
    "cell data array names":["RTData"],
    "cell data array tuple size": 1,
    "do programmable filter": False,
    "data controlled sampling method":"cells within bounding box around min/max highest data value cell",
    "data controlled bounding box": [-1.25, 1.25, -2.25, 1.25, -2.25, 2.25],
    "data controlled sampling use min or max": "min",
    "sampling geometry bounding box":[-7.75, 7.75, -8.25, 9.75, -9.25, 8.215]
    }

    PhactoriGeometricCellSampler1Instance.ParseParametersFromJson(operationParams)
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfDataControlledSampledCellsOnThisProcess()

    #PhactoriGeometricCellSampler1Instance.WriteCellListToFile("test_PhactoriGeometricCellSampler1_6.json",
    #  PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess)

    self.assertEqual(len(PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess), 18)
    firstCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[0]
    firstCellString = firstCell.ToStrTerseOneLineList()
    lastCell = PhactoriGeometricCellSampler1Instance.DataControlledSampledCellsForThisProcess[-1]
    lastCellString = lastCell.ToStrTerseOneLineList()
    goldFirstCellString = "[[-7.5, 7.5, -8.5], [-1, -1, -1], [112.75462341308594], 0, 1, 742, -1, -1]"
    goldLastCellString = "[[-6.5, 9.5, -6.5], [-1, -1, -1], [96.05937194824219], 0, 1, 1583, -1, -1]"
    self.assertEqual(firstCellString, goldFirstCellString)
    self.assertEqual(lastCellString, goldLastCellString)

  def test_CollectDataOnSampledCellsOnThisProcess(self):
    testWavelet2 = Wavelet()
    testWavelet2.UpdatePipeline()
    testWavelet = PointDatatoCellData(Input=testWavelet2)
    testWavelet.UpdatePipeline()
    newOperationBlock = PhactoriOperationBlock()
    newOperationBlock.mName = "phactorigeometriccellsampler1"

    testOutFileBasename = "test_WriteAllDataFromOneProcessUsingMPI_output_cells_"

    operationParams = {
    "type":"geometriccellsampler1",
    "cell data array names":["RTData"],
    "cell data array tuple size": 1,
    "do programmable filter": False,
    "data controlled sampling method":"cells within bounding box around min/max highest data value cell",
    "data controlled bounding box": [-1.25, 1.25, -2.25, 1.25, -2.25, 2.25],
    "data controlled sampling use min or max": "max",
    "sampling geometry bounding box":[-7.75, 7.75, -8.25, 9.75, -9.25, 8.215]
    }

    ParseOneFilterTypeFromViewMapOperation(newOperationBlock,
              'geometriccellsampler1',
              PhactoriGeometricCellSampler1,
              operationParams)

    PhactoriGeometricCellSampler1Instance = newOperationBlock.mOperationSpecifics
    PhactoriGeometricCellSampler1Instance.myCopyOfInputFilter = testWavelet
    PhactoriGeometricCellSampler1Instance.CreateInternalListOfGeometricallySampledCellsOnThisProcess()
    #PhactoriGeometricCellSampler1Instance.WriteCellListToFile("test_PhactoriGeometricCellSampler1_9a.json",
    #  PhactoriGeometricCellSampler1Instance.GeometricallySampledCellsForThisProcess)
    #change data on some cells then recollect it to see if it is done properly
    testData1 = []
    for oneCell in PhactoriGeometricCellSampler1Instance.GeometricallySampledCellsForThisProcess:
      testData1.append(oneCell.dataTuple[0])
      oneCell.dataTuple[0] = -1.0
    PhactoriGeometricCellSampler1Instance.CollectDataOnSampledCellsOnThisProcess()

    testData2 = []
    for oneCell in PhactoriGeometricCellSampler1Instance.GeometricallySampledCellsForThisProcess:
      testData2.append(oneCell.dataTuple[0])
    self.assertEqual(testData1, testData2)
    #PhactoriGeometricCellSampler1Instance.WriteCellListToFile("test_PhactoriGeometricCellSampler1_9b.json",
    #  PhactoriGeometricCellSampler1Instance.GeometricallySampledCellsForThisProcess)

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()
