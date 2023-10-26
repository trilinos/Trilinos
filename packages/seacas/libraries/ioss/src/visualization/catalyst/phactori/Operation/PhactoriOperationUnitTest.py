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

if True:
  from Operation.test_PhactoriVtkDataExportOperation import *
  from Operation.test_PhactoriSegmentCellSampler3 import *
  from Operation.test_PhactoriSegment import *
  from Operation.test_PhactoriPointSourceFromJsonList import *
  from Operation.test_PhactoriSampledCellInfo import *
  from Operation.test_PhactoriGeometricCellSampler1 import *
  from Operation.test_PhactoriVtkCellOperations import *
  from Operation.test_PhactoriVariableInfo import *
  from Operation.test_PhactoriDataArtifactMetaDataControl import *
  from Operation.test_PhactoriCreateSegmentsNormalToCells import *
  from Operation.test_PhactoriMarkCellSurfaceStatus2 import *
  from Operation.test_PhactoriSliceWithPlaneOperation import *
  from Operation.test_PhactoriExtractBlockOperation import *
  from Operation.test_PhactoriThresholdOperation import *
  from Operation.test_PhactoriContourOperation import *
  from Operation.test_PhactoriClipPlaneOperation import *
  from Operation.test_PhactoriExtractSurfaceOperation import *
  from Operation.test_PhactoriCellDataToPointDataOperation import *
  from Operation.test_PhactoriPointDataToCellDataOperation import *
  from Operation.test_PhactoriMergeBlocksOperation import *
if False:
  print("if statement for temporarily moving tests during development")

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()


