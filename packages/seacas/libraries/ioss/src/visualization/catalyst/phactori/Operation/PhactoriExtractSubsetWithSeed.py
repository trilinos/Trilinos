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
from .PhactoriParallelGeometryUtilities import *

#phactori_combine_to_single_python_file_subpiece_begin_1

class PhactoriExtractSubsetWithSeed(PhactoriOperationSpecifics):
  """manages extract block filter, including creating flat block
     indices list from list of block names"""
  def __init__(self):
    PhactoriOperationSpecifics.__init__(self)
    self.Seed = [0.0, 0.0, 0.0]
    self.Direction = "I"
    self.FindSeedNearestTestPoint = False
    self.TestPointToFindSeed = [0.0, 0.0, 0.0]
    self.SeedOffset = [0.0, 0.0, 0.0]
    self.SeedWithOffset = [0.0, 0.0, 0.0]

  def ParseParametersFromJson(self, inJson):
    key1 = "seed"
    key2 = "seed from cell nearest point"
    if "seed" in inJson:
      if key2 in inJson:
        myDebugPrint3AndException(
          "PhactoriExtractSubsetWithSeed.ParseParametersFromJson exception:\n"
          "block may have one of the following two tokens, but not both:\n"
          "'" + key1 + "'\n"
          "'" + key2 + "'\n")
      self.Seed = inJson["seed"]
      if len(self.Seed) != 3:
        myDebugPrint3AndException(
          "PhactoriExtractSubsetWithSeed.ParseParametersFromJson exception:\n"
          "'" + key1 + "' token must be [x,y,z] (list of 3 floats)\n")

    if key2 in inJson:
      self.FindSeedNearestTestPoint = True
      self.TestPointToFindSeed = inJson[key2]
      if len(self.TestPointToFindSeed) != 3:
        myDebugPrint3AndException(
          "PhactoriExtractSubsetWithSeed.ParseParametersFromJson exception:\n"
          "'" + key2 + "' token must be [x,y,z] (list of 3 floats)\n")

    key3 = "seed offset"
    if key3 in inJson:
      self.SeedOffset = inJson[key3]
      if len(self.SeedOffset) != 3:
        myDebugPrint3AndException(
          "PhactoriExtractSubsetWithSeed.ParseParametersFromJson exception:\n"
          "'" + key3 + "' token must be [x,y,z] (list of 3 floats)\n")

    if "direction" in inJson:
      self.Direction = inJson["direction"]
      validDirectionStrings = ["I", "J", "K", "IJ", "JK", "KI"]
      if self.Direction not in validDirectionStrings:
        myDebugPrint3AndException(
          "PhactoriExtractSubsetWithSeed.ParseParametersFromJson exception:\n"
          "illegal 'direction' token: " + str(self.Direction) + "\n"
          "the direction must one of: " + str(validDirectionStrings) + "\n")

  def CreateParaViewFilter(self, inInputFilter):
    """create the ExtractSubsetWithSeed filter for ParaView"""
    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractSubsetWithSeed.CreateParaViewFilter "
          "entered\n", 100)

    savedActiveSource = GetActiveSource()

    UpdatePipelineWithCurrentTimeArgument(inInputFilter)

    newParaViewFilter = ExtractSubsetWithSeed(inInputFilter)

    if self.FindSeedNearestTestPoint:
      ##nearestCellTestPointList = GetListOfCellTestPointsNearestListOfPointsV5(inInputFilter, [self.TestPointToFindSeed])
      nearestCellTestPointList = GetListOfGridPointsNearestListOfPointsV5(inInputFilter, [self.TestPointToFindSeed])
      self.Seed = nearestCellTestPointList[0]
      if PhactoriDbg(100):
        myDebugPrint3("find grid point for seed nearest this point:\n" + str(self.TestPointToFindSeed) + "\n"
          "seed point found:\n" + str(self.Seed) + "\n")

    self.SeedWithOffset[0] = self.Seed[0] + self.SeedOffset[0]
    self.SeedWithOffset[1] = self.Seed[1] + self.SeedOffset[1]
    self.SeedWithOffset[2] = self.Seed[2] + self.SeedOffset[2]
    if PhactoriDbg(100):
      myDebugPrint3("self.Seed:           " + str(self.Seed) + "\n"
                    "self.SeedWithOffset: " + str(self.SeedWithOffset) + "\n")
    newParaViewFilter.Seed = self.SeedWithOffset
    newParaViewFilter.Direction = self.Direction

    SetActiveSource(newParaViewFilter)

    UpdatePipelineWithCurrentTimeArgument(newParaViewFilter)

    SetActiveSource(savedActiveSource)

    if PhactoriDbg(100):
      myDebugPrint3("PhactoriExtractSubsetWithSeed.CreateParaViewFilter "
          "returning\n", 100)

    return newParaViewFilter
#phactori_combine_to_single_python_file_subpiece_end_1
