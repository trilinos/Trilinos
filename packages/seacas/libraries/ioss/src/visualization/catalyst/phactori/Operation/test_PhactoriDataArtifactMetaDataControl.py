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

import unittest
from Operation.PhactoriDataArtifactMetaDataControl import *
import paraview.simple
import subprocess

class TestPhactoriDataArtifactMetaDataControl(unittest.TestCase):

  def test_PhactoriDataArtifactMetaDataControl(self):
    test_pdamdc = PhactoriDataArtifactMetaDataControl(0, 1)
    test_pdamdc.AddImageToDataArtifactOutputList("CatalystOutput/test1.png")
    test_pdamdc.AddDataExportToDataArtifactOutputList("CatalystVtkDataOutput/test2.vtm")
    test_pdamdc.CloseFiles()
    ff = open(test_pdamdc.artifactListCsvFileName, "r")
    self.assertNotEqual(ff, None)
    fflines = ff.readlines()
    ff.close()
    self.assertEqual(len(fflines), 2)
    ln1 = fflines[0].split(",")
    self.assertEqual(len(ln1), 2)
    self.assertEqual(ln1[0], "CatalystOutput/test1.png")
    ln1 = fflines[1].split(",")
    self.assertEqual(len(ln1), 2)
    self.assertEqual(ln1[0], "CatalystVtkDataOutput/test2.vtm")
    subprocess.run(["rm", test_pdamdc.artifactListCsvFileName])

  def test_PhactoriDataArtifactMetaDataControl_nonwrite_process(self):
    test_pdamdc = PhactoriDataArtifactMetaDataControl(0, 8)
    test_pdamdc.AddImageToDataArtifactOutputList("CatalystOutput/test1.png")
    test_pdamdc.AddDataExportToDataArtifactOutputList("CatalystVtkDataOutput/test2.vtm")
    self.assertEqual(test_pdamdc.artifactListCsvFileName, None)
    self.assertEqual(test_pdamdc.artifactListCsvFileHandle, None)

  def test_PhactoriDataArtifactMetaDataControl_output_enable(self):
    test_pdamdc = PhactoriDataArtifactMetaDataControl(0, 1)
    test_pdamdc.EnableOutput(0)
    test_pdamdc.AddImageToDataArtifactOutputList("CatalystOutput/test1.png")
    test_pdamdc.AddDataExportToDataArtifactOutputList("CatalystVtkDataOutput/test2.vtm")
    self.assertEqual(test_pdamdc.artifactListCsvFileName, None)
    self.assertEqual(test_pdamdc.artifactListCsvFileHandle, None)
    test_pdamdc.EnableOutput(1)
    test_pdamdc.AddImageToDataArtifactOutputList("CatalystOutput/test1.png")
    test_pdamdc.AddDataExportToDataArtifactOutputList("CatalystVtkDataOutput/test2.vtm")
    test_pdamdc.CloseFiles()
    ff = open(test_pdamdc.artifactListCsvFileName, "r")
    self.assertNotEqual(ff, None)
    fflines = ff.readlines()
    ff.close()
    self.assertEqual(len(fflines), 2)
    ln1 = fflines[0].split(",")
    self.assertEqual(len(ln1), 2)
    self.assertEqual(ln1[0], "CatalystOutput/test1.png")
    ln1 = fflines[1].split(",")
    self.assertEqual(len(ln1), 2)
    self.assertEqual(ln1[0], "CatalystVtkDataOutput/test2.vtm")
    subprocess.run(["rm", test_pdamdc.artifactListCsvFileName])

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()


