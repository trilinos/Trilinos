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
import unittest
from Operation.PhactoriPointSourceFromJsonList import *
from paraview.simple import *

class Test_PhactoriPointSourceFromJsonList(unittest.TestCase):

  def test_ValidateJsonPointList(self):
    theSource = PhactoriPointSourceFromJsonList()
    theSource.JsonList = [
      [0.0, 0.0, 0.0],
      [1.0, 2.0, 3.0],
      [-1.0, -2.0, -3.0],
      [5.75, -32.125, 17.625]
    ]
    self.assertTrue(theSource.ValidateJsonPointList())

    theSource.JsonList = [
      [0.0, 0.0, 0.0],
      [1.0, 2.0, 3.0],
      [-1.0, -2.0],
      [5.75, -32.125, 17.625]
    ]
    with self.assertRaises(Exception):
      theSource.ValidateJsonPointList()

    theSource.JsonList = []
    with self.assertRaises(Exception):
      theSource.ValidateJsonPointList()


  def test_CreateVtkPolyDataFromJsonList(self):
    theSource = PhactoriPointSourceFromJsonList()
    theSource.JsonList = [
      [0.0, 0.0, 0.0],
      [1.0, 2.0, 3.0],
      [-1.0, -2.0, -3.0],
      [5.75, -32.125, 17.625]
    ]
    theSource.CreateVtkPolyDataFromJsonList()

    vtkPts = theSource.myVtkPolyData.GetPoints()
    self.assertEqual(4, vtkPts.GetNumberOfPoints())
    pt0 = vtkPts.GetPoint(0)
    self.assertEqual(pt0, tuple(theSource.JsonList[0]))
    pt3 = vtkPts.GetPoint(3)
    self.assertEqual(pt3, tuple(theSource.JsonList[3]))

if __name__ == '__main__':
    cc = Cone()
    rr = Show()
    unittest.main()
