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

import math

#phactori_combine_to_single_python_file_subpiece_begin_1
def vecDotProduct(inVecA, inVecB):
  """ returns dotproct inVecZ dot inVecB """
  return inVecA[0] * inVecB[0] + inVecA[1] * inVecB[1] + inVecA[2] * inVecB[2]

def vecCrossProduct(inVecA, inVecB):
  """ returns inVecZ X inVecB """
  return [inVecA[1] * inVecB[2] - inVecA[2] * inVecB[1],
          inVecA[2] * inVecB[0] - inVecA[0] * inVecB[2],
          inVecA[0] * inVecB[1] - inVecA[1] * inVecB[0]]

def vecCrossProduct2(outVec, inVecA, inVecB):
  """ calculates inVecA X inVecB and puts the result in outVec"""
  outVec[0] = inVecA[1] * inVecB[2] - inVecA[2] * inVecB[1]
  outVec[1] = inVecA[2] * inVecB[0] - inVecA[0] * inVecB[2]
  outVec[2] = inVecA[0] * inVecB[1] - inVecA[1] * inVecB[0]

def vecCopy(destinationVec, sourceVec):
  destinationVec[0] = sourceVec[0]
  destinationVec[1] = sourceVec[1]
  destinationVec[2] = sourceVec[2]

def vecMagnitude(inVec):
  xx = inVec[0]
  yy = inVec[1]
  zz = inVec[2]
  return math.sqrt(xx*xx + yy*yy + zz*zz)

def vecMagnitudeSquared(inVec):
  xx = inVec[0]
  yy = inVec[1]
  zz = inVec[2]
  return xx*xx + yy*yy + zz*zz

def vecNormalizeWithSmallCheck(inVec, smallValue, smallNormal):
  xx = inVec[0]
  yy = inVec[1]
  zz = inVec[2]
  mag = math.sqrt(xx*xx + yy*yy + zz*zz)
  if mag <= smallValue:
    return smallNormal
  else:
    return [xx/mag, yy/mag, zz/mag]

def vecNormalize(inVec):
  xx = inVec[0]
  yy = inVec[1]
  zz = inVec[2]
  mag = math.sqrt(xx*xx + yy*yy + zz*zz)
  return [xx/mag, yy/mag, zz/mag]

def vecNormalize2(outVec, inVec):
  xx = inVec[0]
  yy = inVec[1]
  zz = inVec[2]
  mag = math.sqrt(xx*xx + yy*yy + zz*zz)
  outVec[0] = xx/mag
  outVec[1] = yy/mag
  outVec[2] = zz/mag

def vecFromAToB(inVecA, inVecB):
  return [inVecB[0] - inVecA[0], inVecB[1] - inVecA[1], inVecB[2] - inVecA[2]]

def vecDistanceSquared(inPtA, inPtB):
  """calculate the square of the distance between two points"""
  ddx = inPtA[0] - inPtB[0]
  ddy = inPtA[1] - inPtB[1]
  ddz = inPtA[2] - inPtB[2]
  return ddx*ddx + ddy*ddy + ddz*ddz

def vecDistance(inPtA, inPtB):
  """calculate the distance between two points"""
  return math.sqrt(vecDistanceSquared(inPtA, inPtB))

def vecAdd(inVecA, inVecB):
  return [inVecA[0]+inVecB[0],inVecA[1]+inVecB[1],inVecA[2]+inVecB[2]]

def vecAddInPlace(inVecA, inVecB):
  """changes inVecA to inVecA + inVecB"""
  inVecA[0] += inVecB[0]
  inVecA[1] += inVecB[1]
  inVecA[2] += inVecB[2]

def vecScale(inScale, inVec):
  return [inScale*inVec[0], inScale*inVec[1], inScale*inVec[2]]

def vecScaleInPlace(inScale, inVec):
  inVec[0] *= inScale
  inVec[1] *= inScale
  inVec[2] *= inScale

def vecMultiplyAdd(inVecA, inVecB, inMM):
  """returns inVecA + (inMM * inVecB)"""
  return [inVecA[0] + inMM * inVecB[0],
          inVecA[1] + inMM * inVecB[1],
          inVecA[2] + inMM * inVecB[2]]

def GetPhactoriVectorLibraryProgrammableFilterLines(pfLns):
  pfLns.append("def vecDotProduct(inVecA, inVecB):\n")
  pfLns.append("  return inVecA[0] * inVecB[0] + inVecA[1] * inVecB[1] + inVecA[2] * inVecB[2]\n")
  pfLns.append("def vecCrossProduct(inVecA, inVecB):\n")
  pfLns.append("  return [inVecA[1] * inVecB[2] - inVecA[2] * inVecB[1],\n")
  pfLns.append("          inVecA[2] * inVecB[0] - inVecA[0] * inVecB[2],\n")
  pfLns.append("          inVecA[0] * inVecB[1] - inVecA[1] * inVecB[0]]\n")
  pfLns.append("def vecCrossProduct2(outVec, inVecA, inVecB):\n")
  pfLns.append("  outVec[0] = inVecA[1] * inVecB[2] - inVecA[2] * inVecB[1]\n")
  pfLns.append("  outVec[1] = inVecA[2] * inVecB[0] - inVecA[0] * inVecB[2]\n")
  pfLns.append("  outVec[2] = inVecA[0] * inVecB[1] - inVecA[1] * inVecB[0]\n")
  pfLns.append("def vecCopy(destinationVec, sourceVec):\n")
  pfLns.append("  destinationVec[0] = sourceVec[0]\n")
  pfLns.append("  destinationVec[1] = sourceVec[1]\n")
  pfLns.append("  destinationVec[2] = sourceVec[2]\n")
  pfLns.append("def vecMagnitude(inVec):\n")
  pfLns.append("  xx = inVec[0]\n")
  pfLns.append("  yy = inVec[1]\n")
  pfLns.append("  zz = inVec[2]\n")
  pfLns.append("  return math.sqrt(xx*xx + yy*yy + zz*zz)\n")
  pfLns.append("def vecMagnitudeSquared(inVec):\n")
  pfLns.append("  xx = inVec[0]\n")
  pfLns.append("  yy = inVec[1]\n")
  pfLns.append("  zz = inVec[2]\n")
  pfLns.append("  return xx*xx + yy*yy + zz*zz\n")
  pfLns.append("def vecNormalizeWithSmallCheck(inVec, smallValue, smallNormal):\n")
  pfLns.append("  xx = inVec[0]\n")
  pfLns.append("  yy = inVec[1]\n")
  pfLns.append("  zz = inVec[2]\n")
  pfLns.append("  mag = math.sqrt(xx*xx + yy*yy + zz*zz)\n")
  pfLns.append("  if mag <= smallValue:\n")
  pfLns.append("    return smallNormal\n")
  pfLns.append("  else:\n")
  pfLns.append("    return [xx/mag, yy/mag, zz/mag]\n")
  pfLns.append("def vecNormalize(inVec):\n")
  pfLns.append("  xx = inVec[0]\n")
  pfLns.append("  yy = inVec[1]\n")
  pfLns.append("  zz = inVec[2]\n")
  pfLns.append("  mag = math.sqrt(xx*xx + yy*yy + zz*zz)\n")
  pfLns.append("  return [xx/mag, yy/mag, zz/mag]\n")
  pfLns.append("def vecNormalize2(outVec, inVec):\n")
  pfLns.append("  xx = inVec[0]\n")
  pfLns.append("  yy = inVec[1]\n")
  pfLns.append("  zz = inVec[2]\n")
  pfLns.append("  mag = math.sqrt(xx*xx + yy*yy + zz*zz)\n")
  pfLns.append("  outVec[0] = xx/mag\n")
  pfLns.append("  outVec[1] = yy/mag\n")
  pfLns.append("  outVec[2] = zz/mag\n")
  pfLns.append("def vecFromAToB(inVecA, inVecB):\n")
  pfLns.append("  return [inVecB[0] - inVecA[0], inVecB[1] - inVecA[1], inVecB[2] - inVecA[2]]\n")
  pfLns.append("def vecDistanceSquared(inPtA, inPtB):\n")
  pfLns.append("  ddx = inPtA[0] - inPtB[0]\n")
  pfLns.append("  ddy = inPtA[1] - inPtB[1]\n")
  pfLns.append("  ddz = inPtA[2] - inPtB[2]\n")
  pfLns.append("  return ddx*ddx + ddy*ddy + ddz*ddz\n")
  pfLns.append("def vecDistance(inPtA, inPtB):\n")
  pfLns.append("  return math.sqrt(vecDistanceSquared(inPtA, inPtB))\n")
  pfLns.append("def vecAdd(inVecA, inVecB):\n")
  pfLns.append("  return [inVecA[0]+inVecB[0],inVecA[1]+inVecB[1],inVecA[2]+inVecB[2]]\n")
  pfLns.append("def vecScale(inScale, inVec):\n")
  pfLns.append("  return [inScale*inVec[0], inScale*inVec[1], inScale*inVec[2]]\n")
  pfLns.append("def vecMultiplyAdd(inVecA, inVecB, inMM):\n")
  pfLns.append("  return [inVecA[0] + inMM * inVecB[0],\n")
  pfLns.append("          inVecA[1] + inMM * inVecB[1],\n")
  pfLns.append("          inVecA[2] + inMM * inVecB[2]]\n")

#phactori_combine_to_single_python_file_subpiece_end_1
