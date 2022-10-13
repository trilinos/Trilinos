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

#phactori_combine_to_single_python_file_subpiece_begin_1
class PhactoriOutputFileNameAndDirectoryManager():
  def __init__(self, defaultBasename, defaultBasedirectory, defaultNumCounterDigits, defaultExtension):
    self.Basename = defaultBasename
    self.Basedirectory = defaultBasedirectory
    self.NumCounterDigits = defaultNumCounterDigits
    self.Extension = defaultExtension

  def ParseOutputFileNameControl(self, inJson, keyPrefix):
    key1 = keyPrefix + "basename"
    if key1 in inJson:
      self.Basename = inJson[key1]
    key1 = keyPrefix + "basedirectory"
    if key1 in inJson:
      self.Basedirectory = inJson[key1]
    key1 = keyPrefix + "counterdigits"
    if key1 in inJson:
      self.NumCounterDigits = inJson[key1]
    key1 = keyPrefix + "extension"
    if key1 in inJson:
      self.Extension = inJson[key1]

  def GetOutputFilePath(self, counterValue, extraBasename):
    if counterValue < 0:
      counterValue = 0
    retName = self.Basedirectory + "/" + self.Basename + extraBasename + \
      str(counterValue).zfill(self.NumCounterDigits) + self.Extension
    return retName
#phactori_combine_to_single_python_file_subpiece_end_1
