C Copyright (c) 2014, Sandia Corporation.
C Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C the U.S. Governement retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 
C
C     ERROR VECTOR AND FLAGS.
C     THE ERROR PARAMETERS BELONG IN SUBROUTINE VTABLE ALSO.
C     CHANGES IN THE ERROR FLAGS MUST BE ACCOMPANIED BY CHANGES
C     IN THE ERROR FORMAT STATEMENT.
C
      PARAMETER (NERVEC=18)
      PARAMETER (SUCESS=1, NOGET=2, NOINIT=3, REINIT=4, NONAME=5, 
     *   BDNAME=6, BADLEN=7, BDTYPE=8, DFULL=9, VFULL=10, BFULL=11, 
     *   BDVOID=12, BDBLK=13, ILBLK=14, BADCOD=15, BADNAM=16,
     *   DEFRON=17, WRTYPE=18, BADOFF=18)
C
