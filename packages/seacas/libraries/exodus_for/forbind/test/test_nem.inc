C Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C==================================================
C Parameters for the NE_TEST fortran interface test
C program
C==================================================

C Constants for init_global functions
      PARAMETER (NNG = 100)
      PARAMETER (NEG = 50)
      PARAMETER (NEBG = 5)
      PARAMETER (NNSG = 2)
      PARAMETER (NSSG = 3)

C Constants for load-balance functions
      PARAMETER (NPROC = 10)
      PARAMETER (NPROCF = NPROC)
      PARAMETER (NINTN = 200)
      PARAMETER (NBORN = 10)
      PARAMETER (NEXTN = 5)
      PARAMETER (NINTE = 100)
      PARAMETER (NBORE = 20)
      INTEGER NNCMAP, NECMAP
      PARAMETER (NNCMAP = 4)
      PARAMETER (NECMAP = 2)

C Constants for communication map
      PARAMETER (NCNTCM = 20)
      INTEGER ECNTCM
      PARAMETER (ECNTCM = 17)

C Some multiplied values
C   NNCMAP * NPROCF
      PARAMETER (NNCXNPF = 40)
C   NECMAP * NPROCF
      PARAMETER (NECXNPF = 20)
