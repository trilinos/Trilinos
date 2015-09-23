C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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

C $Log: elestf.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:00:19  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:48  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE ELESTF (NLNKF, LINKF1, XN, YN, ZN)
C=======================================================================

C   --*** ELESTF *** (DETOUR) Plot element symbol for state for a face
C   --   Written by Amy Gilkey, revised 10/21/87
C   --   D. P. Flanagan, 12/08/83
C   --
C   --ELESTF plots the element symbol for the state for a face.
C   --
C   --Parameters:
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates

      INTEGER LINKF1(*)
      INTEGER NLNKF
      REAL XN(*), YN(*), ZN(*)

      N2 = LINKF1(NLNKF)
      DO 100 ILINK = 1, NLNKF
         N1 = LINKF1(ILINK)
         CALL MPD2VC (1, XN(N1), YN(N1), XN(N2), YN(N2))
         N2 = N1
  100 CONTINUE
      IF (NLNKF .EQ. 4) THEN
         N1 = LINKF1(1)
         N2 = LINKF1(3)
         CALL MPD2VC (1, XN(N1), YN(N1), XN(N2), YN(N2))
         N1 = LINKF1(2)
         N2 = LINKF1(4)
         CALL MPD2VC (1, XN(N1), YN(N1), XN(N2), YN(N2))
      END IF

      RETURN
      END
