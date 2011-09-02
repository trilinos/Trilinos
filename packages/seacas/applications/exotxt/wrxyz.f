C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
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


C=======================================================================
      SUBROUTINE WRXYZ (NTXT, NDIM, NUMNP, XN, YN, ZN, nameco, namlen)
C=======================================================================

C   --*** WRXYZ *** (EXOTXT) Write database coordinates
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --WRXYZ writes the coordinate array from the database.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NDIM - IN - the number of coordinates per node
C   --   NUMNP - IN - the number of nodes
C   --   XN, YN, ZN - IN - the coordinates
C   --
C   --Database must be positioned at start of coordinates upon entry;
C   --upon exit at end of coordinates.

      REAL XN(*), YN(*), ZN(*)
      character*(namlen) nameco(*)

      WRITE (NTXT, '(A)')
     &   '! Coordinate names'
      write (ntxt, 10020) (nameco(i),i=1,ndim)

      WRITE (NTXT, '(A)') '! Coordinates'

      DO 100 INP = 1, NUMNP
         IF (NDIM .EQ. 2) THEN
            WRITE (NTXT, 10000) XN(INP), YN(INP)
         ELSE IF (NDIM .EQ. 3) THEN
            WRITE (NTXT, 10000) XN(INP), YN(INP), ZN(INP)
         END IF
  100 CONTINUE

      RETURN
10000  FORMAT (5(1pE16.7))
10020  FORMAT (3(A,1x))
       END
