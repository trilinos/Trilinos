C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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

C=======================================================================
      SUBROUTINE EXPXYZN (XN, YN, ZN, XEXPL, YEXPL, ZEXPL,
     &   NUMNPS, IDNPS, NNNPS, IXNNPS, LTNNPS, NUMNP, NDIM,
     $     MODE)
C=======================================================================

C   --*** EXPXYZ *** (GREPOS) Modify coordinates for each block separately
C   --   Written by Amy Gilkey - revised 05/09/88
C   --   Modified by Greg Sjaardema - 02/06/89
C   --
C   --EXPXYZ modifies the coordinate array for the database.
C   --       each block is treated separately if not connected
C   --
C   --Parameters:
C   --   XN, YN, ZN - OUT - the coordinates
C   --   MODE  - 1 = Explode
C   --           2 = Scale
C   --           3 = Randomize
C   --           4 = Node Randomize

      REAL XN(*), YN(*), ZN(*)
      REAL XEXPL(*), YEXPL(*), ZEXPL(*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      INTEGER MODE

      if (MODE .NE. 4) RETURN

C ... Randomize Each Nodeset node
      IDUM = 1
      DO 150 I = 1, NUMNPS
        if (XEXPL(I) .ne. 0.0 .or. YEXPL(I) .ne. 0.0 .or.
     *    ZEXPL(I) .ne. 0.0) THEN
          DO 140 INOD = IXNNPS(I), IXNNPS(I)+nnnps(i)-1
            NODE = LTNNPS(INOD)
            XN(NODE) = (2.0*RAN1(IDUM)-1.0) * XEXPL(I) + XN(NODE)
            YN(NODE) = (2.0*RAN1(IDUM)-1.0) * YEXPL(I) + YN(NODE)
            IF (NDIM .EQ. 3) THEN
               ZN(NODE) = (2.0*RAN1(IDUM)-1.0) * ZEXPL(I) + ZN(NODE)
             END IF
 140       CONTINUE
         END IF
 150   CONTINUE
      RETURN
      END
