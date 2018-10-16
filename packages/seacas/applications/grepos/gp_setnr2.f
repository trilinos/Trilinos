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
      SUBROUTINE SETNR2 (ISPLANE, NUMNP, X, Y, IDIM, VNORM, NRMTYP,
     *  PRJVEC, NEESS, NNESS, LTNESS)
C=======================================================================

      include 'gp_snap.blk'

      LOGICAL ISPLANE
      INTEGER NRMTYP

      REAL X(*), Y(*), VNORM(IDIM,*), PRJVEC(3)
      DIMENSION LTNESS(2,*)

      IF (NRMTYP .EQ. PNORM .OR. ISPLANE) THEN
        DO 20 ISEG = 1, NEESS

            XI = x(LTNESS(1,ISEG))
            YI = y(LTNESS(1,ISEG))

            XJ = x(LTNESS(2,ISEG))
            YJ = y(LTNESS(2,ISEG))

            DX = XI - XJ
            DY = YI - YJ
            RMAG = SQRT ( DX**2 + DY**2)
C
            AI = -dy / rmag
            BJ =  dx / rmag

          if (isplane) then
            vnorm(1, iseg) = ai
            vnorm(2, iseg) = bj
            vnorm(3, iseg) = 0.0
            di = ai * xi + bj * yi
            dj = ai * xj + bj * yj
            vnorm(4,iseg) = (di + dj) / 2.0
          else
            DO 10 INOD = 1, 2
              NODE = LTNESS(INOD, ISEG)
              VNORM(1, NODE) = VNORM(1, NODE) + AI
              VNORM(2, NODE) = VNORM(2, NODE) + BJ
 10         CONTINUE
          end if
 20     CONTINUE
      
        if (isplane) return

C ... NORMALIZE PRJVECS
        do 30 i=1, neess
          do 25 j=1,2
            node = ltness(j, i)
            A = VNORM(1,NODE)
            B = VNORM(2,NODE)
            R = A**2 + B**2
            IF (R .GT. 0.0) THEN
              R = SQRT(R)
              VNORM(1,NODE) = A/R
              VNORM(2,NODE) = B/R
            END IF
 25       continue
 30     CONTINUE

      ELSE IF (NRMTYP .EQ. PRAD) THEN
C ... Radial Projection
        do 70 i=1, neess
          do 60 j=1,2
            node = ltness(j, i)
            a = x(node) - prjvec(1)
            b = y(node) - prjvec(2)
            R = sqrt(a**2 + b**2)
            if (r .eq. 0) r = 1.0
            vnorm(1,node) = a/r
            vnorm(2,node) = b/r
 60       continue
 70     continue

      ELSE IF (NRMTYP .EQ. PVECT) THEN
C ... Use user-supplied normal
        R = SQRT(PRJVEC(1)**2 + PRJVEC(2)**2)
        R = SQRT(R)
        PRJVEC(1) = PRJVEC(1)/R
        PRJVEC(2) = PRJVEC(2)/R
        do 50 i=1, neess
          do 40 j=1,2
            node = ltness(j, i)
            vnorm(1,node) = prjvec(1) 
            vnorm(2,node) = prjvec(2) 
 40       continue
 50     continue
        
      ELSE IF (NRMTYP .EQ. PEDGE) THEN
        
      END IF
      RETURN
      END

