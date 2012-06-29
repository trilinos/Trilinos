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

C=======================================================================
      SUBROUTINE SCAPRT (NAMECO, NAME, IVAR, TIMES,
     &   VALMIN, NUMMIN, XYZMIN, ISTMIN, VALMAX, NUMMAX, XYZMAX, ISTMAX)
C=======================================================================

C   --*** SCAPRT *** (BLOT) Print variable min/max
C   --   Written by Amy Gilkey - revised 04/01/88
C   --
C   --SCAPRT prints the minimum and maximum values for the desired variable.
C   --
C   --Parameters:
C   --   NAMECO - IN - the coordinate names
C   --   NAME - IN - the variable name
C   --   IVAR - IN - the variable index
C   --   TIMES - IN - the database times
C   --   VALMIN, VALMAX - IN - the minimum and maximum value
C   --      (with selected element block and birth/death)
C   --   NUMMIN, NUMMAX - IN - the node or element number of the minimum and
C   --      maximum value (nodal and element variable only)
C   --   XYZMIN, XYZMAX - IN - the coordinates of NUMMIN, NUMMAX
C   --   ISTMIN, ISTMAX - IN - the step number of the minimum and maximum
C   --      value
C   --
C   --Common Variables:
C   --   Uses NDIM of /DBNUMS/

      include 'dbnums.blk'

      CHARACTER*(*) NAMECO(*)
      CHARACTER*(*) NAME
      REAL TIMES(*)
      REAL XYZMIN(3), XYZMAX(3)

      REAL RVAL(2), RXYZ(2,3), RTIM(2)
      CHARACTER*20 SVAL(0:2), SXYZ0(3), SXYZ(2,3), STIM(0:2)
      CHARACTER*4 STR4
      CHARACTER TYP

      IF ((ISTMIN .LE. 0) .OR. (VALMIN .EQ. VALMAX)) THEN
         CALL NUMSTR (1, 4, VALMIN, SVAL(1), LVAL)
         WRITE (*, 10060) NAME(:LENSTR(NAME)),
     &      ' does not vary - all values = ', SVAL(1)(:LVAL)

      ELSE
         CALL DBVTYP (IVAR, TYP, IDUM)

         RVAL(1) = VALMIN
         RVAL(2) = VALMAX
         CALL NUMSTR (2, 4, RVAL, SVAL(1), LVAL)
         CALL CENLAB ('Value', 2, SVAL(1), LVAL, SVAL(0))

         IF ((TYP .EQ. 'N') .OR. (TYP .EQ. 'E')) THEN
            NNDIM = MAX (MIN (NDIM, 3), 2)
            DO 100 I = 1, MIN (NDIM, 3)
               RXYZ(1,I) = XYZMIN(I)
               RXYZ(2,I) = XYZMAX(I)
  100       CONTINUE
            CALL NUMSTR (NNDIM * 2, 3, RXYZ, SXYZ(1,1), LX)
            LY = LX
            LZ = LX
            CALL CENLAB (NAMECO(1), 2, SXYZ(1,1), LX, SXYZ0(1))
            CALL CENLAB (NAMECO(2), 2, SXYZ(1,2), LY, SXYZ0(2))
            IF (NDIM .GE. 3) THEN
               CALL CENLAB (NAMECO(3), 2, SXYZ(1,3), LZ, SXYZ0(3))
            END IF
         END IF

         RTIM(1) = TIMES(ISTMIN)
         RTIM(2) = TIMES(ISTMAX)
         CALL NUMSTR (2, 4, RTIM, STIM(1), LTIM)
         CALL CENLAB ('Time', 2, STIM(1), LTIM, STIM(0))

         IF ((TYP .EQ. 'N') .OR. (TYP .EQ. 'E')) THEN
            IF (TYP .EQ. 'N') THEN
               STR4 = 'Node'
            ELSE
               STR4 = 'Elem'
            END IF
            IF (NDIM .LE. 2) THEN
               WRITE (*, 10000) 'Range: ', SVAL(0)(:LVAL),
     &            SXYZ0(1)(:LX), SXYZ0(2)(:LY),
     &            STIM(0)(:LTIM), STR4, 'Step'
               WRITE (*, 10010, IOSTAT=IDUM) 'Minimum', SVAL(1)(:LVAL),
     &            SXYZ(1,1)(:LX), SXYZ(1,2)(:LY),
     &            STIM(1)(:LTIM), NUMMIN, ISTMIN
               WRITE (*, 10010, IOSTAT=IDUM) 'Maximum', SVAL(2)(:LVAL),
     &            SXYZ(2,1)(:LX), SXYZ(2,2)(:LY),
     &            STIM(2)(:LTIM), NUMMAX, ISTMAX
10000           FORMAT (4X, A7, 1X, A, 3X, 2 (1X, A),
     &            4X, A, 3X, A4, 3X, A4)
10010           FORMAT (4X, A7, 1X, A, 3X, 2 (1X, A), 4X, A, I7, I7)
            ELSE
               WRITE (*, 10020) 'Range: ', SVAL(0)(:LVAL),
     &            SXYZ0(1)(:LX), SXYZ0(2)(:LY), SXYZ0(3)(:LZ),
     &            STIM(0)(:LTIM), STR4, 'Step'
               WRITE (*, 10030, IOSTAT=IDUM) 'Minimum', SVAL(1)(:LVAL),
     &            SXYZ(1,1)(:LX), SXYZ(1,2)(:LY), SXYZ(1,3)(:LZ),
     &            STIM(1)(:LTIM), NUMMIN, ISTMIN
               WRITE (*, 10030, IOSTAT=IDUM) 'Maximum', SVAL(2)(:LVAL),
     &            SXYZ(2,1)(:LX), SXYZ(2,2)(:LY), SXYZ(2,3)(:LZ),
     &            STIM(2)(:LTIM), NUMMAX, ISTMAX
10020           FORMAT (4X, A7, 1X, A, 3X, 3 (1X, A),
     &            4X, A, 3X, A4, 3X, A4)
10030           FORMAT (4X, A7, 1X, A, 3X, 3 (1X, A), 4X, A, I7, I7)
            END IF
         ELSE
            WRITE (*, 10040) 'Range: ', SVAL(0)(:LVAL),
     &         STIM(0)(:LTIM), 'Step'
            WRITE (*, 10050, IOSTAT=IDUM) 'Minimum', SVAL(1)(:LVAL),
     &         STIM(1)(:LTIM), ISTMIN
            WRITE (*, 10050, IOSTAT=IDUM) 'Maximum', SVAL(2)(:LVAL),
     &         STIM(2)(:LTIM), ISTMAX
10040        FORMAT (4X, A7, 1X, A, 4X, A, 3X, A4)
10050        FORMAT (4X, A7, 1X, A, 4X, A, I7)
         END IF
      END IF

      RETURN
10060  FORMAT (' *** WARNING - ', 5A)
      END
