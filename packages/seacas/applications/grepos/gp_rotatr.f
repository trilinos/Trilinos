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
      SUBROUTINE ROTATR (NELBLK, NDIM, IDELB, BLKTYP, NUMATR,
     *  NUMELB, ATRIB)
C=======================================================================

      include 'exodusII.inc'
      include 'gp_attrot.blk'

      CHARACTER*(MXSTLN) BLKTYP(*)
      INTEGER NUMATR(*)
      INTEGER NUMELB(*)
      INTEGER IDELB(*)
      REAL ATRIB(*)

      IEATR = 0
      IAT   = 1

      DO IELB = 1, NELBLK
        ISATR = IEATR + 1
        IEATR = IEATR + NUMATR(IELB) * NUMELB(IELB)
        if (numatr(ielb) .ge. ATTIND+NDIM-1) THEN
          if (rotall .or.
     *      (.NOT. ROTTYP .AND. IDELB(IELB) .EQ. ATTBLK)  .OR.
     *      (ROTTYP .AND. BLKTYP(IELB) .EQ. ROTBLK)) then
            CALL ROTAT1 (NDIM, NUMELB(IELB), NUMATR(IELB), ATTIND,
     *        ATRIB(ISATR), ROTATT)
          end if
        end if
        IAT = IAT + NUMATR(IELB)
      END DO

      RETURN
      END

      SUBROUTINE ROTAT1(NDIM, NUMEL, NUMATR, IDX, ATRIB, ROTMAT)
      REAL ATRIB(*)
      REAL ROTMAT(3,3)

      IBEG = 0
      IF (NDIM .EQ. 3) THEN
        DO IEL = 1, NUMEL
          X = ATRIB(IDX+0 + IBEG)
          Y = ATRIB(IDX+1 + IBEG)
          Z = ATRIB(IDX+2 + IBEG)
          XN = X*ROTMAT(1,1) + Y*ROTMAT(2,1) + Z*ROTMAT(3,1)
          YN = X*ROTMAT(1,2) + Y*ROTMAT(2,2) + Z*ROTMAT(3,2)
          ZN = X*ROTMAT(1,3) + Y*ROTMAT(2,3) + Z*ROTMAT(3,3)
          ATRIB(IDX+0 + IBEG) = XN
          ATRIB(IDX+1 + IBEG) = YN
          ATRIB(IDX+2 + IBEG) = ZN
          IBEG = IBEG + NUMATR
        END DO
      ELSE IF (NDIM .EQ. 2) THEN
        DO IEL = 1, NUMEL
          X = ATRIB(IDX+0 + IBEG)
          Y = ATRIB(IDX+1 + IBEG)
          XN = X*ROTMAT(1,1) + Y*ROTMAT(2,1)
          YN = X*ROTMAT(1,2) + Y*ROTMAT(2,2)
          ATRIB(IDX+0 + IBEG) = XN
          ATRIB(IDX+1 + IBEG) = YN
          IBEG = IBEG + NUMATR
        END DO
      END IF

      RETURN
      END

