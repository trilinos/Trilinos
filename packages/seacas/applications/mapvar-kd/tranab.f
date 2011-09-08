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

      SUBROUTINE TRANAB(IELPT,SOLEA,SOLEB,
     &                  IDBLKA,IDBLKB,
     &                  ITT,iblk,TIMES,CENTER,
     &                  INSUB,ICOMPL,
     &                  XB,YB,ZB,ICONB,DUME)
C
      include 'exodusII.inc'
C
      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'bmesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'steps.blk'
      include 'tapes.blk'
C
      DIMENSION SOLEB(NUMEBB,*),SOLEA(NUMEBA,*),IELPT(*),
     &          ITT(NVAREL,*),TIMES(*),CENTER(NUMEBB,*)
      DIMENSION XB(*),YB(*),ZB(*),XX(27),YY(27),ZZ(27),ICONB(NELNDB,*)
      DIMENSION DUME(*)
C
C
      IF (ISTEP .EQ. -1)THEN
        NTM = NTIMES
      ELSE
        NTM = 1
      END IF
C
      DO 5 IST = 1, NTM
        IF (ISTEP .EQ. -1)THEN
          ISTP = IST
        ELSE
          ISTP = ISTEP
        END IF
C
        DO 10 IVAR = 1, NVAREL
          IF (ITT(IVAR,iblk) .EQ. 0)GO TO 10
C
C If first time into subroutine for this element block,
C initialize the SOLEB array
C If not first time into subroutine for this element block
C retrieve SOLEB from storage in EXODUS
C
          IF (INSUB .EQ. 1)THEN
            CALL INIELT(SOLEB,IVAR,TIMES,ISTP,IDBLKB,CENTER,DUME)
          ELSE
            CALL EXGEV(NTP4EX,IST,IVAR,IDBLKB,NUMEBB,SOLEB(1,IVAR),
     &                 IERR)
          END IF
C
          CALL EXGEV(NTP2EX,ISTP,IVAR,IDBLKA,NUMEBA,SOLEA(1,IVAR),IERR)
          DO 20 IELB = 1, NUMEBB
            IELA = IELPT(IELB)
            IF (IELA .NE. 0)THEN
              SOLEB(IELB,IVAR) = SOLEA(IELA,IVAR)
            END IF
   20     CONTINUE
C
C If there is more searching to do (i.e. many blocks to one)
C use EXODUS as temporary storage
C don't bother to perform needed adjustments yet
C
          IF (ICOMPL .NE. 1)THEN
            CALL EXPEV(NTP4EX,IST,IVAR,IDBLKB,NUMEBB,SOLEB(1,IVAR),
     &                 IERR)
          ELSE
C
C Make needed adjustments to element varaible data and
C write element vars out to EXODUS data base
C
C ELMASS is special
C
            IF (NAMVAR(nvargp+IVAR)(1:6) .EQ. 'ELMASS')THEN
C
C ELMASS was changed to nodal density prior to processing.
C need to go back from density to element mass now
C
C              NNODES=NNELM(ITYPE)
              NNODES = NELNDB
              IF (ITYPE .EQ. 6)NNODES = 4
              DO 100 IEL = 1, NUMEBB
                DO 105 I = 1, NNODES
                  XX(I) = XB(ICONB(I,IEL))
                  YY(I) = YB(ICONB(I,IEL))
                  IF (NDIMB .EQ. 3)THEN
                    ZZ(I) = ZB(ICONB(I,IEL))
                  ELSE
                    ZZ(I) = 0.
                  END IF
 105            CONTINUE
                CALL VOL(ITYPE,XX,YY,ZZ,VOLUME)
                SOLEB(IEL,IVAR) = SOLEB(IEL,IVAR) * VOLUME
 100          CONTINUE
            END IF
            CALL EXPEV(NTP4EX,IST,IVAR,IDBLKB,NUMEBB,SOLEB(1,IVAR),
     &                 IERR)
          END IF
   10   CONTINUE
    5 CONTINUE
      RETURN
      END
