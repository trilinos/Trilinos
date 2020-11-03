C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE TRANAB(IELPT,SOLEA,SOLEB,
     &                  IDBLKA,IDBLKB,
     &                  ITT,iblk,TIMES,CENTER,
     &                  INSUB,ICOMPL,
     &                  XB,YB,ZB,ICONB,DUME)

      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'bmesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'steps.blk'
      include 'tapes.blk'

      DIMENSION SOLEB(NUMEBB,*),SOLEA(NUMEBA,*),IELPT(*),
     &          ITT(NVAREL,*),TIMES(*),CENTER(NUMEBB,*)
      DIMENSION XB(*),YB(*),ZB(*),XX(27),YY(27),ZZ(27),ICONB(NELNDB,*)
      DIMENSION DUME(*)

      IF (ISTEP .EQ. -1)THEN
        NTM = NTIMES
      ELSE
        NTM = 1
      END IF

      DO 5 IST = 1, NTM
        IF (ISTEP .EQ. -1)THEN
          ISTP = IST
        ELSE
          ISTP = ISTEP
        END IF

        DO 10 IVAR = 1, NVAREL
          IF (ITT(IVAR,iblk) .EQ. 0)GO TO 10

C If first time into subroutine for this element block,
C initialize the SOLEB array
C If not first time into subroutine for this element block
C retrieve SOLEB from storage in EXODUS

          IF (INSUB .EQ. 1)THEN
            CALL INIELT(SOLEB,IVAR,TIMES,ISTP,IDBLKB,CENTER,DUME)
          ELSE
            CALL EXGEV(NTP4EX,IST,IVAR,IDBLKB,NUMEBB,SOLEB(1,IVAR),
     &                 IERR)
          END IF

          CALL EXGEV(NTP2EX,ISTP,IVAR,IDBLKA,NUMEBA,SOLEA(1,IVAR),IERR)
          DO 20 IELB = 1, NUMEBB
            IELA = IELPT(IELB)
            IF (IELA .NE. 0)THEN
              SOLEB(IELB,IVAR) = SOLEA(IELA,IVAR)
            END IF
   20     CONTINUE

C If there is more searching to do (i.e. many blocks to one)
C use EXODUS as temporary storage
C don't bother to perform needed adjustments yet

          IF (ICOMPL .NE. 1)THEN
            CALL EXPEV(NTP4EX,IST,IVAR,IDBLKB,NUMEBB,SOLEB(1,IVAR),
     &                 IERR)
          ELSE

C Make needed adjustments to element variable data and
C write element vars out to EXODUS data base

C ELMASS is special

            IF (NAMVAR(nvargp+IVAR)(1:6) .EQ. 'ELMASS')THEN

C ELMASS was changed to nodal density prior to processing.
C need to go back from density to element mass now

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
