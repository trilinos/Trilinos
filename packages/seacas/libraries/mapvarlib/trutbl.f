C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK,TRUTBL
      SUBROUTINE TRUTBL(MP,IMP,IDA,IDB,ITRTA,ITRTB)

C     ******************************************************************

C     SUBROUTINE TO BUILD THE TRUTH TABLE FOR THE RECIPIENT MESH
C     THIS IMPROVES (SPEEDS-UP) THE PROCESS OF WRITING ELEMENT
C     VARIABLES

C     Called by MAPVAR

C     ******************************************************************

C  MP         INT   Array of element block mapping from donor
C                   mesh to recipient mesh (1:2,1:IMP)
C  IMP        INT   Number of entries in MP
C  IDA        INT   Donor mesh element block I.D. (1:nblksa)
C  IDB        INT   Recipient mesh element block I.D. (1:nblksa)
C  ITRTA      INT   Donor mesh truth table (1:nvarel,1:nblksa)
C  ITRTB      INT   Recipient mesh truth table (1:nvarel,1:nblksb)

C     ******************************************************************

      include 'aexds1.blk'
      include 'amesh.blk'
      include 'bmesh.blk'
      include 'ex2tp.blk'
      include 'tapes.blk'

      DIMENSION MP(3,*),IDA(*),IDB(*),ITRTA(NVAREL,*),ITRTB(NVAREL,*)

C     ******************************************************************

C ... For later sanity check
      do 2 iblkb = 1, nblksb
        do 1 ivar = 1, nvarel
          itrtb(ivar,iblkb) = -999
 1      continue
 2    continue

      CALL EXGVTT (NTP2EX,NBLKSA,NVAREL,ITRTA,IERR)

C Match element block I.D.'s from donor and recipient mesh

      DO 10 IBLKB = 1, NBLKSB
        IDBLKB = IDB(IBLKB)
        DO 20 IBLKA = 1, NBLKSA
          IDBLKA = IDA(IBLKA)
          DO 30 IM = 1, IMP
            IF (IDBLKB .EQ. MP(2,IM) .AND. IDBLKA .EQ. MP(1,IM))THEN
              DO 40 IVAR = 1, NVAREL
                ITRTB(IVAR,IBLKB) = ITRTA(IVAR,IBLKA)
   40         CONTINUE
              GO TO 20
            END IF
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE

C Do some error checking - write a warning

      DO 50 IBLKB = 1, NBLKSB
        IFND = 0
        IDBLKB = IDB(IBLKB)
        DO 60 IM = 1, IMP
          IF (IDBLKB .EQ. MP(2,IM))THEN
            IFND = 1
          END IF
   60   CONTINUE
        IF (IFND .EQ. 0)THEN
          DO 70 IVAR = 1, NVAREL
            ITRTB(IVAR,IBLKB) = 0
   70     CONTINUE
          WRITE(NOUT,1000)IBLKB,IDBLKB
          WRITE(NTPOUT,1000)IBLKB,IDBLKB
        END IF
   50 CONTINUE

C ... Sanity check
      do 90 iblkb = 1, nblksb
        do 80 ivar = 1, nvarel
          if (itrtb(ivar,iblkb) .eq. -999) then
            idblkb = idb(iblkb)
            write(nout, 1010) iblkb, idblkb, ivar
            write(ntpout, 1010) iblkb, idblkb, ivar
          end if
 80     continue
 90   continue

      CALL EXPVTT(NTP4EX,NBLKSB,NVAREL,ITRTB,IERR)

 1000 FORMAT(5X,'MAPPING BACK INTO DONOR MESH FOR',/,
     1'          RECIPIENT MESH ELEMENT BLOCK NUMBER',I7,/,
     2'          ELEMENT BLOCK I. D.',I7,/,
     3'          WAS NOT FOUND. THIS ELEMENT BLOCK WILL',/,
     4'          NOT BE MAPPED.')
 1010 FORMAT(5X,'Block Number ', i7, /,
     1  'Element Block ID ', i7,/,
     2  'Has unset truth table entry for variable ',i7)
      RETURN
      END
