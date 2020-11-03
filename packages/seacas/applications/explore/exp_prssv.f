C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRSSV (NOUT, NSTEP, NUMESS, LISESS, LESSEL,
     &     IDESS, NEESS, IXEESS, LTEESS, LTSESS, NAME,
     $     NVAR,  LISVAR, NAMEV, ISVOK, VARS, nvardm,
     $     MAPEL, DOMAP)
C=======================================================================

C     --*** PRSSV *** (EXPLORE) Display database nodal point set
C     --
C     --PRSSV displays the sides set vars
C     --
C     --Parameters:
C     --   NOUT - IN - the output file, <=0 for standard
C     --   NUMESS - IN - the number of nodal point sets
C     --   LISESS - IN - the indices of the selected side sets
C     --   LESSEL - IN - the number of elements for all sets
C     --   IDESS - IN - the set ID for each set
C     --   NEESS - IN - the number of elements for each set
C     --   IXEESS - IN - the index of the first element for each set
C     --   LTEESS - IN - the elements for all sets
C     --   LTSESS - IN - the element sides for all sets
      include 'exodusII.inc'
      include 'exp_dbase.blk'
      INTEGER LISESS(0:*)
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER IXEESS(*)
      INTEGER LTEESS(*)
      INTEGER LTSESS(*)
      CHARACTER*(*) NAME(*)

      INTEGER LISVAR(0:*)
      CHARACTER*(*) NAMEV(*)
      INTEGER ISVOK(nvardm,*)
      REAL VARS(LESSEL, *)

      INTEGER MAPEL(*)
      LOGICAL DOMAP

      CHARACTER*32 CVAL(6)
      CHARACTER*40 FMT20,FMT30, FMT40
      INTEGER PRTLEN
      INTEGER GETPRC

      CHARACTER*20 STRA, STRB

C ... See if need to read the data
      if (nstep .ne. nstepss) then
         nstepss = nstep
         DO 20 IX = 1, LISESS(0)
            IESS = LISESS(IX)
            IS = IXEESS(IESS)
            IE = IS + NEESS(IESS) - 1

            ID   = idess(iess)
            DO 10 IVAR = 1, LISVAR(0)
               IF (ISVOK (LISVAR(IVAR),IESS) .NE. 0) THEN
                  call exgssv(ndb, nstep, lisvar(ivar), id, neess(iess),
     $                 vars(is,lisvar(ivar)), ierr)
               end if
 10         continue
 20      continue
      end if

      PRTLEN = GETPRC() + 7
      WRITE(FMT20,2000) PRTLEN, PRTLEN-7
      WRITE(FMT30,3000) PRTLEN
      WRITE(FMT40,4000) PRTLEN

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      if (domap) then
         if (nout .gt. 0) then
            write (nout, 10005)
         else
            write (*, 10005)
         end if
      end if

      do 90 i=1, lisvar(0)
        irow = ((i-1)/5)+1
        icol = i - (irow-1)*5
        IF (NOUT .GT. 0) THEN
           WRITE (NOUT, 10010) irow, icol, NAMEV(LISVAR(I))
        ELSE
           WRITE (*, 10010) irow, icol, NAMEV(LISVAR(I))
        END IF
 90   continue

      WRITE (STRA, 10001, IOSTAT=IDUM) NUMESS
10001 FORMAT ('(#', I4, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)
      WRITE (STRB, 10002, IOSTAT=IDUM) LESSEL
10002 FORMAT ('(index=', I12, ')')
      CALL PCKSTR (1, STRB)
      LSTRB = LENSTR (STRB)

      DO 200 IX = 1, LISESS(0)
         IESS = LISESS(IX)
         WRITE (STRA, 10001, IOSTAT=IDUM) IESS
         CALL PCKSTR (1, STRA)
         WRITE (STRB, 10002, IOSTAT=IDUM) IXEESS(IESS)
         CALL PCKSTR (1, STRB)
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10030, IOSTAT=IDUM)
     &           IDESS(IESS), STRA(:LSTRA),
     &           NEESS(IESS), STRB(:LSTRB),
     $           NAME(IESS)(:LENSTR(NAME(IESS)))
         ELSE
            WRITE (*, 10030, IOSTAT=IDUM)
     &           IDESS(IESS), STRA(:LSTRA),
     &           NEESS(IESS), STRB(:LSTRB),
     $           NAME(IESS)(:LENSTR(NAME(IESS)))
         END IF

         IF (NEESS(IESS) .GT. 0) THEN
            IS = IXEESS(IESS)
            IE = IS + NEESS(IESS) - 1
            do 130 IN = IS,IE
               if (domap) then
                  ID = MAPEL(LTEESS(IN))
               else
                  ID = LTEESS(IN)
               end if

               DO 110 IVAR = 1, LISVAR(0), 5
                  MINVAL = IVAR
                  MAXVAL = MIN (LISVAR(0), IVAR+5-1)
                  NVAL = MAXVAL - MINVAL + 1
                  DO 100 I = MINVAL, MAXVAL
                     IF (ISVOK (LISVAR(I),IESS) .NE. 0) THEN
                        WRITE (CVAL(I-MINVAL+1), FMT20, IOSTAT=IDUM)
     &                       VARS(IN, LISVAR(I))
                     ELSE
                        CVAL(I-MINVAL+1) = '-----------'
                     END IF
 100              CONTINUE

                  IF (IVAR .LE. 1) THEN
                     IF (NOUT .GT. 0) THEN
                        WRITE (NOUT, FMT30, IOSTAT=IDUM)
     &                       ID, LTSESS(IN), (CVAL(I), I=1,NVAL)
                     ELSE
                        WRITE (*, FMT30, IOSTAT=IDUM)
     &                       ID, LTSESS(IN), (CVAL(I), I=1,NVAL)
                     END IF
                  ELSE
                     IF (NOUT .GT. 0) THEN
                        WRITE (NOUT, FMT40, IOSTAT=IDUM)
     &                       (CVAL(I), I=1,NVAL)
                     ELSE
                        WRITE (*, FMT40, IOSTAT=IDUM)
     &                       (CVAL(I), I=1,NVAL)
                     END IF
                  END IF
 110           CONTINUE
 130        CONTINUE
         end if
 200  CONTINUE

      RETURN

10000 FORMAT (/, 1X, 'SIDESET TIME STEP VARIABLES')
10010 FORMAT (1X, 'Row ',I4,', Column ',I1,' is variable ',A)
 2000 FORMAT('(1PE',I2.2,'.',I2.2,')')
 3000 FORMAT('(1X, ''Element'', I12,''.'',I1, 5(2X,A',I2,'))')
 4000 FORMAT('(20X, 5 (2X, A',I2,'))')

10005 FORMAT (1X, 'Element ids are Global')
10030 FORMAT (/,1X, 'Set', I12, 1X, A, ':',
     &     I12, ' elements', 1X, A, ' name = "',A,'"')
      END
