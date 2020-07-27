C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PREBN (OPTION, NOUT, IELB, NUME, IDELB,
     *  IEL0, NLISEL, LISEL,
     *  NATR, LINK, ATRIB, natrdm,
     *  MAPN, DOMAPN, MAPE, DOMAPE)
C=======================================================================

C     --*** PREB1 *** (BLOT) Display database element block information
C     --
C     --PREB1 displays the element block internal information (either the
C     --connectivity or the attributes.
C     --
C     --Parameters:
C     --   OPTION - IN - '*' to print all, else print options:
C     --      'C' to print connectivity
C     --      'A' to print attributes
C     --   NOUT - IN - the output file, <=0 for standard
C     --   IEL0 - IN - the element offset of the elements in the block
C     --   NLISEL - IN - the number of selected elements
C     --   LISEL - IN - the indices of the selected elements
C     --   NATR - IN - the number of attributes per element
C     --   LINK - IN - the connectivity array for this block
C     --   ATRIB - IN - the attribute array for this block

      include 'exodusII.inc'
      include 'exp_dbase.blk'
      CHARACTER*(*) OPTION
      INTEGER LISEL(*)
      INTEGER LINK(*)
      REAL ATRIB(natrdm,*)
      INTEGER MAPN(*), MAPE(*)
      LOGICAL DOMAPN, DOMAPE

C ... Automatic Array...
      INTEGER NNPE(0:NUME)

      LOGICAL ALLSAM
      LOGICAL DOCONN, DOATR

      IF (NLISEL .LE. 0) RETURN

      DOCONN = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0))
      DOATR  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0))
     &  .AND. (NATR .GT. 0)
      IF (.NOT. (DOCONN .OR. DOATR)) RETURN

      IF (NOUT .GT. 0) THEN
        IF (DOCONN .AND. DOATR) THEN
          WRITE (NOUT, 10000) 'Connectivity and Attributes'
        ELSE IF (DOCONN) THEN
          WRITE (NOUT, 10000) 'Connectivity'
        ELSE IF (DOATR) THEN
          WRITE (NOUT, 10000) 'Attributes'
        END IF
        if (domape) then
          write (nout, 10030)
        end if
      ELSE
        IF (DOCONN .AND. DOATR) THEN
          WRITE (*, 10000) 'Connectivity and Attributes'
        ELSE IF (DOCONN) THEN
          WRITE (*, 10000) 'Connectivity'
        ELSE IF (DOATR) THEN
          WRITE (*, 10000) 'Attributes'
        END IF
        if (domape) then
          write (*, 10030)
        end if
      END IF

      if (doconn) then
C ... Need to read the nodes/element array for this element block...
        CALL EXGECPP(NDB, EXEBLK, IDELB, NNPE(1), IERR)
        nnpe(0) = 0
        do i=1, nume
          nnpe(i) = nnpe(i) + nnpe(i-1)
        end do

        do 100 ix=1, nlisel
          IEL = LISEL(IX)
          NE = IEL - IEL0
          NLINKS = NNPE(NE-1)
          NLINK  = NNPE(NE) - NLINKS
          if (domape) then
            ID = MAPE(IEL)
          else
            id = iel
          end if
          IF (NOUT .GT. 0) THEN
            if (domapn) then
              WRITE (NOUT, 10010, IOSTAT=IDUM)
     &          NE, ID, (MAPN(LINK(NLINKS+I)), I=1,NLINK)
            else
              WRITE (NOUT, 10010, IOSTAT=IDUM)
     &          NE, ID, (LINK(NLINKS+I), I=1,NLINK)
            end if
          ELSE
            if (domapn) then
              WRITE (*, 10010, IOSTAT=IDUM)
     &          NE, ID, (MAPN(LINK(NLINKS+I)), I=1,NLINK)
            else
              WRITE (*, 10010, IOSTAT=IDUM)
     &          NE, ID, (LINK(NLINKS+I), I=1,NLINK)
            end if
          END IF
 100    continue
      END IF

      if (doatr) then
C     ... See if all attributes are the same
        allsam = .true.
        do 180 ia = 1, natr
          iel = lisel(1)-iel0
          aval = atrib(ia, iel)
          do 170 ix = 2, nlisel
            iel = lisel(ix) - iel0
            if (aval .ne. atrib(ia, iel)) then
              allsam = .false.
              go to 190
            end if
 170      continue
 180    continue

C     ... Print either all values (if not all same), or the common values
 190    continue
        if (allsam) then
          ne = lisel(1) - iel0
          IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10025, IOSTAT=IDUM)
     &        (ATRIB(I,NE), I=1,NATR)
          ELSE
            WRITE (*, 10025, IOSTAT=IDUM)
     &        (ATRIB(I,NE), I=1,NATR)
          END IF
        else
          DO 200 IX = 1, NLISEL
            IEL = LISEL(IX)
            NE = IEL - IEL0
            IF (NOUT .GT. 0) THEN
              WRITE (NOUT, 10020, IOSTAT=IDUM)
     &          NE, IEL, (ATRIB(I,NE), I=1,NATR)
            ELSE
              WRITE (*, 10020, IOSTAT=IDUM)
     &          NE, IEL, (ATRIB(I,NE), I=1,NATR)
            END IF
 200      CONTINUE
        end if
      end if

      RETURN

10000 FORMAT (1X, '         #      elem            ', A)
10010 FORMAT (1X, I12, I12, 5X, 8I12, :, /,
     &  (26X, 8I12))
10020 FORMAT (1X, I12, I12, 3X, 4 (1X, 1pE11.4), :, /,
     &  (24X, 4 (1X, 1pE11.4)))
10025 FORMAT (1X, 'All attributes are equal. Values are: ',
     &  (3 (1x, 1pE11.4)), :, /,
     &  (3X, 6 (1X, 1pE11.4)))
10030 FORMAT (1x, 'Ids are global')
      END
