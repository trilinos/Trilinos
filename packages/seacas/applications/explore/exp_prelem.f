C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRELEM (OPTION, NOUT,
     &  NELBLK, NUMEL, NLISEL, LISEL, LENE,
     &  NVAREL, LISEV, NAMEEV, ISEVOK, VAREL, nvardm,
     *  MAP, DOMAP)
C=======================================================================

C   --*** PRELEM *** (BLOT) Display current database element variables
C   --
C   --PRELEM displays the element data for a time step.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NELBLK - IN - the number of element blocks
C   --   NUMEL - IN - the number of elements
C   --   NLISEL - IN - the number of selected elements by block
C   --   LISEL - IN - the indices of the selected elements by block
C   --   LENE - IN - the cumulative element counts by element block
C   --   NVAREL - IN - the number of element variables
C   --   LISEV - IN - the indices of the selected element variables
C   --   NAMEEV - IN - the names of the element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j) is NOT 0
C   --   VAREL - IN - the selected element variables for the time step

      include 'exodusII.inc'
      CHARACTER*(*) OPTION
      INTEGER NLISEL(0:NELBLK)
      INTEGER LISEL(0:*)
      INTEGER LENE(0:NELBLK)
      INTEGER LISEV(0:*)
      CHARACTER*(*) NAMEEV(*)
      INTEGER ISEVOK(nvardm,*)
      REAL VAREL(NUMEL, *)
      INTEGER MAP(*)
      LOGICAL DOMAP
      INTEGER GETPRC

      CHARACTER*32 CVAL(6)
      CHARACTER*40 FMT20,FMT30, FMT40
      INTEGER PRTLEN

      PRTLEN = GETPRC() + 7
      WRITE(FMT20,20) PRTLEN, PRTLEN-7
      WRITE(FMT30,30) PRTLEN
      WRITE(FMT40,40) PRTLEN

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      if (domap) then
        if (nout .gt. 0) then
          write (nout, 10005)
        else
          write (*, 10005)
        end if
      end if

      do 90 i=1, lisev(0)
        irow = ((i-1)/5)+1
        icol = i - (irow-1)*5
        IF (NOUT .GT. 0) THEN
           WRITE (NOUT, 10010) irow, icol, NAMEEV(LISEV(I))
        ELSE
           WRITE (*, 10010) irow, icol, NAMEEV(LISEV(I))
        END IF
 90   continue

      DO 130 IELB = 1, NELBLK
         IF (NLISEL(IELB) .GT. 0) THEN

            IX0 = LENE(IELB-1)
            DO 120 IX = 1, NLISEL(IELB)
               IEL = LISEL(IX0+IX)
               if (domap) then
                 id = map(iel)
               else
                 id = iel
               end if

               DO 110 IVAR = 1, LISEV(0), 5
                  MINVAL = IVAR
                  MAXVAL = MIN (LISEV(0), IVAR+5-1)
                  NVAL = MAXVAL - MINVAL + 1
                  DO 100 I = MINVAL, MAXVAL
                     IF (ISEVOK (LISEV(I),IELB) .NE. 0) THEN
                          WRITE (CVAL(I-MINVAL+1), FMT20, IOSTAT=IDUM)
     &                     VAREL(IEL, LISEV(I))
                     ELSE
                        CVAL(I-MINVAL+1) = '-----------'
                     END IF
  100             CONTINUE

                  IF (IVAR .LE. 1) THEN
                     IF (NOUT .GT. 0) THEN
                        WRITE (NOUT, FMT30, IOSTAT=IDUM)
     &                     ID, (CVAL(I), I=1,NVAL)
                     ELSE
                        WRITE (*, FMT30, IOSTAT=IDUM)
     &                     ID, (CVAL(I), I=1,NVAL)
                     END IF
                  ELSE
                     IF (NOUT .GT. 0) THEN
                        WRITE (NOUT, FMT40, IOSTAT=IDUM)
     &                     (CVAL(I), I=1,NVAL)
                     ELSE
                        WRITE (*, FMT40, IOSTAT=IDUM)
     &                     (CVAL(I), I=1,NVAL)
                     END IF
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
  130 CONTINUE

      RETURN

10000  FORMAT (/, 1X, 'ELEMENT TIME STEP VARIABLES')
10005  format (1X, 'Element Ids are Global')
10010  FORMAT (1X, 'Row ',I4,', Column ',I1,' is variable ',A)
 20   FORMAT('(1PE',I2.2,'.',I2.2,')')
 30   FORMAT('(1X, ''Element'', I12, 5(2X,A',I2,'))')
 40   FORMAT('(18X, 5 (2X, A',I2,'))')
      END
