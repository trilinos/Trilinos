C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PREB1 (OPTION, NOUT, IEL0, NLISEL, LISEL,
     &   NLINK, NATR, LINK, ATRIB, MAPEL, MAPND)
C=======================================================================

C   --*** PREB1 *** (BLOT) Display database element block information
C   --   Written by Amy Gilkey - revised 01/05/88
C   --
C   --PREB1 displays the element block internal information (either the
C   --connectivity or the attributes.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      'C' to print connectivity
C   --      'A' to print attributes
C   --   NOUT - IN - the output file, <=0 for standard
C   --   IEL0 - IN - the element offset of the elements in the block
C   --   NLISEL - IN - the number of selected elements
C   --   LISEL - IN - the indices of the selected elements
C   --   NLINK - IN - the number of nodes per element
C   --   NATR - IN - the number of attributes per element
C   --   LINK - IN - the connectivity array for this block
C   --   ATRIB - IN - the attribute array for this block

      CHARACTER*(*) OPTION
      INTEGER LISEL(*)
      INTEGER LINK(NLINK,*)
      REAL ATRIB(NATR,*)
      INTEGER MAPEL(*), MAPND(*)
      LOGICAL ISABRT
      LOGICAL DOCONN, DOATR

      IF (NLISEL .LE. 0) RETURN

      DOCONN = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0))
     &   .AND. (NLINK .GT. 0)
      DOATR  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0))
     &   .AND. (NATR .GT. 0)
      IF (.NOT. (DOCONN .OR. DOATR)) RETURN

      IF (NOUT .GT. 0) THEN
         IF (DOCONN .AND. DOATR) THEN
            WRITE (NOUT, 10000) 'Connectivity and Attributes'
         ELSE IF (DOCONN) THEN
            WRITE (NOUT, 10000) 'Connectivity'
         ELSE IF (DOATR) THEN
            WRITE (NOUT, 10000) 'Attributes'
         END IF
      ELSE
         IF (DOCONN .AND. DOATR) THEN
            WRITE (*, 10000) 'Connectivity and Attributes'
         ELSE IF (DOCONN) THEN
            WRITE (*, 10000) 'Connectivity'
         ELSE IF (DOATR) THEN
            WRITE (*, 10000) 'Attributes'
         END IF
      END IF

      DO 100 IX = 1, NLISEL
         IEL = LISEL(IX)
         NE = IEL - IEL0
         IF (ISABRT ()) RETURN
         IF (DOCONN) THEN
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10010, IOSTAT=IDUM)
     &          NE, MAPEL(IEL), (MAPND(LINK(I,NE)), I=1,NLINK)
            ELSE
               WRITE (*, 10010, IOSTAT=IDUM)
     &            NE, MAPEL(IEL), (MAPND(LINK(I,NE)), I=1,NLINK)
            END IF
         END IF
         IF (DOATR) THEN
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10020, IOSTAT=IDUM)
     &            NE, MAPEL(IEL), (ATRIB(I,NE), I=1,NATR)
            ELSE
               WRITE (*, 10020, IOSTAT=IDUM)
     &            NE, MAPEL(IEL), (ATRIB(I,NE), I=1,NATR)
            END IF
         END IF
  100 CONTINUE

      RETURN

10000  FORMAT (1X, '   Local ID  Global ID       ', A)
10010  FORMAT (1X, I11, I11, 5X, 8I11, :, /,
     &   (18X, 8I11))
10020  FORMAT (1X, I11, I11, 3X, 4 (2X, 1PE11.4), :, /,
     &   (16X, 4 (2X, 1PE11.4)))
      END
