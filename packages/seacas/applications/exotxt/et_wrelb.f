C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRELB (NTXT, NELBLK, IDELB, NUMELB, NUMLNK, NUMATR,
     &                  NAMELB, LINK, ATRIB)
C=======================================================================
c      SUBROUTINE WRELB (NTXT, IELB, IDELB, NUMELB, NUMLNK, NUMATR,
c     &   LINK, ATRIB, natrdm)

C   --*** WRELB *** (EXOTXT) Write database element blocks
C   --   Written by Amy Gilkey - revised 02/27/86
C   --   Modified for ExodusIIv2 - 10/12/95
C   --
C   --WRELB writes the element block information from the database,
C   --including the element block connectivity and attribute information.
C   --
C   --Parameters:
C   --   NTXT   - IN - the text file
C   --   NELBLK - IN - number of element blocks
C   --   IDELB  - IN - the element block ID
C   --   NUMELB - IN - array: number of elements in each block
C   --   NUMLNK - IN - array: number of nodes per element in each block
C   --   NUMATR - IN - array: number of attributes in each  block
C   --   NAMELB - IN - array: type of element in each block
C   --   LINK   - IN - the element connectivity for this block
C   --   ATRIB  - IN - the attributes for this block

      include 'exodusII.inc'

      INTEGER NTXT, NELBLK
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      CHARACTER*(MXSTLN) NAMELB(*)
      INTEGER LINK(*)
      REAL    ATRIB(*)
      INTEGER IDXCON
      INTEGER IDXATT

      IDXCON = 1
      IDXATT = 1
      NEL  = 0
      NLNK = 0
      NATR = 0
      DO 10 I = 1, NELBLK
         WRITE (NTXT, '(A, I5)') '! Element block', I
         WRITE (NTXT, 10030)IDELB(I),NUMELB(I),
     &        NAMELB(I)(:LENSTR(NAMELB(I))),
     &        '! ID, elements, element type'
         WRITE (NTXT, 10020)NUMLNK(I),NUMATR(I),
     &   '! nodes per element, attributes'

         IDXCON = IDXCON + NEL*NLNK
         IDXATT = IDXATT + NEL*NATR
C        Write connectivity and attributes for this element block
         CALL WRCONAT(NTXT, NUMELB(I), NUMLNK(I), NUMATR(I),
     &                LINK(IDXCON), ATRIB(IDXATT))
         NEL  = NUMELB(I)
         NLNK = NUMLNK(I)
         NATR = NUMATR(I)

   10 CONTINUE

10020 FORMAT (2I10, 6X, A)
10030 FORMAT (2I10, 6X, A, 6X, A)

      RETURN
      END

      SUBROUTINE WRCONAT(NTXT, NUMELB, NUMLNK, NUMATR, LINK, ATRIB)

C     WRCONAT - Write the connectivity and attributes to a text file
C     Written for ExodusIIv2 database format 10/12/95

C     NTXT   - IN - file id
C     NUMELB - IN - number of element in element block
C     NUMLNK - IN - number of nodes per element in element block
C     NUMATR - IN - number of attributes in element block
C     LINK   - IN - connectivity array
C     ATRIB  - IN - attributes array

      INTEGER NTXT, NUMELB, NUMLNK, NUMATR
      INTEGER LINK(NUMLNK,NUMELB)
      REAL    ATRIB(NUMATR,NUMELB)

C     Write attributes

      WRITE (NTXT, '(A)') '! Connectivity'
      DO 100 NE = 1, NUMELB
         WRITE (NTXT, 10010) (LINK(I,NE), I=1,NUMLNK)
  100 CONTINUE

      IF (NUMATR .GT. 0) THEN
         WRITE (NTXT, '(A)') '! Attributes'
         DO 110 NE = 1, NUMELB
            WRITE (NTXT, 10000) (ATRIB(I,NE), I=1,NUMATR)
  110    CONTINUE
      END IF
10000 FORMAT (5(1pE16.7))
10010 FORMAT (8I10)

      RETURN
      END
