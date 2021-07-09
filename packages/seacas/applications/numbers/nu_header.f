C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE HEADER (NDIM, TITLE, NUMEL, NUMNP, AXI, GENFIL)

      include 'nu_io.blk'
      CHARACTER*16 FORM(3)
      CHARACTER*80 TITLE
      CHARACTER*(*) GENFIL
      LOGICAL AXI, FIRST
      CHARACTER*6 DIMEN(3)

      DATA DIMEN/'One',   'Two',   'Three'/
      DATA FORM /', Planar',', Axisymmetric',' '/
      DATA FIRST /.TRUE./

      IO = IHARD
      IF (FIRST) THEN
          WRITE (IO, 10) GENFIL(:LENSTR(GENFIL))
   10     FORMAT (5X,'Genesis: ',A/)
          WRITE (IO, 20) TITLE(:LENSTR(TITLE))
   20     FORMAT (5X,'Title:   ',A/)
          ILAB = 2
          IF (NDIM .EQ. 2 .AND. .NOT. AXI) ILAB = 1
          IF (NDIM .EQ. 3) ILAB = 3
          WRITE (IO, 30) NUMNP,
     *        NUMEL, DIMEN(NDIM)(:LENSTR(DIMEN(NDIM))),
     *        FORM(ILAB)(:LENSTR(FORM(ILAB)))
   30     FORMAT (5X,'Number of Nodes:    ',I10,/
     *        5X,'Number of Elements: ',I10/
     *        5X,A,'-Dimensional Mesh',A)
          FIRST = .FALSE.
      END IF
      RETURN
      END
