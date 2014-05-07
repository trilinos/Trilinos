C $Id: header.f,v 1.2 1998/03/22 05:34:35 gdsjaar Exp $
C $Log: header.f,v $
C Revision 1.2  1998/03/22 05:34:35  gdsjaar
C General cleanp of unused variables. Reordered DATA statements in
C command.f so would compile with f2c.
C
C Revision 1.1.1.1  1991/02/21 15:43:37  gdsjaar
C NUMBERS: Greg Sjaardema, initial Unix release
C
c Revision 1.1  1991/02/21  15:43:36  gdsjaar
c Initial revision
c
      SUBROUTINE HEADER (NDIM, TITLE, NUMEL, NUMNP, AXI)
C
      include 'nu_io.blk'
      CHARACTER*16 FORM(3)
      CHARACTER*80 TITLE
      CHARACTER*256 GENFIL
      LOGICAL AXI, FIRST
      CHARACTER*6 DIMEN(3)

      DATA DIMEN/'One',   'Two',   'Three'/
      DATA FORM /', Planar',', Axisymmetric',' '/
      DATA FIRST /.TRUE./
C
      INQUIRE (NDB,NAME=GENFIL)
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
   30     FORMAT (5X,'Number of Nodes:    ',I6,/
     *        5X,'Number of Elements: ',I6/
     *        5X,A,'-Dimensional Mesh',A)
          FIRST = .FALSE.
      END IF
      RETURN
      END
