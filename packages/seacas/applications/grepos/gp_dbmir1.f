C Copyright(C) 1999-2024 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBMIR1 (IELB, NUMELB, NUMLNK, LINK, TYPE, NDIM, NONQUD)
C=======================================================================

C   --*** DBMIR1 *** (GJOIN) Fixup element connectivity for reflections
C   --   Written by Greg Sjaardema - revised 02/10/89
C   --   Modified from DBOEB1 Written by Amy Gilkey
C   --
C   --Parameters:
C   --   IELB - IN - the element block number
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element
C   --   LINK - IN/OUT - the element connectivity for this block
C   --   TYPE - IN - the element type for this block
C   --   NDIM - IN - the spatial dimension (1,2,3)
C   --

      include 'exodusII.inc'

      INTEGER NUMELB, NUMLNK
      INTEGER LINK(NUMLNK,*)
      CHARACTER*(MXSTLN) TYPE
      LOGICAL NONQUD

      CHARACTER*132 STRING

      IF (NUMELB .GT. 0) THEN

C...8-node Hexes
        IF ((NUMLNK .EQ. 8) .AND. (NDIM .EQ. 3) .AND.
     *    TYPE(:3) .EQ. 'HEX') THEN
          DO NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(4,NE)
            LINK(4,NE) = ILTMP

            ILTMP = LINK (6,NE)
            LINK(6,NE) = LINK(8,NE)
            LINK(8,NE) = ILTMP
          END DO

C...20-node Hexes
        ELSE IF ((NUMLNK .EQ. 20) .AND. (NDIM .EQ. 3) .AND.
     *    TYPE(:3) .EQ. 'HEX') THEN
          DO NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(4,NE)
            LINK(4,NE) = ILTMP

            ILTMP = LINK (6,NE)
            LINK(6,NE) = LINK(8,NE)
            LINK(8,NE) = ILTMP

            ILTMP = LINK ( 9,NE)
            LINK( 9,NE) = LINK(12,NE)
            LINK(12,NE) = ILTMP

            ILTMP = LINK (10,NE)
            LINK(10,NE) = LINK(11,NE)
            LINK(11,NE) = ILTMP

            ILTMP = LINK (14,NE)
            LINK(14,NE) = LINK(16,NE)
            LINK(16,NE) = ILTMP

            ILTMP = LINK (17,NE)
            LINK(17,NE) = LINK(20,NE)
            LINK(20,NE) = ILTMP

            ILTMP = LINK (18,NE)
            LINK(18,NE) = LINK(19,NE)
            LINK(19,NE) = ILTMP
          END DO

C...Quads/Shells
        ELSE IF ((NUMLNK .EQ. 4) .AND.
     *      (TYPE(:4) .EQ. 'QUAD' .OR. TYPE(:5) .EQ. 'SHELL')) THEN
          if (type(:5) .eq. 'SHELL') NONQUD = .TRUE.
          DO NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(4,NE)
            LINK(4,NE) = ILTMP
          END DO

C...four-node tets...
        ELSE IF ((NUMLNK .EQ. 4) .AND.
     *      (TYPE(:3) .EQ. 'TET')) THEN
          NONQUD = .TRUE.
          DO NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(3,NE)
            LINK(3,NE) = ILTMP
          END DO

C...ten-node tets...
        ELSE IF ((NUMLNK .EQ. 10) .AND.
     *      (TYPE(:3) .EQ. 'TET')) THEN
          NONQUD = .TRUE.
          DO NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(3,NE)
            LINK(3,NE) = ILTMP

            ILTMP = LINK (5,NE)
            LINK(5,NE) = LINK(7,NE)
            LINK(7,NE) = ILTMP

            ILTMP = LINK (9,NE)
            LINK(9,NE) = LINK(10,NE)
            LINK(10,NE) = ILTMP
          END DO

C...Bars
        ELSE IF (NUMLNK .EQ. 2) THEN
          NONQUD = .TRUE.
          DO NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(1,NE)
            LINK(1,NE) = ILTMP
          END DO

C...Triangles
        ELSE IF (NUMLNK .EQ. 3) THEN
          NONQUD = .TRUE.
          DO NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(3,NE)
            LINK(3,NE) = ILTMP
          END DO

C...6-node Triangles
        ELSE IF (NUMLNK .EQ. 6 .and. type(:3) .eq. 'TRI') then
          NONQUD = .TRUE.
          DO NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(3,NE)
            LINK(3,NE) = ILTMP
            ILTMP = LINK (4,NE)
            LINK(4,NE) = LINK(6,NE)
            LINK(6,NE) = ILTMP
          END DO
        ELSE IF (NUMLNK .EQ. 6 .and. type(:5) .eq. 'WEDGE') then
          NONQUD = .TRUE.
          DO NE = 1, NUMELB
            ILTMP = LINK (2,NE)
            LINK(2,NE) = LINK(3,NE)
            LINK(3,NE) = ILTMP
            ILTMP = LINK (5,NE)
            LINK(5,NE) = LINK(6,NE)
            LINK(6,NE) = ILTMP
          END DO
        ELSE
          NONQUD = .TRUE.
          WRITE (STRING, 100) IELB, NUMLNK, TYPE
 100      FORMAT('Element block ',I5,' contains ',I2,'-node ',A,
     *      ' elements which are not supported for mirroring by grepos')
          CALL SQZSTR (STRING, LSTR)
          CALL PRTERR ('PROGRAM', STRING(:LSTR))
        END IF
      END IF
      RETURN
      END
