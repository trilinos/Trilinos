C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WCONAT (NDBOUT, NUMELB, NUMEBO, IEL, IXELEM,
     &   NUMLNK, NUMATR, NLNKDM, NATRDM, LINK, ATRIB,
     *   BLKTYP, IDELB, NEWNOD, NODIX, LNSCR, ATRSCR)
C=======================================================================

C   --*** WCONAT *** Write the element block connectivity and
C   --               attribute information to the output database file
C   --   WCONAT modified for EXODUSIIV2 8/29/95
C   --   This file was basically just renamed to try and indicate
C   --   What it actually does
C   --*** WREB1 *** (ALGEBRA) Write database element block misc.
C   --   Written by Amy Gilkey - revised 04/19/88
C   --
C   --Parameters:
C   --   NDBOUT   - IN - the output database file
C   --   NUMELB   - IN - the number of elements in the block
C   --   NUMEBO   - IN - the number of elements in the output block
C   --   IEL      - IN - the offset for the indices in IXELEM
C   --   IXELEM   - IN - the indices of the output elements
C   --                   (iff NUMEBO <> NUMELB)
C   --   NUMLNK   - IN - the number of nodes per element
C   --   NUMATR   - IN - the number of attributes
C   --   LINK     - IN - the element connectivity for this block
C   --   ATRIB    - IN - the attributes for this block
C   --   BLKTYP   - IN - the element block type
C   --   IDELB    - IN - element block ID
C   --   NEWNOD   - IN - true iff nodes are renumbered
C   --   NODIX(i) - IN - the zoom mesh index for each node (iff NEWNOD)
C   --   LNSCR    - IN - connectivity scratch array
C   --   ATRSCR   - IN - attribute scratch array

      include 'exodusII.inc'

      INTEGER NDBOUT
      INTEGER NUMELB
      INTEGER NUMEBO
      INTEGER IEL
      INTEGER IXELEM(*)
      INTEGER NUMLNK
      INTEGER NUMATR
      INTEGER NLNKDM
      INTEGER NATRDM
      INTEGER LINK(NLNKDM,NUMELB)
      REAL ATRIB(NATRDM,NUMELB)
      CHARACTER*(MXSTLN) BLKTYP
      INTEGER IDELB
      LOGICAL NEWNOD
      INTEGER NODIX(*)
      INTEGER LNSCR(NLNKDM,*)
      REAL ATRSCR(NATRDM,*)

      IEL0 = IEL - 1
      IF ((NUMLNK .GT. 0) .AND. (NUMEBO .GT. 0)) THEN
         IF (NUMELB .EQ. NUMEBO) THEN
            IF (.NOT. NEWNOD) THEN
               CALL EXPELC(NDBOUT, IDELB, LINK, IERR)
            ELSE
               DO 20 i=1, numlnk
                  DO 10 ne = 1, numelb
                     lnscr(i, ne) = nodix(link(i,ne))
 10               CONTINUE
 20            CONTINUE
               CALL EXPELC(NDBOUT, IDELB, LNSCR, IERR)
            END IF
         ELSE
            IF (.NOT. NEWNOD) THEN
               do 40 i=1,numlnk
                  do 30 ix=1,numebo
                     lnscr(i,ix) = link(i,ixelem(ix)-iel0)
 30               continue
 40            continue
            ELSE
               do 60 i=1,numlnk
                  do 50 ix=1,numebo
                     lnscr(i,ix) = nodix(link(i,ixelem(ix)-iel0))
 50               continue
 60            continue
            END IF
            CALL EXPELC(NDBOUT, IDELB, LNSCR, IERR)
         END IF
      END IF

      IF ((NUMATR .GT. 0) .AND. (NUMEBO .GT. 0)) THEN
         IF (NUMELB .EQ. NUMEBO) THEN
            call expeat(ndbout, idelb, atrib, ierr)
         ELSE
            DO 80 i=1, numatr
               DO 70 ix=1,numebo
                  atrscr(i,ix) = atrib(i,ixelem(ix)-iel0)
 70            CONTINUE
 80         CONTINUE
            call expeat(ndbout, idelb, atrscr, ierr)
         END IF
      END IF

      RETURN
      END
