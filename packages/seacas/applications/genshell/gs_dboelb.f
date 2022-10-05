C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBOELB (A, NDB,
     &   IDELB, NUMELB, NUMLNK, NUMATR, LINK, NAMELB, ATRIB, ATRIBNW)
C=======================================================================

C   --*** DBOELB *** (EXOLIB) Write database element blocks
C   --   Written by Amy Gilkey - revised 10/12/87
C   --
C   --DBOELB writes the element block information to the database.
C   --Some dynamic dimensioning is done.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NELBS, NELBE - IN - the number of first and last element blocks
C   --      to write
C   --   IDELB - IN - the element block IDs for each block
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   NUMATR - IN - the number of attributes in each block
C   --   LINK - IN - the connectivity for each block
C   --   ATRIB - IN - the attributes for each block
C   --   ATRIBNW - IN - the new attributes if block is a 3D beam
      include 'exodusII.inc'

      REAL A(*)
      INTEGER NDB
      INTEGER IDELB
      INTEGER NUMELB
      INTEGER NUMLNK
      INTEGER NUMATR
      INTEGER LINK(*)
      REAL ATRIB(*)
      REAL ATRIBNW(7)
      CHARACTER*(mxstln) NAMELB

      IELNK = 0
      IEATR = 0

      call expelb(ndb, IDELB, NAMELB, NUMELB, NUMLNK, NUMATR, IERR)
      call expelc(ndb, idelb, link, ierr)

      if (numatr .gt. 0) then
        if (namelb(:4) .eq. 'BEAM') then
C ... A 3D beam needs special treatment since it has 7 attributes
C     and the input 2D beam will only have 1 or 3 attributes.
C     The attributes have been specified in the 'ATRIBNW' array
C     by the user (or the defaults are used).  Need to expand the
C     single value per block values from the ATRIBNW array into
C     7 values per element.
          call mdlong('ATRIBNW', KATRIB, 7*numelb)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) THEN
            CALL MEMERR
            STOP 'memory error'
          END IF
          call geneat(ndb, idelb, a(katrib), atribnw, numelb)
        else
          call expeat(ndb, idelb, atrib, ierr)
        end if
      end if

      RETURN
      END

      subroutine geneat(ndb, idelb, atrib, atribnw, numelb)
      real atrib(*)
      real atribnw(7)

      i = 0
      do 20 ie = 1, numelb
          do 10 ia=1, 7
            i = i + 1
            atrib(i) = atribnw(ia)
 10     continue
 20   continue
      call expeat(ndb, idelb, atrib, ierr)

      return
      end

