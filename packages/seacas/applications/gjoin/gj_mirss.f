C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE MIRSS (NUMESS, LESSEL, LESSDL, IDESS, NEESS, NEDSS,
     *     IXEESS, IXEDSS, LTEESS, LTSSS, LTSNC, FAC, USESDF, NONQUAD,
     *     COMTOP)
C=======================================================================
C   --*** MIRSS *** (GJOIN) Mirror element side sets
C   --
C   --MIRSS mirrors a side set and (if USESDF true) the distribution factors
C   --applied to the nodes.
C   --
C   --Parameters:
C   --
C   --   NUMESS - IN/OUT - the number of element side sets
C   --   LESSEL - IN/OUT - the length of the element side sets element list
C   --   IDESS - IN/OUT - the element side set ID for each set
C   --   NEESS - IN/OUT - the number of elements for each set
C   --   NEDSS - IN/OUT - the number of dist-fac for each set
C   --   IXEESS - IN/OUT - the index of the first element for each set
C   --   IXEDSS - IN/OUT - the index of the first dist-fac for each set
C   --   LTEESS - IN/OUT - the elements for all sets
C   --   LTSSS - IN/OUT - the sides for all sets
C   --   LTSNC - IN/OUT - the face count for each element/side in the list
C   --   FACESS - IN/OUT - the distribution factors for all sets????????????
C   --   USESDF - IN - true if df are non-unity, false if all unity
C   --   NONQUAD - IN - true if model contains non-hex/non-quad elements

      INTEGER IDESS(*)   ! NUMESS
      INTEGER NEESS(*)   ! NUMESS
      INTEGER NEDSS(*)   ! NUMESS
      INTEGER IXEESS(*)  ! NUMESS
      INTEGER IXEDSS(*)  ! NUMESS
      INTEGER LTEESS(*)  ! LESSEL
      INTEGER LTSSS(*)   ! LESSEL
      INTEGER LTSNC(*)   ! LESSEL
      REAL    FAC(*)     ! LESSDL
      LOGICAL USESDF, NONQUAD, shells
      CHARACTER*(*)   COMTOP

C ... This routine was originally written to only handle quads, tris, and hexes.
C     There was no checking of this, it blindly went through the list
C     swapping nodes 1 and 2 (line) or nodes 4 and 2 (quad face).

C     The routine now provides more checking and will return with a warning
C     if applied to elements on other faces...

      IF (NONQUAD .AND. COMTOP(:3) .NE. 'TET' .and.
     *  COMTOP(:5) .ne. 'SHELL') THEN
        CALL PRTERR ('PROGRAM',
     *    'Mirroring of sidesets on non-quad/line sides not supported')
        RETURN
      END IF

      IF (USESDF) THEN
        CALL PRTERR ('WARNING',
     *    'Mirroring of sideset distribution factors not supported')
        RETURN
      END IF

C ... We have a quad/line face; do the mirroring.
      if (comtop(:5) .eq. 'SHELL') then
        shells = .true.
      else
        shells = .false.
      end if

      DO 10 NL = 1, NUMESS
        call mirs1(IDESS(NL), NEESS(NL), NEDSS(NL), LTEESS(IXEESS(NL)),
     *    LTSSS(IXEESS(NL)), LTSNC(IXEESS(NL)), FAC(IXEDSS(NL)),
     *    USESDF, shells)
   10 CONTINUE
      RETURN
      END

      subroutine mirs1(id, numsid, numdis, elem, side, dfcnt, facedf,
     *  usesdf, shells)

      integer id, numsid, numdis
      integer elem(*), side(*), dfcnt(*)
      real    facedf(*)
      logical usesdf
      logical shells
      CHARACTER*132 STRING

      idf = 0
      do 10 i = 1, numsid
        icnt = dfcnt(i)

C ... Bar topology side -- Base element is quad (or we wouldn't be here)
        if (icnt .eq. 2) then
          side(i) = 5 - side(i)

C ... Quad topology side -- Base element is hex (or we wouldn't be here)
        else if (icnt .eq. 4) then
          if (shells) then
            side(i) = side(i)
          else
            if (side(i) .le. 4) then
              side(i) = 5 - side(i)
            else
              side(i) = side(i)
            end if
          end if
        else if (icnt .eq. 3) then
C ... NOTE: This is a NONQUD, but COMTOP has been checked in calling code
C           so we know that these are all tets or shells...
          if (shells) then
            side(i) = side(i)
          else
           if (side(i) .eq. 3) then
              side(i) = 1
           else if (side(i) .eq. 1) then
              side(i) = 3
            end if
          end if
        else
          WRITE (STRING, 100) ID, ICNT
 100      FORMAT('Sideset ',I5,' contains ',I2,'-node ',
     *      ' sides which are not supported for mirroring by gjoin2')
          CALL SQZSTR (STRING, LSTR)
          CALL PRTERR ('PROGRAM', STRING(:LSTR))
          RETURN
        end if
 10   continue
      return
      end
