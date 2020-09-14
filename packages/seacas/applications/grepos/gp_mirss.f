C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MIRSS (IDESS, NEESS, IXEESS, LTEESS, LTSSS,
     $     IBLOCK, BLKTYP, ALLONE, COMTOP)
C=======================================================================

C   --*** MIRSS *** (GJOIN) Mirror element side sets
C   --   Written by Greg Sjaardema
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
C   --   IXEESS - IN/OUT - the index of the first element for each set
C   --   LTEESS - IN/OUT - the elements for all sets
C   --   LTSSS - IN/OUT - the sides for all sets
C   --   ALLONE - IN - true if df are all unity, false if not all unity

      include 'gp_params.blk'
      include 'gp_dbnums.blk'
      INTEGER IDESS(*)   ! NUMESS
      INTEGER NEESS(*)   ! NUMESS
      INTEGER IXEESS(*)  ! NUMESS
      INTEGER LTEESS(*)  ! LESSEL
      INTEGER LTSSS(*)   ! LESSEL
      INTEGER IBLOCK(*)  ! NELBLK - gives index of last element in block
      CHARACTER*(MXSTLN) BLKTYP(*)

      LOGICAL ALLONE
      CHARACTER*(*)   COMTOP

C     The routine now provides more checking and will return with a warning
C     if applied to elements on other faces...

      if (comtop(:8) .eq. 'MULTIPLE') THEN
         DO NL = 1, NUMESS
            iel = IXEESS(NL)
            call mirs2(NEESS(NL), LTEESS(iel),
     *           LTSSS(iel), iblock, blktyp, nelblk, ndim)
         end do
         return
      else
C     ... comtop was not equal to "MULTIPLE_TOPOLOGIES", so
C         at this point, the underlying element topology for all elements
C         is the same....
         DO NL = 1, NUMESS
            iel = IXEESS(NL)
            call mirs1(NEESS(NL), LTSSS(iel), blktyp(1), ndim)
         end do
      end if
      IF (.NOT. ALLONE) THEN
         CALL PRTERR ('WARNING',
     *        'Mirroring of sideset distribution factors not supported')
      END IF

      RETURN
      END

      subroutine mirs1(numsid, side, type, ndim)

      include 'gp_params.blk'

      integer numsid
      integer side(*)
      CHARACTER*(MXSTLN) TYPE

      if (type(:3) .eq. 'HEX') then
         do i = 1, numsid
            if (side(i) .le. 4) then
               side(i) = 5 - side(i)
            else
               side(i) = side(i)
            end if
         end do
      else if (type(:4) .eq. 'QUAD') then
         do i = 1, numsid
            side(i) = 5 - side(i)
         end do
      else if (type(:5) .eq. 'SHELL') then
         return
      else if (type(:5) .eq. 'TRUSS' .or.
     $         type(:4) .eq. 'BEAM') then
         return
      else if (type(:3) .eq. 'TET') then
         do i = 1, numsid
            if (side(i) .eq. 3) then
               side(i) = 1
            else if (side(i) .eq. 1) then
               side(i) = 3
            end if
         end do
      else if (type(:3) .eq. 'TRI') then
C ... 'TRI' can be triangle or trishell
         if (ndim .eq. 2) then
            do i = 1, numsid
               side(i) = 4 - side(i)
            end do
         end if
      else if (type(:5) .eq. 'WEDGE') then
        do i = 1, numsid
          if (side(i) .eq. 3) then
            side(i) = 1
          else if (side(i) .eq. 1) then
            side(i) = 3
          end if
        end do
      end if
      return
      end

      subroutine mirs2(numsid, elems, side,
     $     iblock, blktyp, nelblk, ndim)

      include 'gp_params.blk'
      integer numsid, nelblk, elem
      integer elems(*), side(*), iblock(*)
      CHARACTER*(MXSTLN) BLKTYP(*)

      lstelem = 0
      lstblk = 1
      iblk = 0
      do i = 1, numsid
         elem = elems(i)
         if (elem .gt. lstelem) then
            do j = lstblk, nelblk
               if (elem .le. iblock(j)) then
                  iblk = j;
                  goto 30
               end if
            end do
 30         continue
         else
            do j = 1, nelblk
               if (elem .le. iblock(j)) then
                  iblk = j;
                  goto 20
               end if
            end do
 20         continue
         end if
         lstelem = elem
         lstblk  = iblk

         if (blktyp(iblk)(:3) .eq. 'HEX') then
            if (side(i) .le. 4) then
               side(i) = 5 - side(i)
            else
               side(i) = side(i)
            end if
         else if (blktyp(iblk)(:4) .eq. 'QUAD') then
            side(i) = 5 - side(i)
         else if (blktyp(iblk)(:5) .eq. 'SHELL') then
            side(i) = side(i)
         else if (blktyp(iblk)(:5) .eq. 'TRUSS' .or.
     $            blktyp(iblk)(:4) .eq. 'BEAM') then
            side(i) = side(i)
         else if (blktyp(iblk)(:3) .eq. 'TET') then
            if (side(i) .eq. 3) then
               side(i) = 1
            else if (side(i) .eq. 1) then
               side(i) = 3
            end if
         else if (blktyp(iblk)(:3) .eq. 'TRI') then
C ... 'TRI' can be triangle or trishell
           if (ndim .eq. 2) then
             side(i) = 4 - side(i)
           end if
         else if (blktyp(iblk)(:5) .eq. 'WEDGE') then
           if (side(i) .eq. 3) then
             side(i) = 1
           else if (side(i) .eq. 1) then
             side(i) = 3
           end if
         end if
C ... If not in list above, assume that they need no adjustment.
      end do
      return
      end
