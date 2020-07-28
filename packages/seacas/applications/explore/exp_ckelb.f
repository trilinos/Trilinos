C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE CKELB (NELBLK, NUMEL, NUMNP, EBTYPE,
     &   IDELB, NUMELB, NUMLNK, NUMATR, LINK, ATRIB, ATNAME,
     &   NODUSE, MAPNO)
C=======================================================================

C   --*** CKELB *** (EXPLORE) Check database element blocks
C   --
C   --CKELB checks the element block information.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   NUMEL - IN - the number of elements in all blocks
C   --   NUMNP - IN - the number of nodes
C   --   IDELB - IN - the element block ID for each block
C   --   NUMELB - IN - the number of elements for each block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   NUMATR - IN - the number of attributes for each block
C   --   LINK - IN - the connectivity array
C   --   ATRIB - IN - the attribute array

      include 'exp_errcnt.blk'
      include 'exodusII.inc'
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      INTEGER LINK(*)
      REAL ATRIB(*)
      INTEGER NODUSE(*)
      INTEGER MAPNO(*)
      CHARACTER*(MXSTLN) EBTYPE(*)
      CHARACTER*(*) ATNAME(*)

      CHARACTER*99 STRA
      CHARACTER*99 STRB

      if (numel .eq. 0 .and. nelblk .eq. 0) return

C   --Check for unique identifier

      DO 100 IELB = 1, NELBLK
         IF (LOCINT (IDELB(IELB), IELB-1, IDELB) .GT. 0) THEN
            CALL INTSTR (1, 0, IDELB(IELB), STRA, LSTRA)
            CALL PRTERR ('CMDSPEC',
     &         'Element block ID ' // STRA(:LSTRA) // ' is not unique')
         END IF
  100 CONTINUE

      iatoff = 0
      do ielb = 1, nelblk
        do natr = 1, numatr(ielb)
          if (atname(iatoff+natr)(1:1) .ne. ' ') then
            do jatr = natr+1, numatr(ielb)
              if (atname(iatoff+natr) .eq. atname(iatoff+jatr)) then
                CALL INTSTR (1, 0, IDELB(IELB), STRA, LSTRA)
                CALL PRTERR ('CMDSPEC',
     &            'Attribute name "'//
     *            ATNAME(iatoff+NATR)(:LENSTR(ATNAME(iatoff+NATR))) //
     *            '" in Element block ID ' // STRA(:LSTRA)
     *            // ' is not unique')
              end if
            end do
          end if
        end do
        iatoff = iatoff + numatr(ielb)
      end do

C   --Check that the number of elements in element blocks is equal to the total
C   --number of elements

      NELB = INTADD (NELBLK, NUMELB)
      IF (NELB .NE. NUMEL) THEN
         CALL PRTERR ('CMDSPEC',
     *    'Sum of elements in all element blocks does not match total')
      END IF

C   --Check the connectivity for each element block

C .. Initialize noduse array to 0.
      call iniint (numnp, 0, noduse)

      IEL0 = 0
      ISLNK = 1
      DO 110 IELB = 1, NELBLK
        if (ebtype(ielb) .eq. 'nsided' .or.
     *    ebtype(ielb) .eq. 'NSIDED') THEN
         CALL CKEB1 (IEL0,
     &      IELB, IDELB(IELB), 1, NUMLNK(IELB), NUMNP,
     &      LINK(ISLNK), noduse)
         ISLNK = ISLNK + NUMLNK(IELB)
        else
         CALL CKEB1 (IEL0,
     &      IELB, IDELB(IELB), NUMELB(IELB), NUMLNK(IELB), NUMNP,
     &      LINK(ISLNK), noduse)
         ISLNK = ISLNK + NUMLNK(IELB) * NUMELB(IELB)
       end if
         IEL0 = IEL0 + NUMELB(IELB)
  110 CONTINUE

C ... Check whether all nodes have been specified in an elements connectivity.
      nerr = 0
      DO 120 I = 1, numnp
         if (noduse(i) .ne. 1) then
            if (nerr .lt. maxerrs .or. maxerrs .le. 0) then
               CALL INTSTR (1, 0, I, STRA, LSTRA)
               CALL INTSTR (1, 0, MAPNO(I), STRB, LSTRB)
               CALL PRTERR ('CMDSPEC',
     &              'Local Node ' // STRA(:LSTRA) //
     $              ' (Global Node ' // STRB(:LSTRB) //
     &              ') is not connected to any elements')
            else if (nerr .eq. maxerrs .and. maxerrs .gt. 0) then
               call prterr('CMDSPEC',
     $              '...skipping additional errors...'//
     *          ' (use command "maxerr #" to change)')
            end if
            nerr = nerr + 1
         end if
 120  continue
      if (nerr .gt. 0) then
         write (stra, 10020) nerr
         call sqzstr(stra, lstra)
         CALL PRTERR ('CMDSPEC', STRA(:lstra))
      end if
10020 FORMAT('NODE USE ERROR: Found ',I12,
     $     ' nodes that are not connected to any elements.')
      RETURN
      END
