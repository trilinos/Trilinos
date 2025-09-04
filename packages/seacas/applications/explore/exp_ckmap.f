C    Copyright(C) 1999-2020, 2024 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE CKMAP (ICNT, MAP, INDX, TYPE)
C=======================================================================

C   --*** CKMAP *** (EXPLORE) Check database node/element map
C   --
C   --CKMAP checks the node/element order map.
C   --
C   --Parameters:
C   --   ICNT - IN - the number of nodes/elements
C   --   MAP  - IN - the node/element order map
C   --   INDX - SCRATCH - size = ICNT

      include 'exp_errcnt.blk'
      INTEGER MAP(*)
      INTEGER INDX(*)
      CHARACTER*(*) TYPE

      CHARACTER*1024 STRA
C   --Check that each node/element appears once and only once in the map
C     The 'map(i)' values may be larger than icnt, so we can't do a
C     simple check.  Instead, we do an indexed sort and check for no
C     duplicate adjacent values.

      nerr = 0

      ibad = 0
      do iel = 1, icnt
         if (map(iel) .le. 0) then
            ibad = ibad + 1
            indx(ibad) = iel
         end if
      end do
      if (ibad .gt. 0) then
         write (stra, 10002) ibad, type
10002    FORMAT('MAP ERROR: There are ',I12,' ',A,
     $        ' with invalid global ids (<= 0): ')
         call sqzstr(stra, lstra)
         CALL PRTERR ('CMDSPEC', STRA(:lstra))
         call formlist(indx, ibad)
      end if

      if (icnt .le. 1) return

      CALL INDEXI (MAP, INDX, ICNT, .TRUE.)

C ... There has been a request to show min and max ids to help with
C     debugging potential database corruption issues.  Do it here.

      if (indx(1) .eq. 1 .and. indx(icnt) .eq. icnt .and.
     $     (map(icnt) - map(1) +1 .eq. icnt)) then
         write (stra, 10001) type, map(indx(1)), map(indx(icnt)),
     $        'sequential'
      else
         write (stra, 10001) type, map(indx(1)), map(indx(icnt)),
     $        'NOT sequential'
      endif
10001 FORMAT('INFO: ', A, ' global id range: ',I12, ' to ', I12,
     $     '. Map is ', A)
      call sqzstr(stra, lstra)
      CALL PRTERR ('CMDSPEC', STRA(:lstra))

      ILAST = MAP(INDX(1))
      DO IEL = 2, ICNT
         if (map(indx(iel)) .eq. ilast) then
            if (ilast .ne. 0) then
               if (nerr .lt. maxerrs .or. maxerrs .le. 0) then
                  write (stra, 10000) type, ilast, type,
     *                 indx(iel-1), indx(iel)
10000             FORMAT('MAP ERROR: ',A,' global  id ',I12,
     *                 ' assigned to ',A,'s', I12,' and ',I12,'.')
                  call sqzstr(stra, lstra)
                  CALL PRTERR ('CMDSPEC', STRA(:lstra))
               else if (nerr .eq. maxerrs .and. maxerrs .gt. 0) then
                  call prterr('CMDSPEC',
     $                 '...skipping additional errors...')
               end if
               nerr = nerr + 1
            end if
         end if
         ilast = map(indx(iel))
      END DO
      if (ibad .gt. 0 .or. nerr .gt. 0) then
         write (stra, 10010) ibad, nerr, type
10010    FORMAT('MAP ERROR: Found ',I12, ' bad ids and ',I12,
     $        ' other errors in ',A,' map check.')
         call sqzstr(stra, lstra)
         CALL PRTERR ('CMDSPEC', STRA(:lstra))
      end if
      RETURN
      END
