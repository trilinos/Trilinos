C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRMAP (NTXT, OPTION, NUMNP, NUMEL,
     &                  NPMAP, ELMAP, MAPEL)
C=======================================================================

C   --*** WRMAP *** (EXOTXT) Write database node number map,
C   --                       element number map, and/or element order map
C   --   Written by Amy Gilkey - revised 02/27/86
C   --   Modified for ExodusIIv2 - 10/24/95
C   --
C   --Parameters:
C   --   NTXT   - IN - the database text file
C   --   OPTION - IN - '*' write all
C   --                 'N' write node number map
C   --                 'E' write element number map
C   --                 'O' write element order map
C   --   NUMNP  - IN - number of nodes
C   --   NUMEL  - IN - number of elements
C   --   NPMAP  - IN - node number map (if OPTION)
C   --   ELMAP  - IN - element number map (if OPTION)
C   --   MAPEL  - IN - element order map (if OPTION)

      INTEGER NTXT
      CHARACTER*(*) OPTION
      INTEGER NUMNP
      INTEGER NUMEL
      INTEGER NPMAP(*)
      INTEGER ELMAP(*)
      INTEGER MAPEL(*)
      LOGICAL ALL, NOPT, EOPT, OOPT
      LOGICAL INORDR

      ALL  = (OPTION .EQ. '*')
      NOPT = (INDEX (OPTION, 'N') .GT. 0)
      EOPT = (INDEX (OPTION, 'E') .GT. 0)
      OOPT = (INDEX (OPTION, 'O') .GT. 0)

C     Write node number map
      IF (ALL .OR. NOPT) THEN
         WRITE (NTXT, 1000) '! Node number map'
         if (inordr(npmap, numnp)) then
           write (ntxt, 1000) 'sequence 1..numnp'
         else
           write (ntxt, 1000) 'explicit map'
           WRITE (NTXT, 1010) (NPMAP(I), I = 1, NUMNP)
         end if
      END IF

C     Write element number map
      IF (ALL .OR. EOPT) THEN
         WRITE (NTXT, 1000) '! Element number map'
         if (inordr(elmap, numel)) then
           write (ntxt, 1000) 'sequence 1..numel'
         else
           write (ntxt, 1000) 'explicit map'
           WRITE (NTXT, 1010) (ELMAP(I), I = 1, NUMEL)
         end if
      END IF

C     Write element order map
      IF (ALL .OR. OOPT) THEN
         WRITE (NTXT, 1000) '! Element order map'
         if (inordr(mapel, numel)) then
           write (ntxt, 1000) 'sequence 1..numel'
         else
           write (ntxt, 1000) 'explicit map'
           WRITE (NTXT, 1010) (MAPEL(I), I = 1, NUMEL)
         end if
      END IF

 1000 FORMAT (A)
 1010 FORMAT (8I10)

      RETURN
      END

C=======================================================================
      logical function inordr(MAP, ISIZE)
C=======================================================================

C ... Determine if the passed in map is simply a sequence from 1..isize

      integer map(isize)

      inordr = .FALSE.
      do 10 i=1, isize
        if (map(i) .ne. i) then
          inordr = .FALSE.
          return
        end if
 10   continue
      inordr = .TRUE.
      return
      end
