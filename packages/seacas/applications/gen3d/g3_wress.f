C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRESS (A, IA, IDFRO, IDBCK,
     &     ISSFRO, ISSBCK, NSSUR, NESUR, NSSFRO, NSSBCK,
     &     IDESS, NEES3, NNES3, IXEES3, IXNES3, LTEES3, LTSES3,
     &     LTNES3, FACES3, *)
C=======================================================================

C   --*** WRESS *** (GEN3D) Write 3D node sets
C   --   Written by Amy Gilkey - revised 05/05/86
C   --
C   --WRESS writes the side set information for the 3D database.
C   --Calculations have been done elsewhere.
C   --
C   --Parameters:
C   --   IDFRO - IN - ids for front surface side sets; (0) = length
C   --   IDBCK - IN - ids for back surface side sets; (0) = length
C   --   ISSFRO - IN - the elements in the front surface side set
C   --   ISSBCK - IN - the elements in the back surface side set
C   --   NSSUR - IN - the number of nodes in the surface side set
C   --   NESUR - IN - the number of elements in the surface side set
C   --   NSSFRO - IN - the nodes in the front surface side set
C   --   NSSBCK - IN - the nodes in the back surface side set
C   --   IDESS - IN - the ids for each 2D set
C   --   NEES3 - IN - the number of elements for each 3D set
C   --   NNES3 - IN - the number of nodes for each 3D set
C   --   IXEES3 - IN - the index of the first element for each 3D set
C   --   IXNES3 - IN - the index of the first node for each 3D set
C   --   LTEES3 - IN - the elements for all 3D sets
C   --   LTSES3 - IN - the element sides for all 3D sets
C   --   LTNES3 - IN - the nodes for all 3D sets
C   --   FACES3 - IN - the distribution factors for all 3D sets
C   --
C   --Common Variables:
C   --   Uses NDBOUT of /DBASE/
C   --   Uses NUMESS, LESSEL, LESSNL of /DBNUMS/
C   --   Uses LESSEO, LESSNO of /DBNUM3/
C   --   Uses NNREPL, NEREPL of /PARAMS/

      INCLUDE 'exodusII.inc'
      INCLUDE 'g3_dbase.blk'
      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'

      REAL    A(*)
      INTEGER IA(*)
      INTEGER IDFRO(0:*)
      INTEGER IDBCK(0:*)
      INTEGER ISSFRO(NESUR), ISSBCK(NESUR)
      INTEGER NSSFRO(*), NSSBCK(*)
      INTEGER IDESS(*)
      INTEGER NEES3(*)
      INTEGER NNES3(*)
      INTEGER IXEES3(*)
      INTEGER IXNES3(*)
      INTEGER LTEES3(*)
      INTEGER LTSES3(*)
      INTEGER LTNES3(*)
      REAL FACES3(*)

      LOGICAL ANYESS

      NFRO = IDFRO(0)
      NBCK = IDBCK(0)
      ANYESS = (NFRO .GT. 0) .OR. (NBCK .GT. 0) .OR. (NUMESS .GT. 0)

C   --Write 3D
      iend = 0
      do 10 iss = 1, numess
        istart = iend + 1
        iend   = istart + nnes3(iss) - 1
        call expsp (ndbout, idess(iss), nees3(iss), nnes3(iss), ierr)
        if (ierr .lt. 0) then
           call exerr('gen3d2', 'Error from expsp', exlmsg)
           go to 40
        endif
        call expss (ndbout, idess(iss), ltees3(ixees3(iss)),
     &       ltses3(ixees3(iss)), ierr)
        if (ierr .lt. 0) then
           call exerr('gen3d2', 'Error from expss', exlmsg)
           go to 40
        endif
        call expssd(ndbout, idess(iss), faces3(istart), ierr)
        if (ierr .lt. 0) then
           call exerr('gen3d2', 'Error from expssd', exlmsg)
           go to 40
        endif
 10   continue

C ... Assume distribution factors are 1.0 for all front and back sidesets
C     Need to build a temporary array to hold the '1.0's
C     The size of the array is MAX(NSSFRO, NSSBCK)
      if (nfro .gt. 0 .or. nbck .gt. 0) then
         call mdrsrv('TDIST', ktdist, nssur)
         call mdrsrv('ISIDE', kiside, nesur)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 40
         call inirea(nssur, 1.0, a(ktdist))
C ... Front sidesets are surface 6
         call iniint(nesur, 6, ia(kiside))
C ... If the element number is negative, then we use surface 5
C     See newess for the code that sets the element number negative
         if (nfro .gt. 0) then
           do i = 1, nesur
             if (issfro(i) .lt. 0) then
               ia(kiside+i-1) = 5
               issfro(i) = -issfro(i)
             end if
           end do

         end if
         do iss = 1, nfro
           call expsp (ndbout, idfro(iss), nesur, nssur, ierr)
           if (ierr .lt. 0) then
              call exerr('gen3d2', 'Error from expsp', exlmsg)
              go to 40
           endif
           call expss (ndbout, idfro(iss), issfro, ia(kiside), ierr)
           if (ierr .lt. 0) then
              call exerr('gen3d2', 'Error from expss', exlmsg)
              go to 40
           endif
           call expssd(ndbout, idfro(iss), a(ktdist), ierr)
           if (ierr .lt. 0) then
              call exerr('gen3d2', 'Error from expssd', exlmsg)
              go to 40
           endif
         end do

C ... Back sidesets are surface 5
         call iniint(nesur, 5, ia(kiside))
C ... If the element number is negative, then we use surface 4
C     See newess for the code that sets the element number negative
         if (nbck .gt. 0) then
           do i = 1, nesur
             if (issbck(i) .lt. 0) then
               ia(kiside+i-1) = 4
               issbck(i) = -issbck(i)
             end if
           end do
         end if

         do iss = 1, nbck
            call expsp (ndbout, idbck(iss), nesur, nssur, ierr)
           if (ierr .lt. 0) then
              call exerr('gen3d2', 'Error from expsp', exlmsg)
              go to 40
           endif
           call expss (ndbout, idbck(iss), issbck, ia(kiside), ierr)
           if (ierr .lt. 0) then
              call exerr('gen3d2', 'Error from expss', exlmsg)
              go to 40
           endif
           call expssd(ndbout, idbck(iss), a(ktdist), ierr)
           if (ierr .lt. 0) then
              call exerr('gen3d2', 'Error from expssd', exlmsg)
              go to 40
           endif
         end do
         call mddel('TDIST')
         call mddel('ISIDE')
      end if

      RETURN
C ... Control passes here if any memory or other errors
 40   continue
      return 1
      END
