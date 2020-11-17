C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRESS (A, IA, IDFRO, IDBCK,
     &   ISSFRO, ISSBCK, NSSUR, NSSFRO, NSSBCK,
     &   IDESS, NEES3, NNES3, IXEES3, IXNES3, LTEES3, LTSSS3, FACES3, *)
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
C   --   NSSFRO - IN - the nodes in the front surface side set
C   --   NSSBCK - IN - the nodes in the back surface side set
C   --   IDESS - IN - the ids for each 2D set
C   --   NEES3 - IN - the number of elements for each 3D set
C   --   NNES3 - IN - the number of nodes for each 3D set
C   --   IXEES3 - IN - the index of the first element for each 3D set
C   --   IXNES3 - IN - the index of the first node for each 3D set
C   --   LTEES3 - IN - the elements for all 3D sets
C   --   LTSSS3 - IN - the element sides for all 3D sets
C   --   FACES3 - IN - the distribution factors for all 3D sets
C   --
C   --Common Variables:
C   --   Uses NDBOUT of /DBASE/
C   --   Uses NUMESS, LESSEL, LESSNL of /DBNUMS/
C   --   Uses LESSEO, LESSNO of /DBNUM3/
C   --   Uses NNREPL, NEREPL of /PARAMS/

      include 'exodusII.inc'
      INCLUDE 'gs_dbase.blk'
      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'

      REAL A(*)
      INTEGER IA(*)
      INTEGER IDFRO(0:*)
      INTEGER IDBCK(0:*)
      INTEGER ISSFRO(NUMEL), ISSBCK(NUMEL)
      INTEGER NSSFRO(*), NSSBCK(*)
      INTEGER IDESS(*)
      INTEGER NEES3(*)
      INTEGER NNES3(*)
      INTEGER IXEES3(*)
      INTEGER IXNES3(*)
      INTEGER LTEES3(*)
      INTEGER LTSSS3(*)
      REAL FACES3(*)

      LOGICAL ANYESS

      NFRO = IDFRO(0)
      NBCK = IDBCK(0)
      ANYESS = (NFRO .GT. 0) .OR. (NBCK .GT. 0) .OR. (NUMESS .GT. 0)

C   --Write 3D Sidesets
      IF (.NOT. ANYESS) RETURN

C     Each sideset has numel elements and faces, and NSSUR distribution factors
      if ((nfro .gt. 0) .or. (nbck .gt. 0)) then
C ... Allocate array to store side data for front and/or back sidesets
         call mdrsrv('LSTSID', klst, numel)
C ... Allocate array for distribution factors
         call mdrsrv('DISTF', kdistf, nssur)
         call mdstat(nerr, mem)
         if (nerr .gt. 0) go to 400
         call inirea(nssur, 1.0, a(kdistf))

C ... Write Front Sidesets
         if (nfro .gt. 0) then
C ... Fill the side list
            call iniint(numel, 1, ia(klst))
         end if
         do 100 i=1, nfro
           call expsp  (ndbout, idfro(i), numel, nssur, ierr)
           call expss  (ndbout, idfro(i), issfro, ia(klst), ierr)
           call expssd (ndbout, idfro(i), a(kdistf), ierr)
 100     continue

C ... Write Back  Sidesets
         if (nbck .gt. 0) then
C ... Fill the side list
            call iniint(numel, 2, ia(klst))
         end if
         do 200 i=1, nbck
           call expsp  (ndbout, idbck(i), numel, nssur, ierr)
           call expss  (ndbout, idbck(i), issbck, ia(klst), ierr)
           call expssd (ndbout, idbck(i), a(kdistf), ierr)
 200     continue
         call mddel('LSTSID')
         call mddel('DISTF')
      end if

C ... Write 2D-generated Sidesets
C ... Fix the side list. Assuming input is quads, then the
C     following should hold newside = oldside + 2
      do 300 i=1, lesseo
        ltsss3(i) = ltsss3(i) + 2
 300  continue

      do 310 i=1, numess
        call expsp  (ndbout, idess(i), nees3(i), nnes3(i), ierr)
        call expss  (ndbout, idess(i), ltees3(ixees3(i)),
     &       ltsss3(ixees3(i)), ierr)
        call expssd (ndbout, idess(i), faces3(ixnes3(i)), ierr)
 310  continue

      RETURN

 400  CONTINUE
      RETURN 1
      END
