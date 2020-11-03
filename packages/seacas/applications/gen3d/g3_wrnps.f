C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRNPS (A, IA, IDFRO, IDBCK,
     &   IDNPS, NNNP3, IXNNP3, LTNNP3, FACNP3,
     &   IXNP, NRNP, *)
C=======================================================================

C   --*** WRNPS *** (GEN3D) Write 3D node sets
C   --   Written by Amy Gilkey - revised 05/05/86
C   --
C   --WRNPS writes the node set information for the 3D database.
C   --Calculations have been done elsewhere.
C   --
C   --Parameters:
C   --   IDFRO - IN - ids for front surface node sets; (0) = length
C   --   IDBCK - IN - ids for back surface node sets; (0) = length
C   --   IDNPS - IN - the 2D node sets ids
C   --   NNNP3 - IN - the number of nodes for each 3D set
C   --   IXNNP3 - IN - the index of the first node for each 3D set
C   --   LTNNP3 - IN - the nodes for all 3D sets
C   --   FACNP3 - IN - the distribution factors for all 3D sets
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --
C   --Common Variables:
C   --   Uses NDBOUT of /DBASE/
C   --   Uses NUMNPS, LNPSNL of /DBNUMS/
C   --   Uses LNPSNO of /DBNUM3/

      INCLUDE 'exodusII.inc'
      INCLUDE 'g3_dbase.blk'
      INCLUDE 'g3_dbnums.blk'
      INCLUDE 'g3_dbnum3.blk'

      REAL    A(*)
      INTEGER IA(*)
      INTEGER IDFRO(0:*)
      INTEGER IDBCK(0:*)
      INTEGER IDNPS(*)
      INTEGER NNNP3(*)
      INTEGER IXNNP3(*)
      INTEGER LTNNP3(*)
      REAL FACNP3(*)
      INTEGER IXNP(*), NRNP(*)

      LOGICAL ANYNPS

      NFRO = IDFRO(0)
      NBCK = IDBCK(0)
      ANYNPS = (NFRO .GT. 0) .OR. (NBCK .GT. 0) .OR. (NUMNPS .GT. 0)

C   --Write 3D

call expnp (exoid, 20, 5, 5, ierr)
call expns (exoid, 20, node_list, ierr)
call expnsd (exoid, 20, dist_fact, ierr)
      IF (ANYNPS) THEN
C     ... Output nodeset id, number nodes, number dist factors
C     Assumes that there are the same number of distribution factors
C     as there are nodes in the nodeset.
         DO 10 ins = 1, numnps
           call expnp (ndbout, idnps(ins), nnnp3(ins), nnnp3(ins), ierr)
           if (ierr .lt. 0) then
              call exerr('gen3d2', 'Error from expnp', exlmsg)
              go to 50
           endif
           call expns (ndbout, idnps(ins), LTNNP3(IXNNP3(ins)), ierr)
           if (ierr .lt. 0) then
              call exerr('gen3d2', 'Error from expns', exlmsg)
              go to 50
           endif
           call expnsd(ndbout, idnps(ins), FACNP3(IXNNP3(ins)), ierr)
           if (ierr .lt. 0) then
              call exerr('gen3d2', 'Error from expnsd', exlmsg)
              go to 50
           endif
 10      continue
C     ... Output front and back nodesets (if any)
C     Front and back nodesets contain NUMNP (2D database) nodes
C     If there are any front or back, then create a temporary
C     Array to hold the distribution factors. Defaulted to 1.0
         if (nfro .gt. 0 .or. nbck .gt. 0) then
            call mdrsrv('factorns', knfac, numnp)
            call mdstat(mnerrs, mnused)
            if (mnerrs .gt. 0) goto 50
            call inirea(numnp, 1.0, a(knfac))

            do 20 ins = 1, nfro
              call expnp (ndbout, idfro(ins), numnp, numnp, ierr)
              if (ierr .lt. 0) then
                 call exerr('gen3d2', 'Error from expnp', exlmsg)
                 go to 50
              endif
              call expns (ndbout, idfro(ins), IXNP, ierr)
              if (ierr .lt. 0) then
                 call exerr('gen3d2', 'Error from expns', exlmsg)
                 go to 50
              endif
              call expnsd(ndbout, idfro(ins), a(knfac), ierr)
              if (ierr .lt. 0) then
                 call exerr('gen3d2', 'Error from expnsd', exlmsg)
                 go to 50
              endif
 20         continue
            if (nbck .gt. 0) then
               call mdrsrv('nodelist', knlst, numnp)
               call mdstat(mnerrs, mnused)
               if (mnerrs .gt. 0) goto 50

               do 30 i=1, numnp
                 ia(knlst+i-1) = ixnp(i) + nrnp(i) - 1
 30            continue

               do 40 ins = 1, nbck
                  call expnp (ndbout, idbck(ins), numnp, numnp, ierr)
                  if (ierr .lt. 0) then
                     call exerr('gen3d2', 'Error from expnp', exlmsg)
                     go to 50
                  endif
                  call expns (ndbout, idbck(ins), ia(knlst), ierr)
                  if (ierr .lt. 0) then
                     call exerr('gen3d2', 'Error from expns', exlmsg)
                     go to 50
                  endif
                  call expnsd(ndbout, idbck(ins), a(knfac), ierr)
                  if (ierr .lt. 0) then
                     call exerr('gen3d2', 'Error from expnsd', exlmsg)
                     go to 50
                  endif
 40            continue
            end if
         end if
      end if

      if (nfro .gt. 0 .or. nbck .gt. 0) then
         call mddel('factorns')
         if (nbck .gt. 0) then
            call mddel('nodelist')
         end if
      end if
      call mdstat(mnerrs, mnused)
      if (mnerrs .gt. 0) goto 50
      RETURN

 50   continue
      RETURN 1
      END
