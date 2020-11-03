C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRNPS (A, IA, IDFRO, IDBCK,
     &  IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS, *)
C=======================================================================

C   --*** WRNPS *** (GENSHELL) Write 3D node sets
C   --
C   --WRNPS writes the node set information for the 3D database.
C   --Calculations have been done elsewhere.
C   --
C   --Parameters:
C   --   IDFRO - IN - ids for front surface node sets; (0) = length
C   --   IDBCK - IN - ids for back surface node sets; (0) = length
C   --   IDNPS - IN - the 2D node sets ids
C   --   NNNPS - IN - the number of nodes for each 3D set
C   --   IXNNPS - IN - the index of the first node for each 3D set
C   --   LTNNPS - IN - the nodes for all 3D sets
C   --   FACNPS - IN - the distribution factors for all 3D sets
C   --
C   --Common Variables:
C   --   Uses NDBOUT of /DBASE/
C   --   Uses NUMNPS, LNPSNL of /DBNUMS/
C   --
C   --Database must be positioned at start of node set information
C   --upon entry; upon exit at end of node set information.

      INCLUDE 'exodusII.inc'
      INCLUDE 'gs_dbase.blk'
      INCLUDE 'gs_dbnums.blk'

      REAL A(*)
      INTEGER IA(*)
      INTEGER IDFRO(0:*)
      INTEGER IDBCK(0:*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)

      LOGICAL ANYNPS

      NFRO = IDFRO(0)
      NBCK = IDBCK(0)
      ANYNPS = (NFRO .GT. 0) .OR. (NBCK .GT. 0) .OR. (NUMNPS .GT. 0)

C   --Write 3D

      IF (ANYNPS) THEN
C     ... Output nodeset id, number nodes, number dist factors
C     Assumes that there are the same number of distribution factors
C     as there are nodes in the nodeset.
        DO 10 ins = 1, numnps
          call expnp (ndbout, idnps(ins), nnnps(ins), nnnps(ins), ierr)
          call expns (ndbout, idnps(ins), LTNNPS(IXNNPS(ins)), ierr)
          call expnsd(ndbout, idnps(ins), FACNPS(IXNNPS(ins)), ierr)
 10     continue
C     ... Output front and back nodesets (if any)
C     Front and back nodesets contain NUMNP (2D database) nodes
C     If there are any front or back, then create a temporary
C     Array to hold the distribution factors. Defaulted to 1.0
        if (nfro .gt. 0 .or. nbck .gt. 0) then
          call mdrsrv('factorns', knfac, numnp)
          call mdrsrv('nodlst', knl, numnp)
          call mdstat(mnerrs, mnused)
          if (mnerrs .gt. 0) goto 50
          call inirea(numnp, 1.0, a(knfac))
          do 20 i=1, numnp
            ia(knl+i-1) = i
 20       continue

          do 30 ins = 1, nfro
            call expnp (ndbout, idfro(ins), numnp, numnp, ierr)
            call expns (ndbout, idfro(ins), ia(knl), ierr)
            call expnsd(ndbout, idfro(ins), a(knfac), ierr)
 30       continue

          do 40 ins = 1, nbck
            call expnp (ndbout, idbck(ins), numnp, numnp, ierr)
            call expns (ndbout, idbck(ins), ia(knl), ierr)
            call expnsd(ndbout, idbck(ins), a(knfac), ierr)
 40       continue
          call mddel('factorns')
          call mddel('nodlst')
        end if
      end if

      call mdstat(mnerrs, mnused)
      if (mnerrs .gt. 0) goto 50

      RETURN

 50   continue
      return 1
      END

