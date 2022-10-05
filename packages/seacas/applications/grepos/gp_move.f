C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MOVE(ndbin, A, IA, isnp, X, Y, Z, NDIM, numnp, numel,
     *  numess, idess, neess, ixeess, lteess,
     *  nelb, idelb, numelb, numlnk, link,
     *  issblk, iscrn, iscre)
C=======================================================================

      REAL A(*)
      INTEGER IA(*)
      REAL X(*), Y(*), Z(*)
      INTEGER NUMESS
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER IXEESS(*)
      INTEGER LTEESS(*)
      INTEGER IDELB(NELB)
      INTEGER NUMELB(NELB)
      INTEGER NUMLNK(NELB)
      INTEGER LINK(*)
      INTEGER ISSBLK(NELB), ISCRN(*), ISCRE(*)

      INCLUDE 'gp_snap.blk'

C   --   X,Y,Z -- REAL - IN/OUT - Coordinates of nodes
C   --   NDIM   - IN - the number of dimensions
C   --   NUMESS - IN - the number of side sets
C   --   IDESS  - IN - the side set ID for each set
C   --   NEESS  - IN - the number of elements for each set
C   --   IXEESS - IN - the index of the first element for each set
C   --   LTEESS - IN - the elements for all sets

      call iniint(numel, 0, iscre)
      call iniint(numnp, 0, iscrn)
      call iniint(nelb,  0, issblk)

      call mdrsrv('ISCR',  KISCR,  NUMNP)
      call mdrsrv('VNORM', KVNORM, NUMNP*3)
      call mdrsrv('PLANE', KPLANE, 0)

      indma = 0
      indsl = 0

C ... Find master index
      do 100 i=1, numess
        if (idess(i) .eq. idssma(isnp)) then
          indma = i
          go to 110
        endif
 100  continue
 110  continue

C ... Find slave index
      do 120 i=1, numess
        if (idess(i) .eq. idsssl(isnp)) then
          indsl = i
          go to 130
        endif
 120  continue
 130  continue

C------------------------------------------------------------------------
C ... Get the number of nodes on the slave surface
      call exgsp(ndbin, idess(indsl), nssess, nsdess, ierr)
      if (nsdess .eq. 0) then
C...     Distribution factors not stored, estimate max size of node list
C        based on maximum of 9 nodes/face.
         nsdess = nssess * 9
      end if
C...  Allocate storage for sideset node and count list.
      call mdrsrv('NODSLV', islvnd, nsdess)
      call mdrsrv('NNDSLV', islvnn, nssess)
C...  Get the sideset nodes...
      call exgssn(ndbin, idess(indsl), ia(islvnn), ia(islvnd), ierr)

C ... Get the number of nodes on the master surface
      call exgsp(ndbin, idess(indma), nmsess, nmdess, ierr)
      if (nmdess .eq. 0) then
C...     Distribution factors not stored, estimate max size of node list
C        based on maximum of 9 nodes/face.
         nmdess = nmsess * 9
      end if
C...  Allocate storage for sideset node and count list.
      call mdrsrv('NODMAS', imasnd, nmdess)
      call mdrsrv('NNDMAS', imasnn, nmsess)
C...  Get the sideset nodes...
      call exgssn(ndbin, idess(indma), ia(imasnn), ia(imasnd), ierr)
C------------------------------------------------------------------------

C ... Find element blocks containing elements in sideset 'indsl'
C ... First, use the scratch element storage to map element to block
      ibeg = 1
      do 140 iblk = 1, nelb
        call iniint(numelb(iblk), iblk, iscre(ibeg))
        ibeg = ibeg+numelb(iblk)
 140  continue

C ... Now, for each element in the sideset, find which block it is in
C     and set 'issblk' equal to 1.
      ibeg = ixeess(indsl)
      do 150 isse = ibeg, ibeg+neess(indsl)-1
        iel = lteess(isse)
        if (iscre(iel) .eq. 1) then
          PRINT *, iel
        end if
        issblk(iscre(iel)) = 1
 150  continue

C ... Now, for each block in the sideset, find the nodes contained in
C     that block.
      ielnk = 0
      do 170 iblk = 1, nelb
        ISLNK = IELNK + 1
        IELNK = IELNK + NUMLNK(IBLK) * NUMELB(IBLK)
        if (issblk(iblk) .eq. 1) then
          call exfcon(1, numelb(iblk), numlnk(iblk), link(islnk),
     *      iscrn)
        end if
 170  continue
C ... At this point, iscrn(i) will be 1 if node in block, 0 otherwise

C ... Allocate temporary storage
      nummaf = neess(indma)
      call mdlong('PLANE', KPLANE, nummaf*4)

      call inirea(3*numnp,  0.0, a(kvnorm))
      call inirea(4*nummaf, 0.0, a(kplane))

C ... Calculate the slave surface normals
      if (ndim .eq. 3) then
        call setnor(.FALSE., numnp, x, y, z, 3, a(kvnorm),
     *    usnorm(isnp), vector(1, isnp),
     *    nssess, nsdess, ia(islvnd))

C ... Calculate the master surface planes
        call setnor(.TRUE., numnp, x, y, z, 4, a(kplane),
     *    0, vector(1, isnp),
     *    nmsess, nmdess, ia(imasnd))

        call movnod(numnp, ndim, x, y, z,
     *    a(kvnorm), nssess, nsdess, ia(islvnd),
     *    a(kplane), nmsess, nmdess, ia(imasnd),
     *    snptol(isnp), delmax(isnp), iscrn, vector(1,isnp),
     *    gap(isnp))
      else
        call setnr2 (.FALSE., numnp, x, y, 3, a(kvnorm),
     *    usnorm(isnp), vector(1, isnp),
     *    nssess, nsdess, ia(islvnd))

        call movnd2(numnp, ndim, x, y, a(kvnorm),
     *    nmsess, nmdess, ia(imasnd),
     *    snptol(isnp), delmax(isnp), iscrn, vector(1,isnp),
     *    gap(isnp))
      end if

      call mddel('NNDMAS')
      call mddel('NODMAS')
      call mddel('NNDSLV')
      call mddel('NODSLV')

      call mddel('PLANE')
      call mddel('VNORM')
      call mddel('ISCR')

      RETURN
      END

