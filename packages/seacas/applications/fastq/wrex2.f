C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WREX2 (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &   NPNBC, NPSBC, IUNIT, NNN, KKK, NNXK, NODES, NELEMS, NNFLG,
     &   NNPTR, NNLEN, NSFLG, NSPTR, NSLEN, NVPTR, NVLEN, NSIDEN,
     &   MAPDXG, XN, YN, NXK, MAT, MAPGXD, MATMAP, WTNODE, WTSIDE,
     &   NBCNOD, NNLIST, NBCSID, NSLIST, NVLIST, NUMMAT, LINKM, TITLE,
     &   ERR, EIGHT, NINE, VERSN, A, IA, FILENAME)
C=======================================================================
C************************************************************************

C  SUBROUTINE WREX2 = WRITES GENESIS DATABASE MESH OUTPUT

C***********************************************************************

      include 'exodusII.inc'

      DIMENSION XN (NPNODE), YN (NPNODE), NXK (NNXK, NPELEM)
      DIMENSION MAT (NPELEM)
      DIMENSION NODES (NPNBC), NELEMS (NPSBC), NSIDEN (NPSBC)
      DIMENSION NNFLG (MXNFLG), NNLEN (MXNFLG)
      DIMENSION NNPTR (MXNFLG), WTNODE (NPNBC)
      DIMENSION NSFLG (MXSFLG), NSLEN (MXSFLG)
      DIMENSION NSPTR (MXSFLG), WTSIDE (NPSBC)
      DIMENSION NVLEN (MXSFLG), NVPTR (MXSFLG), LINKM (2,  (MS+MR))
      DIMENSION MAPDXG (NPNODE), MAPGXD (NPNODE), MATMAP (3, NPREGN)

      CHARACTER*(*) FILENAME
      CHARACTER*1 cdum
      CHARACTER*72 TITLE, HOLD*(MXLNLN)
      CHARACTER*(MXSTLN) QAREC(4), CNAME(2)
      CHARACTER*(MXSTLN) ENAME
      CHARACTER*10 VERSN

      REAL A(*)
      INTEGER IA(*)

      LOGICAL ERR, EIGHT, NINE

      integer lcon(9)
      integer lbar(3)
      integer lquad(4)

      data lcon /1,3,5,7,2,4,6,8,9/
      data lbar /1,3,2/
      data lquad /1,2,3,4/

      QAREC(1) = '        '
      QAREC(2) = '        '
      QAREC(3) = '        '
      QAREC(4) = '        '
      CALL EXDATE (QAREC(3))
      CALL EXTIME (QAREC(4))
      QAREC(1) = VERSN (1:5)
      QAREC(2) = VERSN (6:10)
      call expqa (iunit, 1, QAREC, ierr)

      ERR = .TRUE.
      HOLD = TITLE
      CNAME(1) = 'X'
      CNAME(2) = 'Y'

C  WRITE OUT HEADER INFORMATION

      IJK = 2
      call expini (iunit, hold, ijk, nnn, kkk, nummat, nbcnod,
     &     nbcsid, ierr)

C ... WRITE OUT NODE BLOCK
      call expcor (iunit, xn, yn, rdum, ierr)
      call expcon (iunit, cname, ierr)

C ... Write out the element map
      call expmap (iunit, mapdxg, ierr)

C ... Write out element blocks
C ... The fastq connectivity storage is not consistent with the
C     exodusII connectivity storage order. Determine the maximum
C     storage required to store the connectivity for an element
C     block and allocate that space. We duplicate some code here,
C     but we don't have any extra arrays to store the information
      MAXLNK = 0
      MAXATT = 0
      DO 90 I = 1, NUMMAT
        NUMEL = MATMAP (3, I) - MATMAP (2, I)+1
        IF (NXK (3, MATMAP (2, I)) .EQ. 0) THEN
C ...2-NODE BEAM
           INODE = 2
           NATTR = 1
        ELSEIF (NXK (4, MATMAP (2, I)) .EQ. 0) THEN
C ...3-NODE BEAM
           INODE = 3
           NATTR = 1
        ELSEIF (EIGHT) THEN
C ...8-NODE QUAD
           INODE = 8
           NATTR = 0
        ELSEIF (NINE) THEN
C ...9-NODE QUAD
           INODE = 9
           NATTR = 0
        ELSE
C ...4-NODE QUAD
           INODE = 4
           NATTR = 0
        ENDIF
        maxlnk = max(maxlnk, numel * inode)
        maxatt = max(maxatt, numel * nattr)
 90   CONTINUE

C ... Have the maximum size of the link and attribute arrays,
C     Now allocate the space.
      CALL MDRSRV('LINK', KLINK, maxlnk)
      CALL MDRSRV('ATTR', KATRIB, maxatt)
      call mdstat(nerr, mused)
      if (nerr .gt. 0) go to 110

      DO 100 I = 1, NUMMAT
         IF (NXK (3, MATMAP (2, I)) .EQ. 0) THEN
            INODE = 2
            NATTR = 1
            ATTR = 1.
            ENAME = 'BEAM'
         ELSEIF (NXK (4, MATMAP (2, I)) .EQ. 0) THEN
            INODE = 3
            NATTR = 1
            ATTR = 1.
            ENAME = 'BEAM3'
        ELSEIF (EIGHT) THEN
            INODE = 8
            NATTR = 0
            ENAME = 'QUAD8'
         ELSEIF (NINE) THEN
            INODE = 9
            NATTR = 0
            ENAME = 'QUAD9'
         ELSE
            INODE = 4
            NATTR = 0
            ENAME = 'QUAD'
         ENDIF

         call expelb(iunit, matmap(1,i), ename,
     &        MATMAP (3, I) - MATMAP (2, I)+1, INODE, NATTR, ierr)
C... 8 or 9 node quads
         IF (INODE .EQ. 8 .or. inode .eq. 9) THEN
            call trnlnk(ia(klink), nxk, nnxk, lcon, nnxk,
     &           matmap(2,i), matmap(3,i), .TRUE.)
            call expelc(iunit, matmap(1,i), ia(klink), ierr)
C... 3 node beam/truss
         ELSEIF (INODE .EQ. 3) THEN
            call trnlnk(ia(klink), nxk, nnxk, lbar, inode,
     &           matmap(2,i), matmap(3,i), .TRUE.)
            call expelc(iunit, matmap(1,i), ia(klink), ierr)
C... 4 node quad or 2 node beam/truss
         ELSE
            call trnlnk(ia(klink), nxk, nnxk, lquad, inode,
     &           matmap(2,i), matmap(3,i), .FALSE.)
            call expelc(iunit, matmap(1,i), ia(klink), ierr)
         ENDIF
         if (NATTR .gt. 0) then
C ... Initialize attributes
           NUMEL = MATMAP (3, I) - MATMAP (2, I)+1
           do 95 iat = 0, nattr*numel-1
             a(katrib+iat) = attr
 95        continue
           call expeat(iunit, matmap(1,i), a(katrib), ierr)
         end if
c         NLOOP = MAX0 (1, NATTR*KKK)
c         WRITE (IUNIT, ERR = 110) (ATTR, J = 1, NLOOP)
  100 CONTINUE
      call mddel('LINK')
      call mddel('ATTR')
      call mdstat(nerr, mused)
      if (nerr .gt. 0) go to 110

C  WRITE OUT NODAL BOUNDARY FLAGS
      do 200 i = 1, nbcnod
        call expnp (iunit, nnflg(i), nnlen(i), nnlen(i), ierr)
        call expns (iunit, nnflg(i), nodes(nnptr(i)), ierr)
        call expnsd(iunit, nnflg(i), wtnode(nnptr(i)), ierr)
 200  continue

C  WRITE OUT SIDE BOUNDARY FLAGS
C ... Fix up side set nodes and elements for 8 and 9 node elements.
C ... At this point, they are treated as two linear segments,
C ... They should be ends followed by middle

C    1-----3-----2 Now: 1 3 3 2 Correct: 1 2 3

      if (eight .or. nine) then
         ipn = 1
         ipe = 1
         do 130 ibc = 1, nbcsid

            ibee = nsptr(ibc)
            iene = nsptr(ibc) + nslen(ibc) - 1
            nsptr(ibc) = ipe
            nslen(ibc) = max(1, nslen(ibc)/2)
            do 115 iel = ibee, iene, 2
               nelems(ipe)   = nelems(iel)
               ipe = ipe + 1
 115        continue

            ibeg = nvptr(ibc)
            iend = nvptr(ibc) + nvlen(ibc) - 1
            nvptr(ibc) = ipn
            nvlen(ibc) = nvlen(ibc) / 4 * 3

            do 120 inod = ibeg, iend, 4
               nsiden(ipn)   = nsiden(inod)
               nsiden(ipn+2) = nsiden(inod+1)
               nsiden(ipn+1) = nsiden(inod+3)
               wtside(ipn)   = wtside(inod)
               wtside(ipn+2) = wtside(inod+1)
               wtside(ipn+1) = wtside(inod+3)
               ipn = ipn + 3
 120        continue
 130     continue
         nvlist = ipn - 1
         nslist = ipe - 1
      end if

C ... Convert sideset nodes to sideset sides
      if (nbcsid .gt. 0) then
         call mdrsrv('ISIDES', ksides, nslist)
         call EXCN2S(iunit, nslen, nvlen, nsptr, nvptr,
     &        nelems, nsiden, ia(ksides), IERR)

         do 140 i=1, nbcsid
           call expsp (iunit, nsflg(i), nslen(i), nvlen(i), ierr)
           call expss (iunit, nsflg(i), nelems(nsptr(i)),
     &          ia(ksides+nsptr(i)-1), ierr)
           call expssd(iunit, nsflg(i), wtside(nvptr(i)), ierr)
 140     continue
         call mddel('ISIDES')
      end if
C     SUCCESSFUL WRITE COMPLETED
      CALL MESAGE (' ')
      CALL MESAGE ('ExodusII output file successfully written')
C ... Title is char*72, dbpini and exodusII expect char*80
      HOLD = TITLE
      if (nbcnod .gt. 0) then
         call exinq (iunit, EXNSNL, lnsnl, rdum, cdum, ierr)
         call exinq (iunit, EXNSDF, lnsdf, rdum, cdum, ierr)
      else
         lnsnl = 0
         lnsdf = 0
      end if
      if (nbcsid .gt. 0) then
         call exinq (iunit, EXSSEL, lssel, rdum, cdum, ierr)
         call exinq (iunit, EXSSDF, lssdf, rdum, cdum, ierr)
         call exinq (iunit, EXSSNL, lssnl, rdum, cdum, ierr)
      else
         lssel = 0
         lssdf = 0
         lssnl = 0
      end if
      CALL FQDBPINI ('NTIS', HOLD, ijk, nnn, kkk, nummat, nbcnod,
     *  nnlist, nnlist, nbcsid, lssel, lssnl, lssdf, 0, 0, 0,
     *  filename)
      ERR = .FALSE.
      RETURN

C     ERROR DURING WRITE PROBLEMS
 110  CONTINUE
      CALL MESAGE ('ERR DURING WRITE TO OUTPUT FILE')
      CALL MESAGE ('...File may be incomplete...')
      RETURN

      END

      subroutine trnlnk(LINK, NXK, NNXK, INDEX, NNODE, IBEG, IEND,USING)
      integer link(nnode, *)
      integer nxk(nnxk, *)
      integer index(*)
      logical using

      ii = 0
      do 20 i = ibeg, iend
        ii = ii + 1
        do 10 j = 1, nnode
          if (using) then
             link(j,ii) = nxk(index(j), i)
          else
             link(j,ii) = nxk(j, i)
          end if
 10     continue
 20   continue
      return
      end
