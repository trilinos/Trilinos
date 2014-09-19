C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
C    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

C $Id: wrgens.f,v 1.9 1998/03/23 05:17:54 gdsjaar Exp $
C $Log: wrgens.f,v $
C Revision 1.9  1998/03/23 05:17:54  gdsjaar
C Fixed data statement ordering
C
C Revision 1.8  1997/09/02 18:51:45  caforsy
C Increase inital value of NIN in rdfsq.f and IGUESS in wrgens.f
C
C Revision 1.7  1992/11/17 15:06:39  gdsjaar
C Fixed problem with version number output.
C
c Revision 1.6  1992/11/16  23:32:47  gdsjaar
c Fixed problem with outputting too many attributes for some elements,
c upped version number to 2.2X
c
c Revision 1.5  1992/06/09  22:31:18  gdsjaar
c Fixed problem with sideset node and element numbering for 8 and 9 node elements
c
c Revision 1.4  1992/04/15  22:56:38  gdsjaar
c Fixed EXODUS output of 3-node beams--connectivity
c conforms to EXODUS document
c Upped version number to 2.1X to reflect change
c
c Revision 1.3  1991/12/18  21:12:03  gdsjaar
c Removed the info field
c
c Revision 1.2  1990/11/30  11:30:08  gdsjaar
c Rewrote indexing for reads and writes
c
c Revision 1.1.1.1  90/11/30  11:17:52  gdsjaar
c FASTQ Version 2.0X
c 
c Revision 1.1  90/11/30  11:17:51  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]WRGENS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE WRGENS (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &   NPNBC, NPSBC, IUNIT, NNN, KKK, NNXK, NODES, NELEMS, NNFLG,
     &   NNPTR, NNLEN, NSFLG, NSPTR, NSLEN, NVPTR, NVLEN, NSIDEN,
     &   MAPDXG, XN, YN, NXK, MAT, MAPGXD, MATMAP, WTNODE, WTSIDE,
     &   NBCNOD, NNLIST, NBCSID, NSLIST, NVLIST, NUMMAT, LINKM, TITLE,
     &   ERR, EIGHT, NINE, VERSN)
C************************************************************************
C
C  SUBROUTINE WRGENS = WRITES GENESIS DATABASE MESH OUTPUT
C
C***********************************************************************
C
      PARAMETER (IGUESS = 1000)
C
C  IGUESS IS THE NUMBER OF ELEMENT BLOCKS,  FOR USE WITH THE ENAME
C  VARIABLE.  IF THIS VARIABLE IS NOT AS LARGE AS NUMMAT,  IT WILL NOT
C  RESULT IN A FATAL ERROR,  BUT SIMPLY A WARNING,  AND NO ELEMENT
C  NAMES WILL BE WRITTEN.
C
      DIMENSION XN (NPNODE), YN (NPNODE), NXK (NNXK, NPELEM)
      DIMENSION MAT (NPELEM)
      DIMENSION NODES (NPNBC), NELEMS (NPSBC), NSIDEN (NPSBC)
      DIMENSION NNFLG (MXNFLG), NNLEN (MXNFLG)
      DIMENSION NNPTR (MXNFLG), WTNODE (NPNBC)
      DIMENSION NSFLG (MXSFLG), NSLEN (MXSFLG)
      DIMENSION NSPTR (MXSFLG), WTSIDE (NPSBC)
      DIMENSION NVLEN (MXSFLG), NVPTR (MXSFLG), LINKM (2,  (MS+MR))
      DIMENSION MAPDXG (NPNODE), MAPGXD (NPNODE), MATMAP (3, NPREGN)
C
      CHARACTER*72 TITLE, HOLD*80
      CHARACTER*8 DATE, TIME, VERSN1, VERSN2, XNAME, YNAME
      CHARACTER*8 ENAME (IGUESS)
      CHARACTER*10 VERSN
C
      LOGICAL ERR, EIGHT, NINE
C
      integer lcon(9)
      integer lbar(3)

      data lcon /1,3,5,7,2,4,6,8,9/
      data lbar /1,3,2/

      CALL EXDATE (DATE)
      CALL EXTIME (TIME)
      VERSN1 = '        '
      VERSN2 = '        '
      VERSN1 = VERSN (1:5)
      VERSN2 = VERSN (6:10)
      ERR = .TRUE.
      HOLD = TITLE
      XNAME = 'X'
      YNAME = 'Y'
C
C  CHECK TO MAKE SURE THAT THERE IS ENOUGH ROOM FOR ELEMENT NAMES
C
      IF (NUMMAT.GT.IGUESS) THEN
         CALL MESAGE ('WARNING:  THE NUMBER OF ELEMENT BLOCKS EXCEEDS')
         CALL MESAGE ('          THE CAPACITY TO NAME EACH BLOCK.')
         CALL MESAGE ('          NO ELEMENT NAMES WILL BE WRITTEN.')
      ENDIF
C
C  WRITE OUT HEADER INFORMATION
C
      WRITE (IUNIT, ERR = 110)HOLD
      IJK = 2
      IVERS = 1
C ... Fix up side set nodes for 8 and 9 node elements.
C ... At this point, they are treated as two linear segments, 
C ... They should be ends followed by middle

C    1-----3-----2 Now: 1 3 3 2 Correct: 1 2 3

      if (eight .or. nine) then
         nvlst = nvlist / 4 * 3
         if (nslist .gt. 1) then
            nslst = nslist / 2
         else
            nslst = nslist
         end if
      else
         nvlst = nvlist
         nslst = nslist
      end if

      WRITE (IUNIT, ERR = 110)NNN, IJK, KKK, NUMMAT, NBCNOD, NNLIST,
     &   NBCSID, NSLST, NVLST, IVERS
C
C  WRITE OUT NODE BLOCK
C
      WRITE (IUNIT, ERR = 110) (XN (I), I = 1, NNN),
     &   (YN (I), I = 1, NNN)
      WRITE (IUNIT, ERR = 110) (MAPDXG (I), I = 1, KKK)
C
C  WRITE OUT ELEMENT BLOCKS
C
      DO 100 I = 1, NUMMAT
         IF (NXK (3, MATMAP (2, I)) .EQ. 0) THEN
            INODE = 2
            NATTR = 1
            ATTR = 1.
            IF (I.LE.IGUESS)ENAME (I) = 'BEAM'
         ELSEIF (NXK (4, MATMAP (2, I)) .EQ. 0) THEN
            INODE = 3
            NATTR = 1
            ATTR = 1.
            IF (I.LE.IGUESS)ENAME (I) = 'BEAM'
            CALL MESAGE ('NOTE:  The connectivity numbering for 3-node')
            CALL MESAGE ('       beams/trusses has been fixed to')
            CALL MESAGE ('       conform to EXODUS convention (1-3-2)')
         ELSEIF (EIGHT) THEN
            INODE = 8
            NATTR = 0
            IF (I.LE.IGUESS)ENAME (I) = 'QUAD'
         ELSEIF (NINE) THEN
            INODE = 9
            NATTR = 0
            IF (I.LE.IGUESS)ENAME (I) = 'QUAD'
         ELSE
            INODE = 4
            NATTR = 0
            IF (I.LE.IGUESS)ENAME (I) = 'QUAD'
         ENDIF
C
C  NLOOP IS NEEDED TO WRITE SOMETHING OUT THE CURRENT COUNTER IS ZERO.
C  THIS IS DONE TO SOLVE A CRAY OPERATING SYSTEM [CTSS] PROBLEM
C  WHERE NULL RECORD WRITES ARE NOT DONE APPROPRIATELY
C
         WRITE (IUNIT, ERR = 110) MATMAP (1, I),
     &      MATMAP (3, I) - MATMAP (2, I)+1, INODE, NATTR
C... 8 or 9 node quads
         IF (INODE .EQ. 8 .or. inode .eq. 9) THEN
            write (iunit, err = 110) ((nxk(lcon(ii), k), ii=1, nnxk),
     $           k = matmap(2,i), matmap(3,i))
C... 3 node beam/truss
         ELSEIF (INODE .EQ. 3) THEN
            write (iunit, err = 110) ((nxk(lbar(ii), k), ii=1, inode),
     $           k = matmap(2,i), matmap(3,i))
C... 4 node quad or 2 node beam/truss
         ELSE
            WRITE (IUNIT, ERR = 110) ((NXK (J, K), J = 1, INODE),
     &         K = MATMAP (2, I), MATMAP (3, I))
         ENDIF
         NLOOP = MAX0 (1, NATTR*KKK)
         WRITE (IUNIT, ERR = 110) (ATTR, J = 1, NLOOP)
  100 CONTINUE
C
C  WRITE OUT NODAL BOUNDARY FLAGS
C
      NLOOP = MAX0 (1, NBCNOD)
      WRITE (IUNIT, ERR = 110) (NNFLG (I), I = 1, NLOOP)
      WRITE (IUNIT, ERR = 110) (NNLEN (I), I = 1, NLOOP)
      WRITE (IUNIT, ERR = 110) (NNPTR (I), I = 1, NLOOP)
      NLOOP = MAX0 (1, NNLIST)
      WRITE (IUNIT, ERR = 110) (NODES (I), I = 1, NLOOP)
      WRITE (IUNIT, ERR = 110) (WTNODE (I), I = 1, NLOOP)
C
C  WRITE OUT SIDE BOUNDARY FLAGS
C

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

      NLOOP = MAX0 (1, NBCSID)
      WRITE (IUNIT, ERR = 110) (NSFLG (I), I = 1, NLOOP)
      WRITE (IUNIT, ERR = 110) (NSLEN (I), I = 1, NLOOP)
      WRITE (IUNIT, ERR = 110) (NVLEN (I), I = 1, NLOOP)
      WRITE (IUNIT, ERR = 110) (NSPTR (I), I = 1, NLOOP)
      WRITE (IUNIT, ERR = 110) (NVPTR (I), I = 1, NLOOP)
      NLOOP = MAX0 (1, NSLIST)
      WRITE (IUNIT, ERR = 110) (NELEMS (I), I = 1, NLOOP)
      NLOOP = MAX0 (1, NVLIST)
      WRITE (IUNIT, ERR = 110) (NSIDEN (I), I = 1, NLOOP)
      WRITE (IUNIT, ERR = 110) (WTSIDE (I), I = 1, NLOOP)
C
C  WRITE OUT THE QA INFORMATION
C
      IHOLD = 1
      WRITE (IUNIT, ERR = 110)IHOLD
      WRITE (IUNIT, ERR = 110)VERSN1, VERSN2, DATE, TIME
C
C  WRITE THE HEADER INFORMATION
C
      IHOLD = 0
      WRITE (IUNIT, ERR = 110)IHOLD
C
C  WRITE THE COORDINATE NAMES AND ELEMENT NAMES
C
      WRITE (IUNIT)XNAME, YNAME
      IF (NUMMAT.LE.IGUESS)WRITE (IUNIT) (ENAME (I), I = 1, NUMMAT)
C
C  SUCCESSFULL WRITE COMPLETED
C
      CALL MESAGE (' ')
      CALL MESAGE (' ')
      CALL MESAGE ('GENESIS OUTPUT FILE SUCCESSFULLY WRITTEN')
      CALL MESAGE (' ')
      ERR = .FALSE.
      RETURN
C
C  ERR DURING WRITE PROBLEMS
C
  110 CONTINUE
      CALL MESAGE ('ERR DURING WRITE TO OUTPUT FILE')
      CALL MESAGE ('      - NO FILE SAVED -')
      RETURN
C
      END
