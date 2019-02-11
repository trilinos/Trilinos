C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
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

C $Id: wrnast.f,v 1.1 1990/11/30 11:17:57 gdsjaar Exp $
C $Log: wrnast.f,v $
C Revision 1.1  1990/11/30 11:17:57  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]WRNAST.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE WRNAST (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &   NPNBC, NPSBC, IUNIT, NNN, KKK, NNXK, NODES, NELEMS, NNFLG,
     &   NNPTR, NNLEN, NSFLG, NSPTR, NSLEN, NVPTR, NVLEN, NSIDEN,
     &   MAPDXG, XN, YN, NXK, MAT, MAPGXD, MATMAP, NBCNOD, NNLIST,
     &   NBCSID, NSLIST, NVLIST, NUMMAT, LINKM, TITLE, ERR, EIGHT,
     &   NINE, LONG)
C************************************************************************
C
C  SUBROUTINE WRNAST = WRITES NASTRAN DATABASE MESH OUTPUT FILE
C
C***********************************************************************
C
      DIMENSION XN(NPNODE), YN(NPNODE), NXK(NNXK, NPELEM), MAT(NPELEM)
      DIMENSION NODES(NPNBC), NELEMS(NPSBC), NSIDEN(NPSBC)
      DIMENSION NNFLG(MXNFLG), NNLEN(MXNFLG), NNPTR(MXNFLG)
      DIMENSION NSFLG(MXSFLG), NSLEN(MXSFLG), NSPTR(MXSFLG)
      DIMENSION NVLEN(MXSFLG), NVPTR(MXSFLG), LINKM(2, (MS+MR))
      DIMENSION MAPDXG(NPNODE), MAPGXD(NPNODE), MATMAP(3, NPREGN)
C
      CHARACTER*72 TITLE, DUMMY, DUMMY2
C
      LOGICAL ERR, EIGHT, NINE, DEFTYP, LONG
C
      ERR = .TRUE.
C
C  WRITE OUT HEADER TITLE AND INFORMATION
C
      WRITE(IUNIT, 10000, ERR = 180)TITLE
      WRITE(IUNIT, 10010, ERR = 180)NNN, KKK, NBCNOD
C
C  WRITE OUT NODE BLOCK
C
      WRITE(IUNIT, 10020, ERR = 180)
      Z = 0.
      DO 100 I = 1, NNN
         CALL GETDUM(I, DUMMY, LEN)
         IF (LONG) THEN
            WRITE(IUNIT, 10030, ERR = 180)I, XN(I), YN(I), DUMMY(1:6),
     &         DUMMY(1:6), Z
         ELSE
            WRITE(IUNIT, 10040, ERR = 180)I, XN(I), YN(I), Z
         ENDIF
  100 CONTINUE
C
C  QUERY THE USER FOR LOCAL CONTROL OF ELEMENT TYPE
C
      CALL INQTRU('USE DEFAULT ELEMENT TYPES FOR ELEMENT BLOCKS',
     &   DEFTYP)
C
C  WRITE OUT ELEMENT BLOCKS
C
      DO 150 I = 1, NUMMAT
         CALL GETDUM(MATMAP(1, I), DUMMY, LEN)
         IF(NXK(3, MATMAP(2, I)).EQ.0)THEN
            WRITE(IUNIT, 10050, ERR = 180)DUMMY(1:LEN)
            INODE = 2
            IF(DEFTYP)THEN
               DUMMY2 = 'CBAR'
               LEN2 = 4
            ELSE
               WRITE(*, 10060)MATMAP(1, I)
               CALL INQSTR('2 NODE ELEMENT TYPE:  ', DUMMY2)
               CALL STRLNG(DUMMY2, LEN2)
            ENDIF
         ELSEIF (NXK(4, MATMAP(2, I)) .EQ. 0)THEN
            CALL MESAGE('THREE NODE BAR ELEMENTS NOT SUPPORTED')
            CALL MESAGE('THE CENTER NODE WILL BE IGNORED')
            WRITE(IUNIT, 10050, ERR = 180)DUMMY(1:LEN)
            IF (DEFTYP)THEN
               DUMMY2 = 'CBAR'
               LEN2 = 4
            ELSE
               WRITE(*, 10060) MATMAP(1, I)
               CALL INQSTR ('2 NODE ELEMENT TYPE:  ', DUMMY2)
               CALL STRLNG (DUMMY2, LEN2)
            ENDIF
            INODE = 3
         ELSEIF(EIGHT.OR.NINE)THEN
            WRITE(IUNIT, 10070, ERR = 180) DUMMY(1:LEN)
            IF (NINE) THEN
               CALL MESAGE('NINE NODE QUAD ELEMENTS NOT SUPPORTED')
               CALL MESAGE('THE CENTER NODE WILL BE IGNORED')
            ENDIF
            IF(DEFTYP)THEN
               DUMMY2 = 'CQUAD8'
               LEN2 = 6
            ELSE
               WRITE(*, 10060)MATMAP(1, I)
               CALL INQSTR('8 NODE ELEMENT TYPE:  ', DUMMY2)
               CALL STRLNG(DUMMY2, LEN2)
            ENDIF
            INODE = 8
         ELSE
            WRITE(IUNIT, 10080, ERR = 180)DUMMY(1:LEN)
            IF(DEFTYP)THEN
               DUMMY2 = 'CQUAD4'
               LEN2 = 6
            ELSE
               WRITE(*, 10060)MATMAP(1, I)
               CALL INQSTR('4 NODE ELEMENT TYPE:  ', DUMMY2)
               CALL STRLNG(DUMMY2, LEN2)
            ENDIF
            INODE = 4
         ENDIF
         CALL STRIPB(DUMMY2, ILEFT, IRIGHT)
         IRIGHT = ILEFT+7
         IF(NXK(3, MATMAP(2, I)).EQ.0)THEN
            DO 110 K = MATMAP(2, I), MATMAP(3, I)
               WRITE(IUNIT, 10090, ERR = 180)DUMMY2(ILEFT:IRIGHT), K,
     &            MATMAP(1, I), (NXK(J, K), J = 1, INODE)
  110       CONTINUE
         ELSEIF(NXK(4, MATMAP(2, I)).EQ.0)THEN
            DO 120 K = MATMAP(2, I), MATMAP(3, I)
               WRITE(IUNIT, 10090, ERR = 180)DUMMY2(ILEFT:IRIGHT), K,
     &            MATMAP(1, I), (NXK(J, K), J = 1, INODE, 2)
  120       CONTINUE
         ELSEIF(EIGHT.OR.NINE)THEN
            DO 130 K = MATMAP(2, I), MATMAP(3, I)
               WRITE(IUNIT, 10100, ERR = 180)DUMMY2(ILEFT:IRIGHT), K,
     &            MATMAP(1, I), NXK(1, K), NXK(3, K), NXK(5, K),
     &            NXK(7, K), NXK(2, K), NXK(4, K), NXK(6, K), NXK(8, K)
  130       CONTINUE
         ELSE
            DO 140 K = MATMAP(2, I), MATMAP(3, I)
               WRITE(IUNIT, 10110, ERR = 180)DUMMY2(ILEFT:IRIGHT), K,
     &            MATMAP(1, I), (NXK(J, K), J = 1, INODE)
  140       CONTINUE
         ENDIF
  150 CONTINUE
C
C  WRITE OUT THE NODAL BOUNDARY CONDITIONS
C
      IF(NBCNOD.GT.0)THEN
         DO 170 I = 1, NBCNOD
            J1 = NNPTR(I)
            J2 = NNPTR(I) + NNLEN(I)-1
            CALL GETDUM (NNFLG(I), DUMMY, LEN)
            WRITE (IUNIT, 10120) DUMMY(1:LEN)
            WRITE (*, 10130)NNFLG(I)
            CALL INQSTR ('DEGREES OF FREEDOM  RESTRAINED (NO BLANKS): ',
     &         DUMMY)
            DO 160 J = J1, J2
               WRITE(IUNIT, 10140, ERR = 180)NNFLG(I), NODES(J),
     &            DUMMY(1:6)
  160       CONTINUE
  170    CONTINUE
      ENDIF
C
C  NOTIFY USER THAT SIDE BOUNDARY FLAGS ARE NOT SUPPORTED
C
      IF (NBCSID .GT. 0) THEN
         CALL MESAGE('NO SIDE BOUNDARY FLAGS SUPPORTED BY NASTRAN')
      ENDIF
C
C  END THE DATA
C
      WRITE(IUNIT, 10150)
      CALL MESAGE ('NASTRAN OUTPUT FILE SUCCESSFULLY WRITTEN')
      ERR = .FALSE.
      RETURN
C
C  ERR DURING WRITE PROBLEMS
C
  180 CONTINUE
      CALL MESAGE ('ERR DURING WRITE TO ABAQUS OUTPUT FILE')
      CALL MESAGE ('         - NO FILE SAVED -            ')
      RETURN
C
C
10000 FORMAT('$TITLE: ', /, A72)
10010 FORMAT('$', /,
     &   '$     MESH GENERATED USING FASTQ        ', /,
     &   '$  NUMBER OF NODES:                     ', I5, /,
     &   '$  NUMBER OF ELEMENTS:                  ', I5, /,
     &   '$  NUMBER OF NODAL BOUNDARY CONDITIONS: ', I5, /,
     &   '$', /,
     &   'BEGIN BULK')
10020 FORMAT('$ NODE (GRID) DATA FOLLOWS:')
10030 FORMAT('GRID*   ', I16, 16X, 2E16.9, '*N', A6, /, '*N', A6,
     &   E16.9, '345')
10040 FORMAT('GRID    ', I8, 8X, 3F8.4, '345')
10050 FORMAT('$ 2 NODE BAR ELEMENTS FOR BLOCK ID ', A, ' FOLLOW:')
10060 FORMAT(' FOR BLOCK ID:', I7, '  ENTER NEW')
10070 FORMAT('$ 8 NODE QUAD ELEMENTS FOR BLOCK ID ', A, ' FOLLOW:')
10080 FORMAT('$ 4 NODE QUAD ELEMENTS FOR BLOCK ID ', A, ' FOLLOW:')
10090 FORMAT(A8, 4I8)
10100 FORMAT(A8, 9I8, /, 8X, I8)
10110 FORMAT(A8, 6I8)
10120 FORMAT('$ NODAL CONSTRAINTS FOR BOUNDARY FLAG ', A, ' FOLLOW:')
10130 FORMAT(' INPUT THE CONSTRAINTS FOR NODAL BOUNDARY FLAG: ', I5)
10140 FORMAT('SPC     ', 2I8, A8)
10150 FORMAT('ENDDATA')
C
      END
