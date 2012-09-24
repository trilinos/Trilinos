C Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of Sandia Corporation nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
C                                                 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C=======================================================================
      SUBROUTINE DBPELB (OPTION, NELBLK, IDELB, NUMELB, NUMLNK, NUMATR,
     &   BLKTYP, EBNAME, ATNAME, NVAREL, NAMEEV, ISEVOK, LISEV)
C=======================================================================

C   --*** DBPELB *** (EXOLIB) Print database element block summary
C   --
C   --DBPELB displays the database element block summary (block ID,
C   --number of elements, element block name).
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      ' ' to print block summary
C   --      'N' to print element block name
C   --      'V' to print element variable truth table
C   --   NELBLK - IN - the number of element blocks
C   --   IDELB - IN - the element block ID for each block
C   --   NUMELB - IN - the number of elements for each block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   NUMATR - IN - the number of attributes for each block
C   --   BLKTYP - IN - the names of the element block types
C   --   NVAREL - IN - the number of element variables
C   --   NAMEEV - IN - the names of the element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   LISEV - SCRATCH - size = NVAREL (if 'V' in OPTION)

      include 'gp_params.blk'
      include 'gp_namlen.blk'
      CHARACTER*(*) OPTION
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      CHARACTER*(MXSTLN) BLKTYP(*)
      CHARACTER*(maxnam) EBNAME(*)
      CHARACTER*(maxnam) NAMEEV(*)
      CHARACTER*(maxnam) ATNAME(*)
      LOGICAL ISEVOK(*)
      INTEGER LISEV(*)

      LOGICAL ISABRT
      LOGICAL DONAM, DOVTBL
      CHARACTER*20 STRA

      DONAM  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0))
      DOVTBL = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0))

      WRITE (*, *)

      WRITE (STRA, 10000, IOSTAT=IDUM) NELBLK
10000  FORMAT ('(#', I5, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)

      isatn = 1
      DO 110 IELB = 1, NELBLK
         IF (ISABRT ()) RETURN

         WRITE (STRA, 10000, IOSTAT=IDUM) IELB
         CALL PCKSTR (1, STRA)
         WRITE (*, 10010, IOSTAT=IDUM) IDELB(IELB), STRA(:LSTRA),
     &      NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB)

         IF (DONAM) THEN
            WRITE (*, 10020) EBNAME(IELB)(:LENSTR(EBNAME(IELB))),
     $           BLKTYP(IELB)(:LENSTR(BLKTYP(IELB)))
         END IF

         if (numatr(ielb) .gt. 0) then
           mxnam = 0
           DO I = 1, numatr(ielb)
             LNAM = lenstr(atname(isatn+i-1))
             mxnam = max(lnam, mxnam)
           END DO
           if (mxnam .gt. 1) then
             WRITE (*, 10090, IOSTAT=IDUM)
     &         (ATNAME(I)(:mxnam),
     *         I=ISATN,ISATN+NUMATR(IELB)-1)
           end if
           isatn = isatn + numatr(ielb)
         end if
         
         IF (DOVTBL) THEN
            NSEL = 0
            DO 100 I = 1, NVAREL
               IF (ISEVOK( (IELB-1)*NELBLK+I )) THEN
                  NSEL = NSEL + 1
                  LISEV(NSEL) = I
               END IF
  100       CONTINUE

            WRITE (*, 10030, IOSTAT=IDUM) (NAMEEV(LISEV(I)), I=1,NSEL)
         END IF
  110 CONTINUE

      RETURN

10010  FORMAT (1X, 'Block', I11, 1X, A, ':',
     &   I10, ' elements',
     &   I10, '-node', I8, ' attributes')
10020  FORMAT (4X, 'Element block name = "',A,
     $      '", Element type = "', A, '"')
10030  FORMAT (4X, 'Defined variables:', :, 3X, 4 (2X, A), :, /
     &   (4X, 1X, 6 (2X, A)))
10090 FORMAT (4X, 'Attributes: ', 10(2X, A))
      END
