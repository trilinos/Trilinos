C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
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

C=======================================================================
      SUBROUTINE DBPELB (OPTION, NELBLK, IDELB, NUMELB, NUMLNK, NUMATR,
     &   NAMELB, NVAREL, NAMEEV, ISEVOK, LISEV)
C=======================================================================
C$Id: dbpelb.f,v 1.5 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbpelb.f,v $
CRevision 1.5  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.4  2009/03/12 13:27:06  gdsjaar
CWiden output format
C
CRevision 1.3  1998/07/13 21:26:33  gdsjaar
CModified format statement to fit on line better
C
CRevision 1.2  1997/03/20 19:40:13  caforsy
CUpdated Imakefile for Imake 6.1.  Changed printing routines to handle
Clarger problems.
C
CRevision 1.1.1.1  1990/08/14 16:13:45  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:13:44  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:18  gdsjaar
c Initial revision
c 

C   --*** DBPELB *** (EXOLIB) Print database element block summary
C   --   Written by Amy Gilkey - revised 01/19/88
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
C   --   NAMELB - IN - the names of the element block types
C   --   NVAREL - IN - the number of element variables
C   --   NAMEEV - IN - the names of the element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   LISEV - SCRATCH - size = NVAREL (if 'V' in OPTION)

      CHARACTER*(*) OPTION
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      CHARACTER*8 NAMELB(*)
      CHARACTER*8 NAMEEV(*)
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

      DO 110 IELB = 1, NELBLK
         IF (ISABRT ()) RETURN

         WRITE (STRA, 10000, IOSTAT=IDUM) IELB
         CALL PCKSTR (1, STRA)
         WRITE (*, 10010, IOSTAT=IDUM) IDELB(IELB), STRA(:LSTRA),
     &      NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB)

         IF (DONAM) THEN
            WRITE (*, 10020) NAMELB(IELB)
         END IF

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

10010  FORMAT (1X, 'Block', I10, 1X, A, ':',
     &   I10, ' elements',
     &   I6, '-node', I4, ' attributes')
10020  FORMAT (4X, 'Element block type = "', A, '"')
10030  FORMAT (4X, 'Defined variables:', :, 3X, 4 (2X, A), :, /
     &   (4X, 1X, 6 (2X, A)))
      END
