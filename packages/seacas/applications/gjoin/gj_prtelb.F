C Copyright (c) 2008 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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
C 
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
C 

C=======================================================================
      SUBROUTINE PRTELB (ISTAT, NELBL1, NELBL2,
     &   IDELB, NUMELB, NUMLNK, NUMATR, ISCR)
C=======================================================================
C $Id: prtelb.F,v 1.2 1999/03/23 15:56:43 gdsjaar Exp $
C $Log: prtelb.F,v $
C Revision 1.2  1999/03/23 15:56:43  gdsjaar
C Upped format fields to accomodate larger element block ids
C
C Revision 1.1  1999/01/18 19:21:25  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:29  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1  1998/04/09 05:01:59  gdsjaar
C The code which prints out the sidesets/nodesets/blocks for id
C combining/renumbering needed the needsdoubleescape if block for
C certain systems. Needed to move those routines from .f to .F
C
C Revision 1.3  1997/04/09 22:10:07  gdsjaar
C Modified the node matching code to sort the array on the coordinate
C with the largest span. Previously it was done on the x coordinate
C always.
C
C Upped the format string used for the number of elements in the prtelb,
C prtess, and prtnps output to handle large models.
C
C Revision 1.2  1993/07/23 22:34:29  gdsjaar
C Cleaned up unused variables, removed unneeded routine
C
c Revision 1.1.1.1  1990/11/12  14:35:39  gdsjaar
c GJOIN - X1.00.40 - 7/17/90
c
c Revision 1.1  90/11/12  14:35:38  gdsjaar
c Initial revision
c 

C   --*** PRTELB *** (GJOIN) Print summary for element blocks
C   --   Written by Amy Gilkey - revised 02/04/88
C   --
C   --PRTELB prints a summary of the element blocks, including their
C   --status.
C   --
C   --Parameters:
C   --   ISTAT - IN - the status of each element block:
C   --      0 = same
C   --      - = delete
C   --      n = combine with block n
C   --   NELBL1 - IN - the number of element blocks in the first database
C   --   NELBL2 - IN - the number of element blocks in the second database
C   --   IDELB - IN - the element block IDs for each block
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   NUMATR - IN - the number of attributes in each block
C   --   ISCR - SCRATCH - size = NELBL1+NELBL2

      INTEGER ISTAT(*)
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      INTEGER ISCR(*)

      CHARACTER*2 STR2
      CHARACTER*20 STRA, STRB
      CHARACTER CH

      NELBLK = NELBL1 + NELBL2

      IF (NELBLK .LE. 0) RETURN

      WRITE (*, *)
      WRITE (*, *) 'Element blocks:'

      WRITE (STRA, 10000) NELBLK
10000  FORMAT ('(#', I5, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)

      DO 110 IELB = 1, NELBLK
         IF (ISTAT(IELB) .EQ. 0) THEN
            WRITE (STRA, 10000) IELB
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (IELB .GT. NELBL1) CH = '*'
            WRITE (*, 10010) CH, IDELB(IELB), STRA(:LSTRA),
     &         NUMELB(IELB), NUMLNK(IELB)
         ELSE IF (ISTAT(IELB) .GT. 0) THEN
            IF (ISTAT(IELB) .EQ. IELB) THEN
               CALL GETALL (IELB, NELBLK, ISTAT, NINSET, ISCR)
               WRITE (STRB, '(I10)') IDELB(IELB)
               CALL SQZSTR (STRB, LSTRB)
               DO 100 ISET = 1, NINSET
                  IF (ISET .EQ. 1) THEN
#if NeedsDoubleEscape
                     STR2 = '\\ '
#else
                     STR2 = '\ '
#endif
                  ELSE IF (ISET .EQ. NINSET) THEN
                     STR2 = '/ '
                  ELSE
                     STR2 = ' |'
                  END IF
                  I = ISCR(ISET)
                  WRITE (STRA, 10000) I
                  CALL PCKSTR (1, STRA)
                  CH = ' '
                  IF (I .GT. NELBL1) CH = '*'
                  WRITE (*, 10010) CH, IDELB(I), STRA(:LSTRA),
     &               NUMELB(I), NUMLNK(I),
     &               STR2, 'combined into ID ', STRB(:LSTRB)
  100          CONTINUE
            ELSE
               CONTINUE
            END IF
         ELSE IF (ISTAT(IELB) .LT. 0) THEN
            WRITE (STRA, 10000) IELB
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (IELB .GT. NELBL1) CH = '*'
            WRITE (*, 10010) CH, IDELB(IELB), STRA(:LSTRA),
     &         NUMELB(IELB), NUMLNK(IELB), '<deleted>'
         END IF
  110 CONTINUE

      RETURN
10010  FORMAT (2X, A1, 1X, 'Block', I10, 1X, A, ':',
     &   I9, ' elements', I3, '-node', :, 1X, A, :, 1X, 10A)
      END
