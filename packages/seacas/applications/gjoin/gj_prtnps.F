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
      SUBROUTINE PRTNPS (ISTAT, NNPS1, NNPS2, IDNPS, NNNPS, ISCR)
C=======================================================================
C $Id: prtnps.F,v 1.1 1999/01/18 19:21:25 gdsjaar Exp $
C $Log: prtnps.F,v $
C Revision 1.1  1999/01/18 19:21:25  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:29  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1  1998/04/09 05:02:04  gdsjaar
C The code which prints out the sidesets/nodesets/blocks for id
C combining/renumbering needed the needsdoubleescape if block for
C certain systems. Needed to move those routines from .f to .F
C
C Revision 1.3  1997/04/09 22:10:09  gdsjaar
C Modified the node matching code to sort the array on the coordinate
C with the largest span. Previously it was done on the x coordinate
C always.
C
C Upped the format string used for the number of elements in the prtelb,
C prtess, and prtnps output to handle large models.
C
C Revision 1.2  1993/07/23 22:34:31  gdsjaar
C Cleaned up unused variables, removed unneeded routine
C
c Revision 1.1.1.1  1990/11/12  14:35:47  gdsjaar
c GJOIN - X1.00.40 - 7/17/90
c
c Revision 1.1  90/11/12  14:35:46  gdsjaar
c Initial revision
c 

C   --*** PRTNPS *** (GJOIN) Print summary for nodal point sets
C   --   Written by Amy Gilkey - revised 02/04/88
C   --
C   --PRTNPS prints a summary for the nodal point sets, including their
C   --status.
C   --
C   --Parameters:
C   --   ISTAT - IN - the status of each nodal point set:
C   --      0 = same
C   --      - = delete
C   --      n = combine with block n
C   --   NNPS1 - IN - the number of nodal point sets in the first database
C   --   NNPS2 - IN - the number of nodal point sets in the second database
C   --   IDNPS - IN - the nodal point set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   ISCR - SCRATCH - size = NNPS1+NNPS2

      INTEGER ISTAT(*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER ISCR(*)

      CHARACTER*2 STR2
      CHARACTER*20 STRA, STRB
      CHARACTER CH

      NUMNPS = NNPS1 + NNPS2

      IF (NUMNPS .LE. 0) RETURN

      WRITE (*, *)
      WRITE (*, *) 'Nodal point sets:'

      WRITE (STRA, 10000) NUMNPS
10000  FORMAT ('(#', I5, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)

      DO 110 INPS = 1, NUMNPS
         IF (ISTAT(INPS) .EQ. 0) THEN
            WRITE (STRA, 10000) INPS
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (INPS .GT. NNPS1) CH = '*'
            WRITE (*, 10010) CH, IDNPS(INPS), STRA(:LSTRA),
     &         NNNPS(INPS)
         ELSE IF (ISTAT(INPS) .GT. 0) THEN
            IF (ISTAT(INPS) .EQ. INPS) THEN
               CALL GETALL (INPS, NUMNPS, ISTAT, NINSET, ISCR)
               WRITE (STRB, '(I5)') IDNPS(INPS)
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
                  IF (I .GT. NNPS1) CH = '*'
                  WRITE (*, 10010) CH, IDNPS(I), STRA(:LSTRA),
     &               NNNPS(I),
     &               STR2, 'combined into ID ', STRB(:LSTRB)
  100          CONTINUE
            ELSE
               CONTINUE
            END IF
         ELSE IF (ISTAT(INPS) .LT. 0) THEN
            WRITE (STRA, 10000) INPS
            CALL PCKSTR (1, STRA)
            CH = ' '
            IF (INPS .GT. NNPS1) CH = '*'
            WRITE (*, 10010) CH, IDNPS(INPS), STRA(:LSTRA),
     &         NNNPS(INPS), '<deleted>'
         END IF
  110 CONTINUE

      RETURN
10010  FORMAT (2X, A1, 1X, 'Set', I6, 1X, A, ':',
     &   I9, ' nodes', :, 3X, A, :, 3X, 5A)
      END
