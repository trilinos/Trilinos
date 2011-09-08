C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
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
C    * Neither the name of Sandia Corporation nor the names of its
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
C=======================================================================
      SUBROUTINE WRSTEP (NDBOUT, ISTEP, MAXNE, VARVAL, VISELB,
     &   IXNODE, IXELB, IXELBO, IXELEM, IDELB, ISEVOK, GVSCR,
     &   VARSCR, MERR)
C=======================================================================

C   --*** WRSTEP *** (ALGEBRA) Write database variables for one time step
C   --   Written by Amy Gilkey - revised 05/18/88
C   --
C   --WRSTEP writes the database history, global, nodal, and element variables
C   --for one time step.
C   --
C   --Parameters:
C   --   NDBOUT - IN - the output database file
C   --   MAXNE  - IN - the VARVAL dimension (max of NUMEL and NUMNP)
C   --   VARVAL - IN - the output database variables
C   --   VISELB - IN - true iff element block i is to be written
C   --   IXNODE - IN - the indices of the output nodes (iff NUMNPO <> NUMNP)
C   --   IXELB  - IN - the cumulative element counts for each element block
C   --   IXELBO - IN - the cumulative element counts for each output block;
C   --   IXELEM - IN - the indices of the output elements (iff NUMELO <> NUMEL)
C   --   ISEVOK - IN - the element block variable truth table;
C   --                 variable i of block j exists iff ISEVOK(j,i)
C   --   MERR   - OUT - error code
C   --
C   --Database must be positioned in front of time step upon entry;
C   --upon exit positioned after time step.

      PARAMETER (ICURTM = 1, ILSTTM = 2, IONETM = 3)
      include 'namlen.blk'
      include 'var.blk'
      include 'dbnums.blk'
      include 'dbout.blk'
      include 'dbxvar.blk'
      
      REAL VARVAL(MAXNE,*)
      LOGICAL VISELB(NELBLK)
      INTEGER IXNODE(*), IXELEM(*)
      INTEGER IXELB(0:NELBLK)
      INTEGER IXELBO(0:NELBLK)
      INTEGER IDELB(*)
      LOGICAL ISEVOK(NELBLK,*)
      REAL GVSCR(*)
      REAL VARSCR(MAXNE)

      LOGICAL NOIX

C     write time step
      call exptim(ndbout, istep, VARVAL(IDVAR(ITIME),
     &            ISTVAR(ICURTM,ITIME)), ierr)


C      --Write global variables

      IF (JGVBEG .LE. JGVEND) THEN
         i = 0
         do 10 j=jgvbeg, jgvend
            i = i + 1
            gvscr(i) = varval(idvar(j),istvar(icurtm,j))
 10      continue
         call expgv(ndbout, istep, nvargo, gvscr, ierr)
      END IF

C      --Write nodal variables

      DO 100 J = JNVBEG, JNVEND
         NSTO = ISTVAR(ICURTM,J)
         IF (NUMNPO .GT. 0) THEN
            IF (NUMNP .EQ. NUMNPO) THEN
               call expnv(ndbout, istep, j-jnvbeg+1,
     &                    numnp, VARVAL(1,NSTO), ierr)
            ELSE
               do 30 i=1, numnpo
                  varscr(i) = varval(ixnode(i),nsto)
 30            continue
               call expnv(ndbout, istep, j-jnvbeg+1,
     &                    numnpo, varscr, ierr)
            END IF
          END IF
 100  CONTINUE
           
C      --Write element variables

      NOIX = (IXELB(NELBLK) .EQ. IXELBO(NELBLK))
      DO 120 IELB = 1, NELBLK
         IF (VISELB(IELB)) THEN
            DO 110 J = JEVBEG, JEVEND
               IEV = IEVVAR(J)
               IF (ISEVOK(IELB,IEV)) THEN
                  IF (IXELBO(IELB) .GT. IXELBO(IELB-1)) THEN
                     NSTO = ISTVAR(ICURTM,J)
                     IF (NOIX) THEN
                        nelb = ixelb(ielb) - ixelb(ielb-1)
                        call expev(ndbout, istep, j-jevbeg+1,
     &                          idelb(ielb), nelb,
     &                          varval(ixelb(ielb-1)+1,nsto), ierr)
                     ELSE
                        i = 0
                        do 40 n=IXELBO(IELB-1)+1,IXELBO(IELB)
                           i = i + 1
                           varscr(i) = varval(ixelem(n),nsto)
 40                     continue
                        call expev(ndbout, istep, j-jevbeg+1,
     &                             idelb(ielb), i, varscr, ierr)
                     END IF
                  END IF
               END IF
  110       CONTINUE
         END IF
  120  CONTINUE

      RETURN
      END
