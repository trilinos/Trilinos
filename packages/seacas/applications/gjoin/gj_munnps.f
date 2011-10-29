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
      SUBROUTINE MUNNPS (NUMNPS, ISTAT, LNPSNL,
     &   IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS,
     &   LTNX, FACX, IXNPS, NNX, ISCR, NODSCR, NUMNP)
C=======================================================================
C $Id: munnps.f,v 1.2 2008/07/31 20:15:56 gdsjaar Exp $
C $Log: munnps.f,v $
C Revision 1.2  2008/07/31 20:15:56  gdsjaar
C Change the way the nodal point node membership is calculated. For some
C reason, the locint calls used in the old method became a bottleneck on
C some compilers (runtime went from ~30 seconds to multiple hours) with
C no other changes.  Recompiling the same code with 32-bit gcc gives
C short runtime; 32-bit or 64-bit intel gives long runtimes.
C
C Modified the routine to index all active nodes in a NUMNP-long array
C and then check whether a node is indexed instead of doing a locint
C call (constant time vs linear time).
C
C Revision 1.1  1999/01/18 19:21:23  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:26  a294617
c Initial import == gjoin 1.36
c
C Revision 1.2  1997/11/18 15:55:29  gdsjaar
C Changes to improve efficiency (factor of 2 to 15 on some Goodyear Treads)
C
C Redid nodeset munching to reduce number of locint calls. First do a
C quick scan on array to get maximum node id in searched array. Then,
C before doing 'locint' call, ensure that id being searched for is less
C than maximum (most times it won't be, saves searching entire list).
C
C Modified the node matching code to also index the nodes within the
C overlap region if doing a nodeset match. This has greatest benefit if
C not all nodes in the nodeset will be matched.
C
C Minor change in offset f -- 'dimension' to 'real'
C
C Revision 1.1.1.1  1990/11/12 14:35:15  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:35:14  gdsjaar
c Initial revision
c 

C   --*** MUNNPS *** (GJOIN) Compress and rearrange nodal point sets
C   --   Written by Amy Gilkey - revised 02/25/88
C   --
C   --MUNNPS processes the nodal point sets according to the set status.
C   --Sets may be combined or deleted.
C   --
C   --Parameters:
C   --   NUMNPS - IN/OUT - the number of nodal point sets
C   --   ISTAT - IN - the status of each set:
C   --      0 = same
C   --      - = delete
C   --      n = combine with set n
C   --   LNPSNL - IN/OUT - the length of the nodal point sets node list
C   --   IDNPS - IN/OUT - the nodal point set ID for each set
C   --   NNNPS - IN/OUT - the number of nodes for each set
C   --   IXNNPS - IN/OUT - the index of the first node for each set
C   --   LTNNPS - IN/OUT - the nodes for all sets
C   --   FACNPS - IN/OUT - the distribution factors for all sets
C   --   LTNX - SCRATCH - sized to hold the set nodes
C   --   FACX - SCRATCH - sized to hold the set factors
C   --   IXNPS - SCRATCH - size = NUMNPS
C   --   NNX - SCRATCH - size = NUMNPS
C   --   ISCR - SCRATCH - size = NUMNPS
C   --   NUMNP - IN -- number of nodes in model
C   --   NODSCR - SCRATCH - size = NUMNOD

      INTEGER ISTAT(*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*), LTNX(*)
      REAL FACNPS(*), FACX(*)
      INTEGER IXNPS(*)
      INTEGER NNX(*)
      INTEGER ISCR(*)
      INTEGER NODSCR(*)
      
      IF (NUMNPS .LE. 0) RETURN

      JNPS = 0
      JNN0 = 0
      DO 120 INPS = 1, NUMNPS

         IF (ISTAT(INPS) .EQ. 0) THEN
            NINSET = 1
            ISCR(NINSET) = INPS
         ELSE IF (ISTAT(INPS) .EQ. INPS) THEN
            CALL GETALL (INPS, NUMNPS, ISTAT, NINSET, ISCR)
         ELSE
            NINSET = 0
         END IF

         IF (NINSET .GT. 0) THEN
            JNPS = JNPS + 1
            IXNPS(JNPS) = INPS
            NNX(JNPS) = 0
            JOLD = JNN0+1
         END IF

         DO 110 ISET = 1, NINSET
            N = ISCR(ISET)
            JLEN = JNN0+1 - JOLD
            INN0 = IXNNPS(N) - 1

            call iniint(numnp, 0, nodscr)
            do 90 i=jold, jold+jlen-1
              nodscr(ltnx(i)) = 1
 90         continue

            NNEW = 0
            DO 100 I = 1, NNNPS(N)
               IF (nodscr(ltnnps(inn0+i)) .eq. 0) then
                  NNEW = NNEW + 1
                  LTNX(JNN0+NNEW) = LTNNPS(INN0+I)
                  FACX(JNN0+NNEW) = FACNPS(INN0+I)
               END IF
  100       CONTINUE
            NNX(JNPS) = NNX(JNPS) + NNEW
            JNN0 = JNN0 + NNEW
  110    CONTINUE
  120 CONTINUE

      CALL ORDIX (JNPS, IXNPS, NUMNPS, IDNPS, ISCR, IDNPS)
      CALL MOVINT (JNPS, NNX, NNNPS)
      NUMNPS = JNPS
      JNN = 1
      DO 130 INPS = 1, NUMNPS
         IXNNPS(INPS) = JNN
         JNN = JNN + NNNPS(INPS)
  130 CONTINUE
      LNPSNL = JNN - 1

      CALL MOVINT (LNPSNL, LTNX, LTNNPS)
      CALL MOVREA (LNPSNL, FACX, FACNPS)

      RETURN
      END
