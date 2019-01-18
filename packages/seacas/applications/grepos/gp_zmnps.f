C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C * Neither the name of NTESS nor the names of its
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

C=======================================================================
      SUBROUTINE ZMNPS (NUMNPS, ISTAT, LNPSNL, LNPSDF,
     &   IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS, NAME)
C=======================================================================
C   --*** ZMNPS *** (GJOIN) Compress nodal point sets
C   --   Written by Amy Gilkey - revised 01/20/88
C   --
C   --ZMNPS compresses the nodal point sets by removing deleted nodes.
C   --Assumes that the nodes are already renumbered and ordered.
C   --
C   --Parameters:
C   --   NUMNPS - IN - the number of nodal point sets
C   --   LNPSNL - IN/OUT - the length of the nodal point sets node list
C   --   LNPSDF - IN/OUT - the length of the dist.fact. node list
C   --   IDNPS - IN - the nodal point set ID for each set
C   --   NNNPS - IN/OUT - the number of nodes for each set
C   --   IXNNPS - IN/OUT - the index of the first node for each set
C   --   LTNNPS - IN/OUT - the nodes for all sets
C   --   FACNPS - IN/OUT - the distribution factors for all sets

      include 'gp_namlen.blk'

      INTEGER ISTAT(*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)
      CHARACTER*(maxnam) NAME(*)

      IF (NUMNPS .LE. 0) RETURN

      IF (LNPSDF .NE. LNPSNL .AND. LNPSDF .NE. 0) then
        call prterr('ERROR',
     *    'Length of nodeset dist. fact. list must be zero or '//
     *    'equal to length of nodeset node list. It is neither.')
      end if

      JNPS = 0
      JNN = 0
      DO 110 INPS = 1, NUMNPS
         JNNLST = JNN
         DO 100 NN = IXNNPS(INPS), IXNNPS(INPS)+NNNPS(INPS)-1
            IF (LTNNPS(NN) .GT. 0) THEN
               JNN = JNN + 1
               LTNNPS(JNN) = LTNNPS(NN)
               if (LNPSDF .NE. 0) then
                 FACNPS(JNN) = FACNPS(NN)
               end if
            END IF
  100    CONTINUE
         N = JNN - JNNLST
         IF (N .GT. 0) THEN
            JNPS = JNPS + 1
            IDNPS(JNPS) = IDNPS(INPS)
            NNNPS(JNPS) = N
            IXNNPS(JNPS) = JNNLST + 1
            NAME(JNPS) = NAME(INPS)
          ELSE
            ISTAT(INPS) = -IDNPS(INPS)
         END IF
  110 CONTINUE

      NUMNPS = JNPS
      LNPSNL = JNN
      if (lnpsdf .ne. 0) lnpsdf = jnn
        
      RETURN
      END
