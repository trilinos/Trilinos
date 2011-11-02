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
      SUBROUTINE CKNPS (NUMNPS, LNPSNL, NUMNP,
     &   IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS, ICHECK)
C=======================================================================

C   --*** CKNPS *** (GROPE) Check database nodal point sets
C   --
C   --CKNPS checks the nodal point set information.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NUMNPS - IN - the number of nodal point sets
C   --   LNPSNL - IN - the number of nodes for all sets
C   --   NUMNP - IN - the number of nodes
C   --   IDNPS - IN - the nodal point set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   IXNNPS - IN - the index of the first node for each set
C   --   LTNNPS - IN - the nodes for all sets
C   --   FACNPS - IN - the distribution factors for all sets
C   --   ICHECK - SCRATCH - size = MAX (NUMNPS, LNPSNL)

      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)
      INTEGER ICHECK(*)

      CHARACTER*5 STRA

C   --Check for unique identifier

      DO 100 INPS = 1, NUMNPS
         IF (LOCINT (IDNPS(INPS), INPS-1, IDNPS) .GT. 0) THEN
            CALL INTSTR (1, 0, IDNPS(INPS), STRA, LSTRA)
            CALL PRTERR ('CMDSPEC', 'Nodal point set ID '
     &         // STRA(:LSTRA) // ' is not unique')
         END IF
  100 CONTINUE

C   --Check number of nodes in nodal point sets

      NNPS = 0
      DO 110 INPS = 1, NUMNPS
         NNPS = MAX (NNPS, IXNNPS(INPS) + NNNPS(INPS) - 1)
  110 CONTINUE

      IF (NNPS .NE. LNPSNL) THEN
         CALL PRTERR ('CMDSPEC', 'Maximum node index' // 
     &      ' in all nodal point sets does not match total')
      END IF

C   --Check all nodes in node point sets are within node range

      CALL CHKRNG (LTNNPS, LNPSNL, NUMNP, NZERO, NERR)
      IF (NERR .GT. 0) THEN
        CALL PRTERR ('CMDSPEC',
     &    'Nodal point set node ids are out of range')
      END IF
      IF (NZERO .GT. 0) THEN
        CALL PRTERR ('CMDSPEC',
     &    'Nodal point set node ids are zero')
      END IF
      
      RETURN
      END
