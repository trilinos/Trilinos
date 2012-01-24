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
      SUBROUTINE NEWNPS (IDFRO, IDBCK,
     &   IDNPS, NNNPS, NNNP3, IXNNPS, IXNNP3,
     &   LTNNPS, LTNNP3, FACNPS, FACNP3,
     &   IXNP, NRNP)
C=======================================================================

C   --*** NEWNPS *** (GEN3D) Calculate 3D node sets
C   --   Written by Amy Gilkey - revised 05/05/86
C   --
C   --NEWNPS calculates the node set information for the 3D database.
C   --The front and back node sets are not calculated (since they
C   --are easily derived), and they are not included in the length tally.
C   --
C   --Parameters:
C   --   IDFRO - IN - ids for front surface node sets; (0) = length
C   --   IDBCK - IN - ids for back surface node sets; (0) = length
C   --   IDNPS - IN - the 2D node sets ids
C   --   NNNPS - IN - the number of nodes for each 2D set
C   --   NNNP3 - OUT - the number of nodes for each 3D set
C   --   IXNNPS - IN - the index of the first node for each 2D set
C   --   IXNNP3 - OUT - the index of the first node for each 3D set
C   --   LTNNPS - IN - the nodes for all 2D sets
C   --   LTNNP3 - OUT - the nodes for all 3D sets
C   --   FACNPS - IN - the distribution factors for all 2D sets
C   --   FACNP3 - OUT - the distribution factors for all 3D sets
C   --   IXNP - IN - the new index for each node
C   --   NRNP - IN - the number of new nodes generated for each node
C   --
C   --Common Variables:
C   --   Uses NDBOUT of /DBASE/
C   --   Uses NUMNPS, LNPSNL of /DBNUMS/
C   --   Sets LNPSNO of /DBNUM3/
C   --   Uses NNREPL of /PARAMS/

      include 'exodusII.inc'
      INCLUDE 'gs_dbtitl.blk'
      INCLUDE 'gs_dbnums.blk'
      INCLUDE 'gs_dbnum3.blk'
      INCLUDE 'gs_params.blk'

      INTEGER IDFRO(0:*)
      INTEGER IDBCK(0:*)
      INTEGER IDNPS(*)
      INTEGER NNNPS(*), NNNP3(*)
      INTEGER IXNNPS(*), IXNNP3(*)
      INTEGER LTNNPS(*), LTNNP3(*)
      REAL FACNPS(*), FACNP3(*)
      INTEGER IXNP(*), NRNP(*)

C   --Node set ids - unchanged

      CONTINUE

C   --Number of nodes in all sets

      RETURN
      END
