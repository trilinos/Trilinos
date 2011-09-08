C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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

C $Id: rdnps.f,v 1.4 2007/10/17 18:47:22 gdsjaar Exp $
C=======================================================================
      SUBROUTINE RDNPS (NTXT, NUMNPS, LNPSNL, LNPSDF,
     &   IDNPS, NNNPS, NDNPS, IXNNPS, IXDNPS, LSTNPS, FACNPS, *)
C=======================================================================

C   --*** RDNPS *** (TXTEXO) Read database node sets
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --RDNPS reads the node set information from the text file.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NUMNPS - IN - the number of node sets
C   --   LNPSNL - IN - the length of the node sets node list
C   --   IDNPS - OUT - the node set ID for each set
C   --   NNNPS - OUT - the number of nodes for each set
C   --   IXNNPS - OUT - the index of the first node for each set
C   --   LSTNPS - OUT - the nodes for all sets
C   --   FACNPS - OUT - the distribution factors for all sets
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of node set information
C   --upon entry; upon exit at end of node set information.

      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER NDNPS(*)
      INTEGER IXNNPS(*)
      INTEGER IXDNPS(*)
      INTEGER LSTNPS(*)
      REAL FACNPS(*)

      CHARACTER*5 STRA, STRB

      NN = 0
      ND = 0
      READ (NTXT, *, END=120, ERR=120)
      READ (NTXT, *, END=120, ERR=120)
      READ (NTXT, *, END=120, ERR=120)
      DO 110 INPS = 1, NUMNPS
         READ (NTXT, *, END=120, ERR=120)
         READ (NTXT, *, END=120, ERR=120) IDNPS(INPS), NNNPS(INPS),
     &        NDNPS(INPS)

         IXNNPS(INPS) = NN + 1
         NN = NN + NNNPS(INPS)
         IXDNPS(INPS) = ND + 1
         ND = ND + NDNPS(INPS)
         
         DO 100 NL = IXNNPS(INPS), NN
            READ (NTXT, *, END=130, ERR=130) LSTNPS(NL), FACNPS(NL)
  100    CONTINUE
  110 CONTINUE

      IF (NN .NE. LNPSNL) THEN
         CALL INTSTR (1, 0, NN, STRA, LSTRA)
         CALL INTSTR (1, 0, LNPSNL, STRB, LSTRB)
         CALL PRTERR ('WARNING',
     &      'NODE SET NUMBER OF NODES = ' // STRA(:LSTRA)
     &      // ' does not match TOTAL = ' // STRB(:LSTRB))
      END IF

      RETURN

  120 CONTINUE
      CALL INTSTR (1, 0, INPS, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &   'Reading HEADER DATA for NODE SET ' // STRA(:LSTRA))
      GOTO 140
  130 CONTINUE
      CALL INTSTR (1, 0, NL-IXNNPS(INPS)+1, STRA, LSTRA)
      CALL INTSTR (1, 0, INPS, STRB, LSTRB)
      CALL PRTERR ('FATAL',
     &   'Reading NODES and FACTORS for node ' // STRA(:LSTRA)
     &   // ' for NODE SET ' // STRB(:LSTRB))
      GOTO 140
  140 CONTINUE
      RETURN 1
      END
