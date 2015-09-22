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
      SUBROUTINE DBONPS (NDB, NUMNPS, LNPSNL,
     &   IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS)
C=======================================================================
C$Id: dbonps.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbonps.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:13:34  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:13:33  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:17  gdsjaar
c Initial revision
c 

C   --*** DBONPS *** (EXOLIB) Write database node sets
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --DBONPS writes the node set information to the database.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NUMNPS - IN - the number of node sets
C   --   LNPSNL - IN - the length of the node sets node list
C   --   IDNPS - IN - the node set ID for each set
C   --   NNNPS - IN - the number of nodes for each set
C   --   IXNNPS - IN - the index of the first node for each set
C   --   LTNNPS - IN - the nodes for all sets
C   --   FACNPS - IN - the distribution factors for all sets
C   --
C   --Database must be positioned at start of node set information
C   --upon entry; upon exit at end of node set information.

      INTEGER NDB
      INTEGER NUMNPS, LNPSNL
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)

      IF (NUMNPS .GT. 0) THEN
         WRITE (NDB) (IDNPS(INPS), INPS=1,NUMNPS)
         WRITE (NDB) (NNNPS(INPS), INPS=1,NUMNPS)
         WRITE (NDB) (IXNNPS(INPS), INPS=1,NUMNPS)
      ELSE
         WRITE (NDB) 0
         WRITE (NDB) 0
         WRITE (NDB) 0
      END IF
      IF (LNPSNL .GT. 0) THEN
         WRITE (NDB) (LTNNPS(NL), NL=1,LNPSNL)
         WRITE (NDB) (FACNPS(NL), NL=1,LNPSNL)
      ELSE
         WRITE (NDB) 0
         WRITE (NDB) 0
      END IF

      RETURN
      END
