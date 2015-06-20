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
      SUBROUTINE DBOESS (NDB, NUMESS, LESSEL, LESSNL,
     &   IDESS, NEESS, NNESS, IXEESS, IXNESS, LTEESS, LTNESS, FACESS)
C=======================================================================
C$Id: dboess.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dboess.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:13:21  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:13:20  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:15  gdsjaar
c Initial revision
c 

C   --*** DBOESS *** (EXOLIB) Write database node sets
C   --   Written by Amy Gilkey - revised 02/27/86
C   --
C   --DBOESS writes the side set information to the database.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NUMESS - IN - the number of side sets
C   --   LESSEL - IN - the length of the side sets element list
C   --   LESSNL - IN - the length of the side sets node list
C   --   IDESS - IN - the side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   NNESS - IN - the number of nodes for each set
C   --   IXEESS - IN - the index of the first element for each set
C   --   IXNESS - IN - the index of the first node for each set
C   --   LTEESS - IN - the elements for all sets
C   --   LTNESS - IN - the nodes for all sets
C   --   FACESS - IN - the distribution factors for all sets
C   --
C   --Database must be positioned at start of side set information
C   --upon entry; upon exit at end of side set information.

      INTEGER NDB
      INTEGER NUMESS, LESSEL, LESSNL
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER IXEESS(*)
      INTEGER IXNESS(*)
      INTEGER LTEESS(*)
      INTEGER LTNESS(*)
      REAL FACESS(*)

      IF (NUMESS .GT. 0) THEN
         WRITE (NDB) (IDESS(IESS), IESS=1,NUMESS)
         WRITE (NDB) (NEESS(IESS), IESS=1,NUMESS)
         WRITE (NDB) (NNESS(IESS), IESS=1,NUMESS)
         WRITE (NDB) (IXEESS(IESS), IESS=1,NUMESS)
         WRITE (NDB) (IXNESS(IESS), IESS=1,NUMESS)
      ELSE
         WRITE (NDB) 0
         WRITE (NDB) 0
         WRITE (NDB) 0
         WRITE (NDB) 0
         WRITE (NDB) 0
      END IF
      IF (LESSEL .GT. 0) THEN
         WRITE (NDB) (LTEESS(NL), NL=1,LESSEL)
      ELSE
         WRITE (NDB) 0
      END IF
      IF (LESSNL .GT. 0) THEN
         WRITE (NDB) (LTNESS(NL), NL=1,LESSNL)
         WRITE (NDB) (FACESS(NL), NL=1,LESSNL)
      ELSE
         WRITE (NDB) 0
         WRITE (NDB) 0
      END IF

      RETURN
      END
