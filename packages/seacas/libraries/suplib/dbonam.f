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
      SUBROUTINE DBONAM (NDB,
     &   NDIM, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL,
     &   NAMECO, NAMELB, NAMEHV, NAMEGV, NAMENV, NAMEEV, ISEVOK)
C=======================================================================
C$Id: dbonam.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbonam.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:13:29  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:13:28  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:16  gdsjaar
c Initial revision
c 

C   --*** DBONAM *** (EXOLIB) Write database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBONAM writes the names of the coordinates, the element block types,
C   --and the database variables to the database.  The element block variable
C   --truth table is also written.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NDIM - IN - the number of coordinates per node; written only if >= 0
C   --   NELBLK - IN - the number of element blocks; written only if >= 0
C   --   NVARHI - IN - the number of history variables; written only if >= 0
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables
C   --   NAMECO - IN - the names of the coordinates
C   --   NAMELB - IN - the names of the element block types
C   --   NAMEHV - IN - the names of the history variables
C   --   NAMEGV - IN - the names of the global variables
C   --   NAMENV - IN - the names of the nodal variables
C   --   NAMEEV - IN - the names of the element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --
C   --Database must be positioned in front of coordinate names upon entry;
C   --upon exit positioned after element block variable truth table.

      INTEGER NDB
      INTEGER NDIM, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL
      CHARACTER*8 NAMECO(*)
      CHARACTER*8 NAMELB(*)
      CHARACTER*8 NAMEHV(*)
      CHARACTER*8 NAMEGV(*)
      CHARACTER*8 NAMENV(*)
      CHARACTER*8 NAMEEV(*)
      LOGICAL ISEVOK(*)

C   --Write coordinate names

      IF (NDIM .LT. 0) GOTO 100

      IF (NDIM .GT. 0) THEN
         WRITE (NDB) (NAMECO(I), I=1,NDIM)
      ELSE
         WRITE (NDB) 0
      END IF

C   --Write element block type names

      IF (NELBLK .LT. 0) GOTO 100

      IF (NELBLK .GT. 0) THEN
         WRITE (NDB) (NAMELB(I), I=1,NELBLK)
      ELSE
         WRITE (NDB) 0
      END IF

C   --Write the variable names

      IF (NVARHI .LT. 0) GOTO 100

      WRITE (NDB) NVARHI, NVARGL, NVARNP, NVAREL

      IF (NVARHI + NVARGL + NVARNP + NVAREL .GT. 0) THEN
         WRITE (NDB)
     &      (NAMEHV(I), I=1,NVARHI),
     &      (NAMEGV(I), I=1,NVARGL),
     &      (NAMENV(I), I=1,NVARNP),
     &      (NAMEEV(I), I=1,NVAREL)
      ELSE
         WRITE (NDB) 0
      END IF

C   --Write the element block variable truth table

      CALL DBONM1 (NDB, NELBLK, NVAREL, ISEVOK, ISEVOK, MAX(NELBLK,1))

  100 CONTINUE
      RETURN
      END
