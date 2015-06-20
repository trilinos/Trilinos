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
      SUBROUTINE DBOELB (NDB, NELBS, NELBE,
     &   IDELB, NUMELB, NUMLNK, NUMATR, LINK, ATRIB)
C=======================================================================
C$Id: dboelb.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dboelb.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:13:19  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:13:17  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:14  gdsjaar
c Initial revision
c 

C   --*** DBOELB *** (EXOLIB) Write database element blocks
C   --   Written by Amy Gilkey - revised 10/12/87
C   --
C   --DBOELB writes the element block information to the database.
C   --Some dynamic dimensioning is done.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NELBS, NELBE - IN - the number of first and last element blocks
C   --      to write
C   --   IDELB - IN - the element block IDs for each block
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   NUMATR - IN - the number of attributes in each block
C   --   LINK - IN - the connectivity for each block
C   --   ATRIB - IN - the attributes for each block
C   --
C   --Database must be positioned at start of element block information
C   --upon entry; upon exit at end of element block information.

      INTEGER NDB
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      INTEGER LINK(*)
      REAL ATRIB(*)

      IELNK = 0
      IEATR = 0

      DO 100 NELB = NELBS, NELBE
         IELB = NELB-NELBS+1

         WRITE (NDB)
     &      IDELB(IELB), NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB)

         ISLNK = IELNK + 1
         IELNK = IELNK + NUMLNK(IELB) * NUMELB(IELB)
         ISATR = IEATR + 1
         IEATR = IEATR + NUMATR(IELB) * NUMELB(IELB)

         CALL DBOEB1 (NDB, NELB,
     &      NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB),
     &      LINK(ISLNK), ATRIB(ISATR),
     $      MAX(NUMLNK(IELB),1), MAX(NUMATR(IELB),1))
  100 CONTINUE

      RETURN
      END
