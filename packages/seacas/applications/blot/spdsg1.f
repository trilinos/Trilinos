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

C $Log: spdsg1.f,v $
C Revision 1.2  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:14:10  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:58:15  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SPDSG1 (NNENUM, NENUM, IE2ELB, ISEVOK, NEL, ISEGEL)
C=======================================================================

C   --*** SPDSG1 *** (SPLOT) Find all defined element values
C   --   Written by Amy Gilkey - revised 11/05/87
C   --
C   --SPDSG1 returns a list of the defined element indices for an
C   --element variable.
C   --
C   --Parameters:
C   --   NNENUM - IN - the number of selected nodes/elements
C   --   NENUM - IN - the selected nodes/elements
C   --   IE2ELB - IN - the element block for each element
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable of block j exists iff ISEVOK(j)
C   --   NEL - OUT - the number of defined elements
C   --   ISEGEL - OUT - the NENUM indices of the defined elements;
C   --      ISEGEL(0) = the number of defined elements

      INTEGER NENUM(*)
      INTEGER IE2ELB(*)
      LOGICAL ISEVOK(*)
      INTEGER ISEGEL(0:*)

      NEL = 0
      DO 100 I = 1, NNENUM
         IF (ISEVOK(IE2ELB(NENUM(I)))) THEN
            NEL = NEL + 1
            ISEGEL(NEL) = I
         END IF
  100 CONTINUE
      ISEGEL(0) = NEL

      RETURN
      END
