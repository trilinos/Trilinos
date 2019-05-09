C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C $Log: cntelb.f,v $
C Revision 1.2  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 19:56:43  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:48:43  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CNTELB (IELBST, NELBLK, NUMON, NUMSEL)
C=======================================================================

C   --*** CNTELB *** (MESH) Counts selected element blocks
C   --   Written by Amy Gilkey - revised 01/23/87
C   --
C   --CNTELB counts the number of ON element blocks and the number of selected
C   --element blocks.
C   --
C   --Parameters:
C   --   IELBST - IN - the element block status:
C   --      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   NELBLK - IN - the number of element blocks
C   --   NUMON - OUT - the number of ON element blocks
C   --   NUMSEL - OUT - the number of selected element blocks

      INTEGER IELBST(NELBLK)

C   --Count the number of ON element blocks

      NUMON = 0
      DO 100 IELB = 1, NELBLK
         IF (IELBST(IELB) .GE. 0) NUMON = NUMON + 1
  100 CONTINUE

C   --Count the number of selected element blocks

      NUMSEL = 0
      DO 110 IELB = 1, NELBLK
         IF (IELBST(IELB) .GT. 0) NUMSEL = NUMSEL + 1
  110 CONTINUE

      RETURN
      END
