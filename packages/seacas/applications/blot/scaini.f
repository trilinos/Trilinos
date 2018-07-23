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

C $Log: scaini.f,v $
C Revision 1.2  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:10:45  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:56:56  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SCAINI (NELB, NVAR, ISTMN)
C=======================================================================

C   --*** SCAINI *** (BLOT) Initialize variable min/max
C   --   Written by Amy Gilkey - revised 04/01/88
C   --
C   --SCAINI initializes the variable minimums and maximums to uncalculated.
C   --ISTMN is the "state" of the minimums and maximums: <0 if not calculated,
C   --0 if no elements in element block, +n if calculated
C   --
C   --Parameters:
C   --   NELB - IN - the number of element blocks (NELBLK for element,
C   --      else 0)
C   --   NVAR - IN - the number of variables of this type
C   --   ISTMN - OUT - initialized to -999

      INTEGER ISTMN(0:NELB,*)

      DO 110 IXVAR = 1, NVAR
         DO 100 IELB = 0, NELB
            ISTMN(IELB,IXVAR) = -999
  100    CONTINUE
  110 CONTINUE

      RETURN
      END
