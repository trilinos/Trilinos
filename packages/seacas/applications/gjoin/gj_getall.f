C Copyright (c) 2008-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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
C 
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
C 

C=======================================================================
      SUBROUTINE GETALL (MATCH, LLIST, LIST, NINSET, INSET)
C=======================================================================
C $Id: getall.f,v 1.1 1999/01/18 19:21:21 gdsjaar Exp $
C $Log: getall.f,v $
C Revision 1.1  1999/01/18 19:21:21  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:25  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:34:44  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:34:43  gdsjaar
c Initial revision
c 

C   --*** GETALL *** (GJOIN) Get all items that match
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --GETALL checks a list of values and retrieves all items which
C   --match the given values.
C   --
C   --Parameters:
C   --   MATCH - IN - the value to match
C   --   LLIST - IN - the length of LIST
C   --   LIST - IN - the list of values
C   --   NINSET - OUT - the number of matching values
C   --   INSET - OUT - the indices of the matching items

      INTEGER LIST(*)
      INTEGER INSET(*)

      NINSET = 0
      DO 100 I = 1, LLIST
         IF (LIST(I) .EQ. MATCH) THEN
            NINSET = NINSET + 1
            INSET(NINSET) = I
         END IF
  100 CONTINUE

      RETURN
      END
