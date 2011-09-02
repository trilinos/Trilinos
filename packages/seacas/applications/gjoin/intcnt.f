C Copyright (c) 2008 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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

C=======================================================================
      INTEGER FUNCTION INTCNT (INT, LIST, LLIST)
C=======================================================================
C $Id: intcnt.f,v 1.1 1999/01/18 19:21:22 gdsjaar Exp $
C $Log: intcnt.f,v $
C Revision 1.1  1999/01/18 19:21:22  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:26  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:34:54  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:34:53  gdsjaar
c Initial revision
c 

C   --*** INTCNT *** (GJOIN) Returns the number of occurances in list
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --INTCNT returns the number of times an integer is found in a list.
C   --
C   --Parameters:
C   --   INT - IN - the integer to count
C   --   LIST - IN - the list of integers
C   --   LLIST - IN - the length of LIST

      INTEGER LIST(*)

      INTCNT = 0
      DO 100 I = 1, LLIST
         IF (LIST(I) .EQ. INT) THEN
            INTCNT = INTCNT + 1
         END IF
  100 CONTINUE

      RETURN
      END
