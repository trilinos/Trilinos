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
      SUBROUTINE ORDSTR (NORD, IXORD, LOLD, IOLD, ISCR, INEW)
C=======================================================================
C $Id: ordstr.f,v 1.1 1999/01/18 19:21:24 gdsjaar Exp $
C $Log: ordstr.f,v $
C Revision 1.1  1999/01/18 19:21:24  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:27  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:35:24  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:35:23  gdsjaar
c Initial revision
c 

C   --*** ORDSTR *** (GJOIN) Order a list of strings according to indices
C   --   Written by Greg Sjaardema - revised 07/11/90
C   --   Modified from ORDIX Written by Amy Gilkey 
C   --
C   --ORDSTR orders a list of strings according to a list of indices.
C   --
C   --Parameters:
C   --   NORD - IN - the number of indices
C   --   IXORD - IN - the indices of the ordered items
C   --   LOLD - IN - the length of IOLD
C   --   IOLD - IN - the unordered string list
C   --   ISCR - SCRATCH - size = LOLD
C   --   INEW - OUT - the ordered string list

      include 'exodusII.inc'

      INTEGER IXORD(*)
      character*(MXSTLN) iold(*)
      character*(MXSTLN) iscr(*)
      character*(MXSTLN) inew(*)

      DO 100 I = 1, LOLD
         ISCR(I) = IOLD(I)
  100 CONTINUE
      DO 110 I = 1, NORD
         INEW(I) = ISCR(IXORD(I))
  110 CONTINUE

      RETURN
      END
