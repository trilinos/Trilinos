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
      SUBROUTINE RENELB (NELBLK, IOFFNP, IXNP, NUMELB, NUMLNK, LINK)
C=======================================================================
C $Id: renelb.f,v 1.1 1999/01/18 19:21:26 gdsjaar Exp $
C $Log: renelb.f,v $
C Revision 1.1  1999/01/18 19:21:26  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:27  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:36:00  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:35:58  gdsjaar
c Initial revision
c 

C   --*** RENELB *** (GJOIN) Renumber connectivity in element blocks
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --RENELB renumbers the nodes in the connectivity arrays.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   IOFFNP - IN - the nodal offset: if positive, add to node number;
C   --      if negative, use IXNP
C   --   IXNP - IN - the new node number for each node
C   --   NUMELB - IN - the number of elements in each block
C   --   NUMLNK - IN - the number of nodes per element in each block
C   --   LINK - IN - the connectivity for each block

      INTEGER IXNP(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER LINK(*)

      ILNK = 1

      DO 100 IELB = 1, NELBLK
         CALL RENIX (NUMELB(IELB) * NUMLNK(IELB),
     &      IOFFNP, IXNP, LINK(ILNK))
         ILNK = ILNK + NUMLNK(IELB) * NUMELB(IELB)
  100 CONTINUE

      RETURN
      END
