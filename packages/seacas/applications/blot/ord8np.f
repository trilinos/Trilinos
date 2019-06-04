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

C $Log: ord8np.f,v $
C Revision 1.2  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:06:27  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:54:25  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE ORD8NP (NELBLK, LENE, NLNKE, LINKE)
C=======================================================================

C   --*** ORD8NP *** (BLOT) Order 8-node 2D elements
C   --   Modified by John H. Glick -- 10/25/88
C   --   Written by Amy Gilkey, revised 10/29/87
C   --
C   --ORD8NP orders the nodes of 8 and 9 node 2D elements so that they are
C   --connected.  Normally nodes are ordered 1-2-3-4-5-6-7-8-9, where nodes
C   --1 to 4 are corner nodes, 5 to 8 are side nodes, and 9 is the center
C   --node.  The nodes are returned as 1-5-2-6-3-7-4-8-9.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LINK - IN/OUT - the element connectivity

      INTEGER LENE(0:*), LINKE(*)
      INTEGER NLNKE(*)

      INTEGER L(9)

      DO 140 IELB = 1,  NELBLK
         IF ((NLNKE(IELB) .EQ. 8)) THEN
            IXL0 = IDBLNK (IELB, 0, LENE, NLNKE) - 1
            DO 110 IEL = LENE(IELB-1)+1, LENE(IELB)
               L(1) = LINKE(IXL0+1)
               L(2) = LINKE(IXL0+5)
               L(3) = LINKE(IXL0+2)
               L(4) = LINKE(IXL0+6)
               L(5) = LINKE(IXL0+3)
               L(6) = LINKE(IXL0+7)
               L(7) = LINKE(IXL0+4)
               L(8) = LINKE(IXL0+8)
               DO 100 K = 1, NLNKE(IELB)
                  LINKE(IXL0+K) = L(K)
  100          CONTINUE
               IXL0 = IXL0 + NLNKE(IELB)
  110       CONTINUE
         ELSE IF ((NLNKE(IELB) .EQ. 9)) THEN
            IXL0 = IDBLNK (IELB, 0, LENE, NLNKE) - 1
            DO 130 IEL = LENE(IELB-1)+1, LENE(IELB)
               L(1) = LINKE(IXL0+1)
               L(2) = LINKE(IXL0+5)
               L(3) = LINKE(IXL0+2)
               L(4) = LINKE(IXL0+6)
               L(5) = LINKE(IXL0+3)
               L(6) = LINKE(IXL0+7)
               L(7) = LINKE(IXL0+4)
               L(8) = LINKE(IXL0+8)
               L(9) = LINKE(IXL0+9)
               DO 120 K = 1, NLNKE(IELB)
                  LINKE(IXL0+K) = L(K)
  120          CONTINUE
               IXL0 = IXL0 + NLNKE(IELB)
  130       CONTINUE
         END IF
  140 CONTINUE

      RETURN
      END
