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

C $Log: inpcr.f,v $
C Revision 1.2  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:03:48  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:52:39  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE INPCR (NUMNPF, HIDENP, XN, YN, ZN,
     &   COLMIN, ROWMIN, CRDELT, NUMCOL, NUMROW,
     &   IXNCRS, IXNCRE, NPCR, LNPCR, NPX, NPY)
C=======================================================================

C   --*** INPCR *** (MESH) Initialize the NPCR array
C   --   Written by Amy Gilkey - revised 04/04/88
C   --
C   --INPCR initializes the IXNCRS, IXNCRE, and NPCR arrays to point to the
C   --set of all visible nodes in a each block.
C   --
C   --Parameters:
C   --   NUMNPF - IN - the number of nodes
C   --   HIDENP - IN/OUT - node status (as in HIDDEN)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   COLMIN, ROWMIN - IN - the minimum column and row value and 1 / interval
C   --   CRDELT - the column and row interval reciprical (1 / interval)
C   --   NUMCOL, NUMROW - IN - the number of columns and rows
C   --   IXNCRS, IXNCRE - OUT - the starting and ending NPCR index
C   --      for each column and row
C   --   NPCR - OUT - the nodes indexed by IXNCRS, IXNCRE
C   --   LNPCR - OUT - the length of the NPCR array
C   --   NPX, NPY - SCRATCH - size = NUMNPF

      PARAMETER (KNVIS=0, KNFOVR=10, KNHID=100)

      INTEGER HIDENP(NUMNPF)
      REAL XN(NUMNPF), YN(NUMNPF), ZN(NUMNPF)
      INTEGER IXNCRS(0:NUMCOL,0:NUMROW), IXNCRE(0:NUMCOL,0:NUMROW)
      INTEGER NPCR(NUMNPF)
      INTEGER NPX(NUMNPF), NPY(NUMNPF)

      CALL INIINT ((1+NUMCOL) * (1+NUMROW), 0, IXNCRS)

      DO 100 INP = 1, NUMNPF
         IF (HIDENP(INP) .LE. KNVIS) THEN
            ICOL = INT ((XN(INP) - COLMIN) * CRDELT)
            IROW = INT ((YN(INP) - ROWMIN) * CRDELT)
            NPX(INP) = ICOL
            NPY(INP) = IROW
            IXNCRS(ICOL,IROW) = IXNCRS(ICOL,IROW) + 1
         END IF
  100 CONTINUE

      NHEAP = 1
      DO 120 IROW = 0, NUMROW
         DO 110 ICOL = 0, NUMCOL
            L = IXNCRS(ICOL,IROW)
            IXNCRS(ICOL,IROW) = NHEAP
            IXNCRE(ICOL,IROW) = NHEAP - 1
            NHEAP = NHEAP + L
  110    CONTINUE
  120 CONTINUE

      LNPCR = NHEAP - 1

      DO 130 INP = 1, NUMNPF
         IF (HIDENP(INP) .LE. KNVIS) THEN
            ICOL = NPX(INP)
            IROW = NPY(INP)
            NHEAP = IXNCRE(ICOL,IROW) + 1
            IXNCRE(ICOL,IROW) = NHEAP
            NPCR(NHEAP) = INP
         END IF
  130 CONTINUE

      RETURN
      END
