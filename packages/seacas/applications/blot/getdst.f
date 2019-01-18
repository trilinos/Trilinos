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

C $Log: getdst.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:01:37  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:05  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GETDST( NODE1, NODE2, XN, YN, ZN, DIST)
C=======================================================================

C   --*** GETDST *** (MESH) GETS DISTANCE BETWEEN TWO NODES
C   --   Written by RAY J. Meyers 21 June, 1990
C   --
C   -- calculates the distance between two nodes
C   --    (written because of referencing problem using a(xn) directly )
C   --
C   --Parameters:
C   --
C   --   NODE1 - IN - THE INDEX OF FIRST NODE
C   --   NODE2 - IN - THE INDEX OF FIRST NODE
C   --   XN - IN - THE ROTATED, DEFORMED X COORDINATES
C   --   YN - IN - THE ROTATED, DEFORMED Y COORDINATES
C   --   ZN - IN - THE ROTATED, DEFORMED Z COORDINATES
C   --   DIST - OUT - THE DISTANCE BETWEEN NODE1 AND NODE2

C=======================================================================

      REAL XN(*), YN(*), ZN(*)

      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      XDIST = XN(NODE2) - XN(NODE1)
      YDIST = YN(NODE2) - YN(NODE1)
      IF(IS3DIM) THEN
         ZDIST = ZN(NODE2) - ZN(NODE1)
      ELSE
         ZDIST = 0
      END IF

      DIST = SQRT( XDIST*XDIST + YDIST*YDIST + ZDIST*ZDIST)

      RETURN
      END
