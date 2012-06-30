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

C $Log: sqrlim.f,v $
C Revision 1.2  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:15:17  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:58:37  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SQRLIM (RDMESH, SQMESH)
C=======================================================================

C   --*** SQRLIM *** (MESH) Expand 2D mesh limits to square
C   --   Written by Amy Gilkey - revised 05/23/86
C   --
C   --SQRLIM expands the mesh limits to form a square.
C   --
C   --Parameters:
C   --   RDMESH - IN - the mesh limits
C   --      (left, right, bottom, top)
C   --   SQMESH - OUT - the square mesh limits (may be RDMESH)
C   --      (left, right, bottom, top)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      REAL RDMESH(KTOP), SQMESH(KTOP)

      RNG2 = .5 * MAX (RDMESH(KRGT) - RDMESH(KLFT)
     &   ,             RDMESH(KTOP) - RDMESH(KBOT))
      XCEN = .5 * (RDMESH(KRGT) + RDMESH(KLFT))
      YCEN = .5 * (RDMESH(KTOP) + RDMESH(KBOT))
      SQMESH(KLFT) = XCEN - RNG2
      SQMESH(KRGT) = XCEN + RNG2
      SQMESH(KBOT) = YCEN - RNG2
      SQMESH(KTOP) = YCEN + RNG2

      RETURN
      END
