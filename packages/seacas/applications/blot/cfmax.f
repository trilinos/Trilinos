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

C $Log: cfmax.f,v $
C Revision 1.2  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 19:55:20  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:47:55  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CFMAX (VARNP, NLNKF, LINKF1, FEMIN, FEMAX)
C=======================================================================

C   --*** CFMAX *** (DETOUR) Compute min/max values for face
C   --   Written by Amy Gilkey - revised 03/14/88
C   --
C   --CFMAX computes the minimum and maximum value of the nodal variable
C   --for the face.
C   --
C   --Parameters:
C   --   VARNP - IN - the contour function values
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   FEMIN, FEMAX - OUT - minimum and maximum value

      REAL VARNP(*)
      INTEGER LINKF1(NLNKF)

C   --Compute the minimum and maximum values for the face

      FEMAX = VARNP(LINKF1(1))
      FEMIN = VARNP(LINKF1(1))
      DO 100 K = 2, NLNKF
         FEMAX = MAX (FEMAX, VARNP(LINKF1(K)))
         FEMIN = MIN (FEMIN, VARNP(LINKF1(K)))
  100 CONTINUE

      RETURN
      END
