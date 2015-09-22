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

C $Log: vec.f,v $
C Revision 1.2  2009/03/25 12:36:49  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:17:27  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:59:29  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE VEC (IS3DIM, X0, Y0, Z0, XVAR, YVAR, ZVAR,
     &   VECSCL, VWSCL)
C=======================================================================

C   --*** VEC *** (DETOUR) Plot vector
C   --   Written by Amy Gilkey - revised 04/16/85
C   --
C   --VEC displays a vector starting at the given coordinates.  The vector
C   --represents the given X, Y, and Z values.
C   --
C   --Parameters:
C   --   IS3DIM - IN - true iff 3D versus 2D
C   --   X0, Y0, Z0 - IN - the vector coordinates
C   --   XVAR, YVAR, ZVAR - IN - the vector components
C   --   VECSCL - IN - the vector scale factor
C   --   VWSCL - IN - 1.0 if single view, 0.5 if multiple

      LOGICAL IS3DIM

      LOGICAL EXISTS

      EXISTS (M) = (MOD(M,2) .NE. 0)

      IF ((XVAR .EQ. 0.0) .AND. (YVAR .EQ. 0.0)) GOTO 100
      X1 = X0 + XVAR * VECSCL
      Y1 = Y0 + YVAR * VECSCL
      Z1 = Z0 + ZVAR * VECSCL
      CALL MP2VC (1, X0, Y0, X1, Y1, DX0, DY0, DX1, DY1, MASK)
      RAT = .0075 * VWSCL
      IF (EXISTS (MASK)) THEN
         CALL PLTARR (DX0, DY0, DX1, DY1, .5, RAT)
      END IF

  100 CONTINUE
      RETURN
      END
