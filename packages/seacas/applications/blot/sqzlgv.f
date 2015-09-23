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

C $Log: sqzlgv.f,v $
C Revision 1.2  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:15:34  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:58:41  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SQZLGV (NPTIMS, IPTIMS, WHOTIM, VALIN, NPTOUT, VALOUT)
C=======================================================================

C   --*** SQZLGV *** (TPLOT) Compress time-dependent values
C   --   Written by Amy Gilkey - revised 11/06/87
C   --
C   --SQZLGV compresses time-dependent curve values so that only the
C   --values for whole time steps are in the curve.
C   --
C   --Parameters:
C   --   NPTIMS - IN - the number of selected time steps
C   --   IPTIMS - IN - the selected time steps
C   --   WHOTIM - IN - true iff whole (versus history) time step
C   --   VALIN - IN - the input data
C   --   NPTOUT - OUT - the number of selected whole time steps
C   --   VALOUT - OUT - the output (compressed) data

      INTEGER IPTIMS(*)
      LOGICAL WHOTIM(*)
      REAL VALIN(*)
      REAL VALOUT(*)

      NPTOUT = 0
      DO 100 I = 1, NPTIMS
         IF (WHOTIM(IPTIMS(I))) THEN
            NPTOUT = NPTOUT + 1
            VALOUT(NPTOUT) = VALIN(I)
         END IF
  100 CONTINUE

      RETURN
      END
