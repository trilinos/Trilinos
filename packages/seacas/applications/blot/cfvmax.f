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

C $Log: cfvmax.f,v $
C Revision 1.3  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1997/03/20 21:23:52  caforsy
C Update Imakefile for Imake 6.1
C
C Revision 1.1  1994/04/07 19:55:26  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:47:58  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE CFVMAX (NUMFAC, VARFAC, NMIN, NMAX, FMIN, FMAX)
C=======================================================================

C   --*** CFVMAX *** (DETOUR) Calculate face variable min/max and count number
C   --   Written by Amy Gilkey - revised 03/09/88
C   --
C   --CFVMAX calculates the minimum and maximum value of a face variable
C   --and counts the number of occurances of the minimum and maximum.
C   --
C   --Parameters:
C   --   NUMFAC - IN - the number of surface faces
C   --   VARFAC - IN - the face variable values
C   --   NMIN, NMAX - OUT - the number of variables values matching the
C   --      minimum and the maximum
C   --   FMIN, FMAX - OUT - the face variable minimums and maximums

      REAL VARFAC(*)

      LOGICAL INIT

      NMIN = 0
      NMAX = 0
      FMIN = 0.0
      FMAX = 0.0
      INIT = .TRUE.
      DO 100 IFAC = 1, NUMFAC
         IF (INIT) THEN
            FMIN = VARFAC(IFAC)
            NMIN = 1
            FMAX = VARFAC(IFAC)
            NMAX = 1
            INIT = .FALSE.
         ELSE
            IF (FMIN .GE. VARFAC(IFAC)) THEN
               IF (FMIN .EQ. VARFAC(IFAC)) THEN
                  NMIN = NMIN + 1
               ELSE
                  FMIN = VARFAC(IFAC)
                  NMIN = 1
               END IF
            END IF
            IF (FMAX .LE. VARFAC(IFAC)) THEN
               IF (FMAX .EQ. VARFAC(IFAC)) THEN
                  NMAX = NMAX + 1
               ELSE
                  FMAX = VARFAC(IFAC)
                  NMAX = 1
               END IF
            END IF
         END IF
  100 CONTINUE

      RETURN
      END

