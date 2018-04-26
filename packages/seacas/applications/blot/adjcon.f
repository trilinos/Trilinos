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

C $Log: adjcon.f,v $
C Revision 1.2  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 19:54:29  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:47:26  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE ADJCON (DELCOK)
C=======================================================================

C   --*** ADJCON *** (DETOUR) Adjust the contour parameters
C   --   Written by Amy Gilkey - revised 07/10/87
C   --
C   --ADJCON adjusts the contour parameters to be consistent with each
C   --other.  If CINTOK is false and the parameter DELCOK is true, the
C   --contour maximum (CMAX) is set; if DELCOK is false, the contour
C   --interval (DELC) is set.  If CINTOK is true, the contour minimum and
C   --maximum (CMIN and CMAX) and interval (DELC) are set.
C   --
C   --Parameters:
C   --   DELCOK - IN - true if CMAX is to be set, false if DELC is to be set
C   --
C   --Common Variables:
C   --   Uses CINTOK, LINCON, NCNTR, CMIN, CMAX, DELC, CINTV of /CNTR/
C   --   Sets CMIN, CMAX, DELC of /CNTR/

      include 'cntr.blk'

      LOGICAL DELCOK

      IF (LINCON) THEN
         NC = NCNTR
      ELSE
         NC = NCNTR+1
      END IF

      IF (.NOT. CINTOK) THEN
         IF (DELCOK) THEN
            CMAX = CMIN + (NC-1) * DELC
         ELSE
            DELC = CMAX - CMIN
            IF (NC .GT. 1) DELC = (CMAX - CMIN) / (NC-1)
         END IF

      ELSE
         CMIN = CINTV(NC)
         CMAX = CINTV(NC)
         DELC = CMAX - CMIN
         IF (NC .GT. 1) DELC = (CMAX - CMIN) / (NC-1)
      END IF

      RETURN
      END
