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

C $Log: grview.f,v $
C Revision 1.2  2009/03/25 12:36:45  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:58  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:52:05  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GRVIEW (ASPECT, DVIEW0, DVIEW)
C=======================================================================

C   --*** GRVIEW *** (GRPLIB) Set graphics window
C   --   Written by Amy Gilkey - revised 02/20/86
C   --
C   --GRVIEW sets the device coordinates for the window given the X:Y
C   --axis length ratio.
C   --
C   --Parameters:
C   --   ASPECT - IN - the ratio of the X axis length to the Y axis length
C   --      = (X axis length / Y axis length) {in window units}
C   --      / (X axis length / Y axis length) {in device units}
C   --   DVIEW0 - IN - the maximum window (left, right, bottom, top)
C   --      (in device coordinates)
C   --   DVIEW - OUT - the actual window (left, right, bottom, top)
C   --      (in device coordinates)

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      REAL ASPECT
      REAL DVIEW0(KTOP), DVIEW(KTOP)

      IF (ASPECT .GE. 1.0) THEN
         DVIEW(KLFT) = DVIEW0(KLFT)
         DVIEW(KRGT) = DVIEW0(KRGT)
      ELSE
         OS = .5 * (1.0 - ASPECT) * (DVIEW0(KRGT) - DVIEW0(KLFT))
         DVIEW(KLFT) = DVIEW0(KLFT) + OS
         DVIEW(KRGT) = DVIEW0(KRGT) - OS
      END IF
      IF (ASPECT .LE. 1.0) THEN
         DVIEW(KBOT) = DVIEW0(KBOT)
         DVIEW(KTOP) = DVIEW0(KTOP)
      ELSE
         OS = .5 * (1.0 - 1.0/ASPECT) * (DVIEW0(KTOP) - DVIEW0(KBOT))
         DVIEW(KBOT) = DVIEW0(KBOT) + OS
         DVIEW(KTOP) = DVIEW0(KTOP) - OS
      END IF

      RETURN
      END
