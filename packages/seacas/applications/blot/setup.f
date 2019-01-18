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

C $Log: setup.f,v $
C Revision 1.2  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:11:45  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:57:16  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE SETUP (MSHDEF, ZMMESH)
C=======================================================================

C   --*** SETUP *** (MESH) Setup for plot set
C   --   Written by Amy Gilkey - revised 05/26/88
C   --   D. P. Flanagan, 07/27/82
C   --
C   --SETUP sets graphics parameters common to the plot sequence.  It
C   --computes the display area and window for each view.  It also computes
C   --and sets the axis lengths and limits (with "nice" tick intervals).
C   --The only information that is returned about the axis is the
C   --starting locations and the label size factor (through COMMON).
C   --
C   --Parameters:
C   --   MSHDEF - IN - the display modes for all views (none/defined)
C   --      (as in /MSHOPT/)
C   --   ZMMESH - IN - the zoomed mesh coordinates
C   --
C   --Common Variables:
C   --   Uses XISSYM, YISSYM, XAXSYM, YAXSYM, LFTSYM, BOTSYM of /VIEWS/
C   --   Uses ZMMESH of /MSHLIM/
C   --   Uses DVIEW0 of /LAYOUT/
C   --   Sets DVIEW, WVIEW of /LAYOUT/
C   --   Uses DTW, VWSCL of /DEVDAT/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      COMMON /VIEWS/  MULTIM,
     &   XISSYM, YISSYM, XAXSYM, YAXSYM, LFTSYM, BOTSYM
      LOGICAL MULTIM, XISSYM, YISSYM, LFTSYM, BOTSYM
      COMMON /LAYOUT/ CHLSIZ, DBORD0(KTOP), DVIEW0(KTOP)
      COMMON /LAYOUD/ DVIEW(KTOP,4), WVIEW(KTOP,4)
      COMMON /DEVDAT/ DTW, VWSCL

      CHARACTER*(*) MSHDEF(4)
      REAL ZMMESH(KTOP)

      INTEGER NDEFVW, IXVW
      LOGICAL RGTVW, TOPVW
      LOGICAL XISDIV, YISDIV

      RGTVW(I) = (I .EQ. 2) .OR. (I .EQ. 4)
      TOPVW(I) = (I .EQ. 1) .OR. (I .EQ. 2)

      XISDIV = (MSHDEF(1) .NE. 'NONE')
      YISDIV = (MSHDEF(4) .NE. 'NONE')

C   --WXRNG and WYRNG are the window limits for a single view
      WXRNG = ZMMESH(KRGT) - ZMMESH(KLFT)
      WYRNG = ZMMESH(KTOP) - ZMMESH(KBOT)

C   --Set view device boundaries

      IF (XISDIV .EQV. YISDIV) THEN
         RAT = WXRNG / WYRNG
      ELSE IF (XISDIV) THEN
         RAT = (2 * WXRNG) / WYRNG
      ELSE IF (YISDIV) THEN
         RAT = WXRNG / (2 * WYRNG)
      END IF

C   --DXRNG2 and DYRNG2 are the device limits for all defined views / 2
      DXRNG2 = 0.5 * (DVIEW0(KRGT) - DVIEW0(KLFT))
      DYRNG2 = 0.5 * (DVIEW0(KTOP) - DVIEW0(KBOT))
      IF (RAT .LT. 0.99) THEN
         DXRNG2 = DXRNG2 * RAT
      ELSE IF (RAT .GT. 1.01) THEN
         DYRNG2 = DYRNG2 / RAT
      END IF

      DXMID = 0.5 * (DVIEW0(KRGT) + DVIEW0(KLFT))
      DYMID = 0.5 * (DVIEW0(KTOP) + DVIEW0(KBOT))
      DO 100 IVW = 1, NDEFVW (.TRUE.)
         IVIEW = IXVW (.TRUE., IVW)
         IF (XISDIV) THEN
            IF (RGTVW(IVIEW)) THEN
               DVIEW(KLFT,IVIEW) = DXMID
               DVIEW(KRGT,IVIEW) = DXMID + DXRNG2
            ELSE
               DVIEW(KLFT,IVIEW) = DXMID - DXRNG2
               DVIEW(KRGT,IVIEW) = DXMID
            END IF
         ELSE
            DVIEW(KLFT,IVIEW) = DXMID - DXRNG2
            DVIEW(KRGT,IVIEW) = DXMID + DXRNG2
         END IF
         IF (YISDIV) THEN
            IF (TOPVW(IVIEW)) THEN
               DVIEW(KBOT,IVIEW) = DYMID
               DVIEW(KTOP,IVIEW) = DYMID + DYRNG2
            ELSE
               DVIEW(KBOT,IVIEW) = DYMID - DYRNG2
               DVIEW(KTOP,IVIEW) = DYMID
            END IF
         ELSE
            DVIEW(KBOT,IVIEW) = DYMID - DYRNG2
            DVIEW(KTOP,IVIEW) = DYMID + DYRNG2
         END IF
  100 CONTINUE

C   --Set window limits

      IF (XISSYM) THEN
         IF (LFTSYM) THEN
            WLFT = XAXSYM
            WRGT = XAXSYM + WXRNG
         ELSE
            WLFT = XAXSYM
            WRGT = XAXSYM - WXRNG
         END IF
      ELSE
         WLFT = ZMMESH(KLFT)
         WRGT = ZMMESH(KRGT)
      END IF
      IF (YISSYM) THEN
         IF (BOTSYM) THEN
            WBOT = YAXSYM
            WTOP = YAXSYM + WYRNG
         ELSE
            WBOT = YAXSYM
            WTOP = YAXSYM - WYRNG
         END IF
      ELSE
         WBOT = ZMMESH(KBOT)
         WTOP = ZMMESH(KTOP)
      END IF

C   --Set view window boundaries

      DO 110 IVW = 1, NDEFVW (.TRUE.)
         IVIEW = IXVW (.TRUE., IVW)
         IF (XISSYM .AND. (.NOT. RGTVW(IVIEW))) THEN
            WVIEW(KLFT,IVIEW) = WRGT
            WVIEW(KRGT,IVIEW) = WLFT
         ELSE
            WVIEW(KLFT,IVIEW) = WLFT
            WVIEW(KRGT,IVIEW) = WRGT
         END IF
         IF (YISSYM .AND. (.NOT. TOPVW(IVIEW))) THEN
            WVIEW(KBOT,IVIEW) = WTOP
            WVIEW(KTOP,IVIEW) = WBOT
         ELSE
            WVIEW(KBOT,IVIEW) = WBOT
            WVIEW(KTOP,IVIEW) = WTOP
         END IF
  110 CONTINUE

C   --Set window to device units ratio

      DTW = MIN (ABS (WRGT - WLFT), ABS (WTOP - WBOT))
     &   / (DVIEW0(KRGT) - DVIEW0(KLFT))

      IF (XISDIV .OR. YISDIV) THEN
         VWSCL = 0.5
      ELSE
         VWSCL = 1.0
      END IF

      RETURN
      END
