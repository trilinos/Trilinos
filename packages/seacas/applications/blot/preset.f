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

C $Log: preset.f,v $
C Revision 1.2  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:07:29  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:55:00  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE PRESET (LINPLT, QAPLT, DBORD0, DVIEW0)
C=======================================================================

C   --*** PRESET *** (BLOT) Initializes special graphics parameters
C   --   Written by Amy Gilkey - revised 02/20/87
C   --
C   --PRESET graphics parameters specific to the plot format:
C   --   o axes (set larger for line plots, same for mesh plots)
C   --   o tick marks (set larger for line plots, smaller for mesh plots)
C   --   o grid lines (set smaller for line plots)
C   --   o normal lines (set larger for line plots, same for mesh plots
C   --   o curve lines (set larger for line plots)
C   --
C   --It also sets the display area (dependent on QAPLT).
C   --
C   --Parameters:
C   --   LINPLT - IN - true if line versus mesh plot
C   --   QAPLT - IN - if true, set up as standard QA plot, else set up to
C   --      plot the entire screen (minus border)
C   --   DBORD0 - IN - the display and label area boundary
C   --   DVIEW0 - OUT - the display area boundary

C   --Routines Called:
C   --   PLTGTG - (PLTLIB) Get graph parameters (see PLTSTG)
C   --   PLTSTG - (PLTLIB) Set graph parameters:
C   --      28 = (KAXESZ) axis size
C   --      33 = (KTICSZ) tick-mark size
C   --      34 = (KGRISZ) grid-line size
C   --      29 = (KCRVSZ) line size
C   --      15 = (KONGRI) no grid
C   --      37 = (KZLINE) delete dashed line at 0

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4)

      PARAMETER (KAXESZ=28, KTICSZ=33, KGRISZ=34, KCRVSZ=29)
      PARAMETER (KZLINE=37)
      PARAMETER (KONGRI=15)

      LOGICAL LINPLT, QAPLT
      REAL DBORD0(KTOP), DVIEW0(KTOP)

      LOGICAL LDUM, PLTGTG, PLTSTG

      LOGICAL FIRST
      SAVE FIRST
      SAVE SZAXES, SZTICK, SZGRID, SZCURV

      DATA FIRST / .TRUE. /

      IF (QAPLT) THEN
         DVIEW0(KLFT) = DBORD0(KLFT) + 0.09
         DVIEW0(KRGT) = DVIEW0(KLFT) + 0.62
         DVIEW0(KBOT) = DBORD0(KBOT) + 0.10
         DVIEW0(KTOP) = DVIEW0(KBOT) + 0.62
      ELSE
         BORDER = 0.1
         DVIEW0(KLFT) = DBORD0(KLFT) + BORDER
         DVIEW0(KRGT) = DBORD0(KRGT) - BORDER
         DVIEW0(KBOT) = DBORD0(KBOT) + BORDER
         DVIEW0(KTOP) = DBORD0(KTOP) - BORDER
      END IF

C   --Save default device settings

      IF (FIRST) THEN
         LDUM = PLTGTG (KAXESZ, SZAXES)
         LDUM = PLTGTG (KTICSZ, SZTICK)
         LDUM = PLTGTG (KGRISZ, SZGRID)
         LDUM = PLTGTG (KCRVSZ, SZCURV)
         FIRST = .FALSE.
      END IF

      IF (LINPLT) THEN
         LDUM = PLTSTG (KAXESZ, 1.75*SZAXES)
         LDUM = PLTSTG (KTICSZ, 1.75*SZTICK)
         LDUM = PLTSTG (KGRISZ, 0.50*SZGRID)
         LDUM = PLTSTG (KCRVSZ, 1.50*SZCURV)
      ELSE
         LDUM = PLTSTG (KAXESZ, SZAXES)
         LDUM = PLTSTG (KTICSZ, 0.50*SZTICK)
      END IF

C   --Set no line at zero
      LDUM = PLTSTG (KZLINE, 0.0)

C   --Set no grid
      LDUM = PLTSTG (KONGRI, 0.0)

      RETURN
      END
