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

C $Log: grpbeg.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:39  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:48  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GRPBEG
C=======================================================================

C   --*** GRPBEG *** (GRPLIB) Begin a plot by setting defaults (PLT)
C   --   Written by Amy Gilkey - revised 03/01/88
C   --
C   --GRPBEG begins a plot by enabling the interrupt flag, clearing the
C   --screen, and setting standard graphics parameters (such as character
C   --sizes).  The graphics parameters are saved the first time through
C   --this routine and reset to the saved values on following calls.

C   --Routines Called:
C   --   CPUIFC - (PLTLIB) Initialize interrupt
C   --   PLTBGN - (PLTLIB) Erase display surface
C   --   PLTGTG - (PLTLIB) Get graph parameter (see PLTSTG)
C   --   PLTSTG - (PLTLIB) Set graph parameter
C   --      22, 47 = (KXNUMS, KYNUMS) X, Y axis number size
C   --      23, 48 = (KXLABS, KYLABS) X, Y axis label size
C   --   GRCOLR - (GRPLIB) Set color
C   --   GRSNAP - (GRPLIB) Handle device frame snapping

      PARAMETER (KXNUMS=22, KXLABS=23, KYNUMS=47, KYLABS=48)
      PARAMETER (KTICSZ=33)

      LOGICAL CPUIFC
      LOGICAL LDUM, PLTGTG, PLTSTG

      LOGICAL FIRST
      SAVE FIRST
      SAVE SZNUM, SZLAB, SZTIC

      DATA FIRST / .TRUE. /

C   --Enable interrupt flag
      LDUM = CPUIFC (.TRUE.)

C   --Erase display surface
      CALL PLTBGN

C   --Start plot for snapping
      CALL GRSNAP ('START', 0)

C   --Save standard graphics parameters
      IF (FIRST) THEN
         LDUM = PLTGTG (KXNUMS, SZNUM)
         LDUM = PLTGTG (KXLABS, SZLAB)
         LDUM = PLTGTG (KTICSZ, SZTIC)
         FIRST = .FALSE.
      END IF

C   --Set up label/numbering size
      LDUM = PLTSTG (KXNUMS, SZNUM)
      LDUM = PLTSTG (KXLABS, SZLAB)
      LDUM = PLTSTG (KYNUMS, SZNUM)
      LDUM = PLTSTG (KYLABS, SZLAB)
      LDUM = PLTSTG (KTICSZ, SZTIC)

C   --Set to standard color
      CALL GRCOLR (0)

      RETURN
      END
