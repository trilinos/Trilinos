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

C $Log: nummod.f,v $
C Revision 1.2  2009/03/25 12:36:46  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:06:16  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:54:18  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      INTEGER FUNCTION NUMMOD (MODE, MODSUB, MATMOD, MATSUB)
C=======================================================================

C   --*** NUMMOD *** (MESH) Count the number of matching modes
C   --   Written by Amy Gilkey - revised 10/05/87
C   --
C   --NUMMOD returns the number of defined modes that match both mode
C   --and sub-mode (if given).
C   --
C   --Parameters:
C   --   MODE - the modes
C   --   MODSUB - the sub-modes
C   --   MATMOD - IN - the mode to match
C   --   MATSUB - IN - the sub-mode to match (if not ' ')

      CHARACTER*(*) MODE(4), MODSUB(4), MATMOD, MATSUB

      NUMMOD = 0

      IF (MATSUB .EQ. ' ') THEN
         DO 100 IVIEW = 1, 4
            IF (MODE(IVIEW) .EQ. MATMOD) NUMMOD = NUMMOD + 1
  100    CONTINUE
      ELSE
         DO 110 IVIEW = 1, 4
            IF ((MODE(IVIEW) .EQ. MATMOD)
     &         .AND. (MODSUB(IVIEW) .EQ. MATSUB))
     &         NUMMOD = NUMMOD + 1
  110    CONTINUE
      END IF

      RETURN
      END
