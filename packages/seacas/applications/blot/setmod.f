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

C=======================================================================
      SUBROUTINE SETMOD (IVIEW, MMOD, MTYP)
C=======================================================================

C   --*** SETMOD *** (DETOUR) Set display mode and type
C   --   Written by Amy Gilkey - revised 04/06/88
C   --
C   --SETMOD sets the display mode and mode type of one or all views.
C   --Counts the number of variables needed, if changed
C   --
C   --Parameters:
C   --   IVIEW - IN - the view to be set, 0 for all
C   --   MMOD - IN - the display mode to be set
C   --   MTYP - IN - the display mode type to be set
C   --
C   --Common Variables:
C   --   Uses MSHDEF of /MSHOPT/
C   --   Sets MODDET, MODTYP, NNDVAR, NEDVAR of /DETOPT/

      include 'mshopt.blk'
      COMMON /DETOPT/ IDTVAR(4), NNDVAR, NEDVAR
      COMMON /DETOPC/ MODDET(4), MODTYP(4)
      CHARACTER*8 MODDET, MODTYP

      CHARACTER*(*) MMOD, MTYP

      IF (IVIEW .GE. 1) THEN
         IF (MMOD .EQ. 'NONE') THEN
            MODDET(IVIEW) = 'NONE'
            MODTYP(IVIEW) = ' '
         ELSE
            MODDET(IVIEW) = MMOD
            MODTYP(IVIEW) = MTYP
         END IF
      ELSE
         IF (MMOD .EQ. 'NONE') THEN
            DO 100 I = 1, 4
               MODDET(I) = 'NONE'
               MODTYP(I) = ' '
  100       CONTINUE
         ELSE
            DO 110 I = 1, 4
               IF ((MSHDEF(I) .EQ. 'NONE')
     &            .OR. (MSHDEF(I) .EQ. 'EMPTY')) THEN
                  MODDET(I) = 'NONE'
                  MODTYP(I) = ' '
               ELSE
                  MODDET(I) = MMOD
                  MODTYP(I) = MTYP
               END IF
  110       CONTINUE
         END IF
      END IF

C   --Calculate the number of variables needed
      CALL CNTVAR (MODDET, MODTYP, IDTVAR, NNDVAR, NEDVAR)

      RETURN
      END
