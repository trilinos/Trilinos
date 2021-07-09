C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRTEXT (DX, DY, STRING)
C=======================================================================

C   --*** GRTEXT *** (GRPLIB) Write text (PLT)
C   --   Written by Amy Gilkey - revised 02/14/86
C   --
C   --GRTEXT writes a software or hardware character string at a location
C   --(left-justified).
C   --
C   --Parameters:
C   --   DX, DY - IN - the horizontal and vertical string location
C   --      (in device coordinates)
C   --   STRING - IN - the string to be written, may be truncated
C   --
C   --Common Variables:
C   --   Uses ICURDV, SOFTCH of /GRPCOM/

C   --Routines Called:
C   --   PLTXTH - (PLTLIB) Display a hardware string
C   --   PLTXTS - (PLTLIB) Display a software string
C   --   LENSTR - (STRLIB) Find string length

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL

      REAL DX, DY
      CHARACTER*(*) STRING

      LSTR = LENSTR(STRING)
      IF (STRING(LSTR:LSTR) .EQ. ' ') RETURN

      IF (SOFTCH(ICURDV)) THEN
         CALL PLTXTS (DX, DY, STRING(:LSTR))
      ELSE
         CALL PLTXTH (DX, DY, STRING(:LSTR))
      END IF

      RETURN
      END
