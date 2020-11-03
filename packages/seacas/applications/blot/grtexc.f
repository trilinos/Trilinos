C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRTEXC (DX, DY, STRING)
C=======================================================================

C   --*** GRTEXC *** (GRPLIB) Write centered text (PLT)
C   --   Written by Amy Gilkey - revised 03/22/88
C   --
C   --GRTEXC writes a software or hardware character string centered on
C   --a location.
C   --
C   --Parameters:
C   --   DX, DY - IN - the horizontal and vertical string location
C   --      (in device coordinates)
C   --   STRING - IN - the string to be written, may be truncated
C   --
C   --Common Variables:
C   --   Uses ICURDV, SOFTCH of /GRPCOM/

C   --Routines Called:
C   --   PLTXSL - (PLTLIB) Find the software string length
C   --   PLTXTC1 - (PLTLIB) Display a centered software string
C   --   PLTXTH - (PLTLIB) Display a hardware string
C   --   PLTXTL - (PLTLIB) Find the hardware string length
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
         CALL PLTXSL (STRING(:LSTR), SLEN)
         XLEFT = DX
         IF (XLEFT*2 .LT. SLEN) XLEFT = 0.5*SLEN
         CALL PLTXTC1 (XLEFT, DY, STRING(:LSTR))
      ELSE
         CALL PLTXHL (STRING(:LSTR), SLEN)
         XLEFT = DX - .5*SLEN
         IF (XLEFT .LT. 0.0) XLEFT = 0.0000
         CALL PLTXTH (XLEFT, DY, STRING(:LSTR))
      END IF

      RETURN
      END
