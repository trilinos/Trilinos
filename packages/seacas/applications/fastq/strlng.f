C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: strlng.f,v 1.1 1990/11/30 11:16:45 gdsjaar Exp $
C $Log: strlng.f,v $
C Revision 1.1  1990/11/30 11:16:45  gdsjaar
C Initial revision
C
C
      SUBROUTINE STRLNG (STRING, LEN)
C***********************************************************************
C
C  SUBROUTINE STRLNG = FINDS THE NO. OF NONBLANK STRING CHARACTERS
C
C***********************************************************************
C
      CHARACTER * ( * ) STRING
      CALL STRIPB (STRING, ILEFT, LEN)
      IF (LEN .EQ. 0)LEN = 1
      RETURN
      END
