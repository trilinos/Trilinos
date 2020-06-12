C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: sflush.f,v 1.1 1990/11/30 11:15:36 gdsjaar Exp $
C $Log: sflush.f,v $
C Revision 1.1  1990/11/30 11:15:36  gdsjaar
C Initial revision
C
C
      SUBROUTINE SFLUSH
C***********************************************************************
C
C  SUBROUTINE SFLUSH = SCREEN FLUSH (DUMPS GRAPHICS BUFFER TO THE SCREEN)
C
C***********************************************************************
C
      CALL PLTFLU
      RETURN
C
      END
