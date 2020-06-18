C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltcpy.f,v 1.1 1993/07/16 16:47:52 gdsjaar Exp $
C $Log: pltcpy.f,v $
C Revision 1.1  1993/07/16 16:47:52  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE PLTCPY
      INTEGER SUPPO

      CALL PLTFLU
      CALL VDIQES(100,SUPPO)
      IF (SUPPO.NE.0) THEN
         CALL VDESCP(100,0)
      END IF

      RETURN

      END
