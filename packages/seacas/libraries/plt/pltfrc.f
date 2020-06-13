C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltfrc.f,v 1.1 1993/07/16 16:48:10 gdsjaar Exp $
C $Log: pltfrc.f,v $
C Revision 1.1  1993/07/16 16:48:10  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      FUNCTION PLTFRC(REALN)

      PLTFRC = REALN - INT(REALN)
      RETURN

      END
