C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION CPUIFC (LDUM)
C=======================================================================
C$Id: cpuifc.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: cpuifc.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:12:19  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:17  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:06  gdsjaar
c Initial revision
c

C   --*** CPUIFC *** Dummy cancel function
C   --   Written by Amy Gilkey - revised 02/11/88
C   --
C   --CPUIFC returns the cancel flag as false.

      LOGICAL LDUM

      CPUIFC = .FALSE.

      RETURN
      END
