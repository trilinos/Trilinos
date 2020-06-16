C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details
c $Id: myamod.f,v 1.3 2008/03/14 13:45:28 gdsjaar Exp $
C=======================================================================
      REAL FUNCTION MYAMOD (PARM1, PARM2)
C=======================================================================
c

C   --*** MYAMOD *** (ALGEBRA) amod function
C   --
C   --   PARM1, PARM2 - IN - the maximum parameters

      MYAMOD = MOD(PARM1, PARM2)

      RETURN
      END
