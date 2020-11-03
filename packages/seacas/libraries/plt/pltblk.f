C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      BLOCK DATA PLTBLK
      REAL SAVLEN
      INTEGER IDSHSV
      COMMON /PLTSTY/SAVLEN,IDSHSV
      COMMON /PSAVE/TDEVP(5,10),TTEXTP(40,10),TVECTP(5,10),
     *       TGRAPH(100,10),TMAPP(11,10),IPOPD,IPOPT,IPOPV,IPOPG,IPOPM
      COMMON /MPSTCK/SVMAP(195,10),MAPDEP
      DATA IPOPD/0/,IPOPT/0/,IPOPV/0/,IPOPG/0/,IPOPM/0/
      DATA SAVLEN/0./,IDSHSV/0/
      DATA MAPDEP/0/
      end
