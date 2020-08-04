C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE DISCLM (NCOLS)
C***********************************************************************

C  SUBROUTINE DISCLM = PRINTS THE SANDIA DISCLAIMER

C***********************************************************************

      CHARACTER*29 BLANK

      DATA BLANK/' '/

      NSHIFT = MAX ( (NCOLS-76)/2, 1) + 1
      NSHIFT = MIN (NSHIFT, 29)
      WRITE (*, 10000) (BLANK (1:NSHIFT), I = 1, 9)
      WRITE (*, 10010) (BLANK (1:NSHIFT), I = 1, 5)
      RETURN

10000 FORMAT (' ', A,
     &   '**************************************',
     &   '**************************************'//, A,
     &   '                  ISSUED BY SANDIA NAT',
     &   'IONAL LABORATORIES, '/, A,
     &   '                         A PRIME CONTR',
     &   'ACTOR TO THE'/, A,
     &   '                     UNITED STATES DEP',
     &   'ARTMENT OF ENERGY'//, A,
     &   '**************************************',
     &   '**************************************'//, A,
     &   'THIS CODE WAS PREPARED IN THE COURSE O',
     &   'F WORK SPONSORED BY THE UNITED STATES'/, A,
     &   'GOVERNMENT.  NEITHER THE UNITED STATES',
     &   ',  NOR THE UNITED STATES DEPARTMENT OF'/, A,
     &   'ENERGY,  NOR THE UNITED STATES NUCLEAR ',
     &   'REGULATORY COMMISSION,  NOR ANY OF'/, A,
     &   'THEIR EMPLOYEES,  NOR ANY OF THEIR CONT',
     &   'RACTORS,  SUBCONTRACTORS,  OR THEIR')
10010 FORMAT (A,
     &   'EMPLOYEES,  MAKES ANY WARRANTY,  EXPRESS',
     &   ' OR IMPLIED,  OR ASSUMES ANY LEGAL'/A,
     &   'LIABILITY OR RESPONSIBILITY FOR THE AC',
     &   'CURACY,  COMPLETENESS OR USEFULNESS OF'/, A,
     &   'ANY INFORMATION,  APPARATUS,  PRODUCT OR',
     &   ' PROCESS DISCLOSED,  OR REPRESENTS THAT'/, A,
     &   'ITS USE WOULD NOT INFRINGE PRIVATELY O',
     &   'WNED RIGHTS.'//, A,
     &   '**************************************',
     &   '**************************************'/)

      END
