C $Id: disclm.f,v 1.1 1990/11/30 11:06:05 gdsjaar Exp $
C $Log: disclm.f,v $
C Revision 1.1  1990/11/30 11:06:05  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]DISCLM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DISCLM (NCOLS)
C***********************************************************************
C
C  SUBROUTINE DISCLM = PRINTS THE SANDIA DISCLAIMER
C
C***********************************************************************
C
      CHARACTER*29 BLANK
C
      DATA BLANK/' '/
C
      NSHIFT = MAX ( (NCOLS-76)/2, 1) + 1
      NSHIFT = MIN (NSHIFT, 29)
      WRITE (*, 10000) (BLANK (1:NSHIFT), I = 1, 9)
      WRITE (*, 10010) (BLANK (1:NSHIFT), I = 1, 5)
      RETURN
C
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
C
      END
