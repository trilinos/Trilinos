C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C    
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C    
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    

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
