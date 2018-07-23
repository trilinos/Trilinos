C    Copyright(C) 1988-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

C $Id: selssn.f,v 1.2 1998/03/22 05:34:45 gdsjaar Exp $
C $Log: selssn.f,v $
C Revision 1.2  1998/03/22 05:34:45  gdsjaar
C General cleanp of unused variables. Reordered DATA statements in
C command.f so would compile with f2c.
C
C Revision 1.1.1.1  1991/02/21 15:45:33  gdsjaar
C NUMBERS: Greg Sjaardema, initial Unix release
C
c Revision 1.1  1991/02/21  15:45:32  gdsjaar
c Initial revision
c
C=======================================================================
      SUBROUTINE SELSSN (SELECT, NUMNP, NLIST, LIST,
     *   IDSS, NNSS, IPNSS, LTNSS, NUMSS, NUMSEL)
C=======================================================================
C    IDSS  (NUMSS) SIDE SET IDS
C    NNSS  (NUMSS) SIDE SET NODE    COUNTS
C    IPNSS (NUMSS) SIDE SET NODE    POINTERS
C    LTNSS (LSSNL) SIDE SET NODE    LIST
C
      LOGICAL SELECT(*)
      INTEGER LIST(*), IDSS(*), NNSS(*), IPNSS(*), LTNSS(*)
      CHARACTER*80 STRA

      CALL INILOG (NUMNP, .FALSE., SELECT)

      DO 20 II = 1, NLIST
         IFLG = LOCINT (LIST(II), NUMSS, IDSS)
         IF (IFLG .EQ. 0) THEN
            WRITE (STRA, 50) LIST(II)
            CALL SQZSTR (STRA, LSTR)
            CALL PRTERR ('ERROR', STRA(:LSTR))
         ELSE
            IBEG = IPNSS(IFLG)
            IEND = IBEG + NNSS(IFLG) - 1
            DO 10 I=IBEG, IEND
               SELECT(LTNSS(I)) = .TRUE.
   10       CONTINUE
         END IF
   20 CONTINUE

      NUMSEL = NUMEQL (.TRUE., NUMNP, SELECT)
      WRITE (STRA, 30) NUMSEL
      CALL SQZSTR(STRA, LSTR)
      WRITE (*, 40) STRA(:LSTR)
   30 FORMAT (I10,' nodes selected')
   40 FORMAT (/5X,A)

   50 FORMAT (' Set Flag ',I10,' not found. ')
      RETURN
      END
