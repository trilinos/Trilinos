C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

C $Id: getcav.f,v 1.1 1991/02/21 15:43:18 gdsjaar Exp $
C $Log: getcav.f,v $
C Revision 1.1  1991/02/21 15:43:18  gdsjaar
C Initial revision
C
      SUBROUTINE GETCAV (ERROR, IDESS, NUMESS)
      DIMENSION IDESS(*)
      include 'nu_cav.blk'
      LOGICAL ERROR
      CHARACTER*80 STRA

      PARAMETER (MXFLD = 12)
      DIMENSION RV(MXFLD), KV(MXFLD)
      CHARACTER*32 CV(MXFLD), PRMPT
      PRMPT = '  Cavity Side Set > '
      MAXF = MIN (MXFLD, MAXCAV)
      IF (ICAV(1) .EQ. 0 .OR. NUMCAV .EQ. 0) THEN
         CALL FREFLD (0, 0, PRMPT(:LENSTR(PRMPT)+1), MAXF, IOS,
     *      NUMCAV, KV, CV, ICAV, RV)
      END IF
C
      ERROR = .FALSE.
      DO 10 NCAV = 1, NUMCAV
         IFND(NCAV) = LOCINT (ICAV(NCAV), NUMESS, IDESS)
         IF (IFND(NCAV) .EQ. 0) THEN
            WRITE (STRA, 20) ICAV(NCAV)
            CALL SQZSTR (STRA, LSTR)
            CALL PRTERR ('ERROR', STRA(:LSTR))
            ERROR = .TRUE.
         END IF
   10 CONTINUE
   20 FORMAT (' Cavity Flag ',I5,' not found. ')
C
      RETURN
      END
