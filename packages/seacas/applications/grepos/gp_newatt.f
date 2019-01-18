C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
C                                                 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C=======================================================================
      SUBROUTINE NEWATT (IDLST, ID, NUM, IDATT, NUMATR, ATTNAM, NEWNAM)
C=======================================================================
C
      include 'gp_namlen.blk'
      INTEGER IDLST(*)
      INTEGER NUMATR(*)
      CHARACTER*(maxnam) ATTNAM(*), NEWNAM
      CHARACTER*1024 STRING
      CHARACTER*16 STRA

      STRA = 'Element Block'

      IF (NUM .LE. 0) RETURN

C ... Determine location of ID to be changed

      IMAT = LOCINT (ID, NUM, IDLST)
      IF (IMAT .EQ. 0) THEN
        WRITE (STRING, 90) STRA, ID
 90     FORMAT (A,1X,I5,' does not exist')
        CALL SQZSTR (STRING, LSTR)
        CALL PRTERR ('ERROR', STRING(:LSTR))
        RETURN
      END IF

C ... Determine beginning of attribute names for this block
      IOFF = 0
      DO I=1, IMAT-1
        IOFF = IOFF + NUMATR(I)
      END DO
      
      ATTNAM(IOFF + IDATT) = NEWNAM
      ATTNAM(IOFF + IDATT) = NEWNAM 
      WRITE (STRING, 100) IDATT, STRA, ID, NEWNAM
 100  FORMAT ('Name of Attribute ', I9, ' on ', A,1X,I9,
     *  ' changed to ',A)
         CALL SQZSTR (STRING, LSTR)
         CALL PRTERR ('CMDSPEC', STRING(:LSTR))

      RETURN
      END
