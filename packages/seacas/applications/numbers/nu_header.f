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
      SUBROUTINE HEADER (NDIM, TITLE, NUMEL, NUMNP, AXI, GENFIL)
C
      include 'nu_io.blk'
      CHARACTER*16 FORM(3)
      CHARACTER*80 TITLE
      CHARACTER*(*) GENFIL
      LOGICAL AXI, FIRST
      CHARACTER*6 DIMEN(3)

      DATA DIMEN/'One',   'Two',   'Three'/
      DATA FORM /', Planar',', Axisymmetric',' '/
      DATA FIRST /.TRUE./
C
      IO = IHARD
      IF (FIRST) THEN
          WRITE (IO, 10) GENFIL(:LENSTR(GENFIL))
   10     FORMAT (5X,'Genesis: ',A/)
          WRITE (IO, 20) TITLE(:LENSTR(TITLE))
   20     FORMAT (5X,'Title:   ',A/)
          ILAB = 2
          IF (NDIM .EQ. 2 .AND. .NOT. AXI) ILAB = 1
          IF (NDIM .EQ. 3) ILAB = 3
          WRITE (IO, 30) NUMNP,
     *        NUMEL, DIMEN(NDIM)(:LENSTR(DIMEN(NDIM))),
     *        FORM(ILAB)(:LENSTR(FORM(ILAB)))
   30     FORMAT (5X,'Number of Nodes:    ',I10,/
     *        5X,'Number of Elements: ',I10/
     *        5X,A,'-Dimensional Mesh',A)
          FIRST = .FALSE.
      END IF
      RETURN
      END
