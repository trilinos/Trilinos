C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
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
C    * Neither the name of Sandia Corporation nor the names of its
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
C=======================================================================
      SUBROUTINE PRGLOB (OPTION, NOUT, NVARGL, LISGV, NAMEGV, VARGL)
C=======================================================================

C   --*** PRGLOB *** (GROPE) Display current database global variables
C   --
C   --PRGLOB displays the global data for a time step.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NVARGL - IN - the number of global variables
C   --   LISGV - IN - the indices of the selected global variables
C   --   NAMEGV - IN - the names of the global variables
C   --   VARGL - IN - the global variables for the time step

      include 'params.blk'
      CHARACTER*(*) OPTION
      INTEGER LISGV(0:*)
      CHARACTER*(*) NAMEGV(*)
      REAL VARGL(*)
      INTEGER GETPRC, PRTLEN
      CHARACTER*128 FMT1, FMT

      PRTLEN = GETPRC() + 7
      WRITE(FMT1,20) PRTLEN, PRTLEN-7
      CALL SQZSTR(FMT1, LFMT)
      WRITE(FMT, 30) FMT1(:LFMT)

      LNAM = 0
      do i=1, lisgv(0)
        L = lenstr(namegv(lisgv(i)))
        lnam = max(l, lnam)
      end do
      
      IF (NOUT .GT. 0) THEN
        WRITE (NOUT, 10000)
      ELSE
        WRITE (*, 10000)
      END IF

      do 100 i=1, lisgv(0)
        IF (NOUT .GT. 0) THEN
           WRITE (NOUT, FMT) NAMEGV(LISGV(I))(:LNAM), VARGL(LISGV(I))
        ELSE
           WRITE (*, FMT)    NAMEGV(LISGV(I))(:LNAM), VARGL(LISGV(I))
        END IF
 100  continue

      RETURN

 20   FORMAT('1PE',I2.2,'.',I2.2)
 30   FORMAT ('(1X, A, '' = '',', A,')')

10000  FORMAT (/, 1X, 'Global Time Step Variables')
10010  FORMAT (1X, A)
      END
