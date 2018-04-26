C Copyright(C) 2011-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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
      SUBROUTINE FELCEN (NUMCOL, NEROW, IELROW, IXROW, IELCEN)
C=======================================================================

C   --*** FELCEN *** (GEN3D) Put center elements into row by column array
C   --   Written by Amy Gilkey - revised 04/26/88
C   --
C   --It puts the center elements into a row by column array.
C   --
C   --Parameters:
C   --   NUMCOL - IN - the number of columns in the center blocks
C   --   NEROW - IN - the number of element rows in the center blocks
C   --   IELROW - IN - the element numbers of the rows of center elements
C   --   IXROW - IN - the IELROW index of the starting column for each row
C   --   IELCEN - OUT - the element numbers of the center elements
C   --      by column and row (column 1 is not necessarily the center)

      INTEGER IELROW(*)
      INTEGER IXROW(NEROW+1)
C...Assert NUMCOL > 0
      INTEGER IELCEN(NUMCOL,*)

C   --Put rows into row by column array

      DO 20 IROW = 1, NEROW
         IX = IXROW(IROW)
         DO 10 ICOL = 1, NUMCOL
            IF (IX .LE. IXROW(IROW+1)-1) THEN
               IELCEN(ICOL,IROW) = IELROW(IX)
               IX = IX + 1
            ELSE
               IELCEN(ICOL,IROW) = 0
            END IF
   10    CONTINUE
   20 CONTINUE

      RETURN
      END
