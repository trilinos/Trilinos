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
      SUBROUTINE DBINM1 (NDB, TYPE, OPTION, NBLK, NVAR, ISEVOK, IEVOK,
     &                   ITMP, IERR, NBDM, *)
C=======================================================================
C   --*** DBINM1 *** (EXOLIB) Internal to DBINAM
C   --   Written by Amy Gilkey - revised 02/18/88
C   --
C   --DBINM1 reads the element block variable truth table.
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   OPTION - IN  - ' ' to not store, '*' to store all, else store options:
C   --                  'T' to store element block variable truth table
C   --   NELBLK - IN  - the number of element blocks
C   --   NVAREL - IN  - the number of element variables
C   --   ISEVOK - OUT - the element block variable truth table;
C   --                  variable i of block j exists iff ISEVOK(j,i)
C   --   IEVOK  - OUT - the element block variable truth table;
C   --                  variable i of block j exists iff ISEVOK(j,i) is NOT 0
C   --   IERR   - OUT - the returned read error flag
C   --   *      - OUT - return statement if error encountered
C   --                  NO message is printed
C   --
C   --Database must be positioned in front of truth table upon entry;
C   --upon exit positioned after table.

      include 'exodusII.inc'
      INTEGER NDB
      CHARACTER*1 TYPE
      CHARACTER*(*) OPTION
      INTEGER NBLK, NVAR
      LOGICAL ISEVOK(NBDM,*)
      INTEGER IEVOK(NBDM,*)
      INTEGER ITMP(NVAR,NBDM)
      INTEGER IERR

      if (nvar .le. 0) return
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
         if (type .eq. 'E') then
            call exgvtt (ndb, nblk, nvar, itmp, ierr)
         else if (type .eq. 'M') then
            call exgnstt (ndb, nblk, nvar, itmp, ierr)
         else if (type .eq. 'S') then
            call exgsstt (ndb, nblk, nvar, itmp, ierr)
         end if
         IF (ierr .eq. 17) then
            DO 20 I = 1, NVAR
               DO 10 IELB = 1, NBLK
                  ISEVOK(IELB,I) = .true.
 10            CONTINUE
 20         CONTINUE
         ELSE
            DO 110 I = 1, NVAR
               DO 100 IELB = 1, NBLK
                  ISEVOK(IELB,I) = (ITMP(I,IELB) .NE. 0)
 100           CONTINUE
 110        CONTINUE
         ENDIF

      END IF
      
      RETURN
      
      END
