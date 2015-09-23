C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
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
      SUBROUTINE DBINM1 (NDB, OPTION, NELBLK, NVAREL, ISEVOK, IEVOK,
     &                   ITMP, IERR, NELBDM, IDELB, ISHEX, HEXID, A,
     &                   IA, *)
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
      DIMENSION A(*), IA(*)
      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NELBLK, NVAREL
      LOGICAL ISEVOK(NELBDM,*)
      INTEGER IEVOK(NELBDM,*)
      INTEGER ITMP(NVAREL,NELBDM)
      INTEGER IERR
      INTEGER IDELB(*)
      INTEGER ISHEX
      INTEGER HEXID(*)

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
         CALL EXGVTT(NDB, NELBLK-ISHEX, NVAREL, ITMP, IERR)
C        If ishex>0 then new element blocks have been created for
C        HEXSHELL elements
         IF (ISHEX .GT. 0) THEN
C           Reserve scratch space to create new element variable truth table
            CALL MDRSRV ('IETSCR', IETSCR, NELBLK*NVAREL)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) RETURN 1
            CALL MODEVTT(NELBLK, NVAREL, ISHEX, IDELB, HEXID, ITMP,
     &                   IA(IETSCR))
            CALL MDDEL ('IETSCR')
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) RETURN 1
         ENDIF
         IF (ierr .eq. 17) then
            DO 20 I = 1, NVAREL
               DO 10 IELB = 1, NELBLK
                  ISEVOK(IELB,I) = .true.
 10            CONTINUE
 20         CONTINUE
         ELSE
            DO 110 I = 1, NVAREL
               DO 100 IELB = 1, NELBLK
                  ISEVOK(IELB,I) = (ITMP(I,IELB) .NE. 0)
 100           CONTINUE
 110        CONTINUE
         ENDIF

      END IF
      
      RETURN
      
      END

      SUBROUTINE MODEVTT(NELBLK, NVAREL, ISHEX, IDELB, HEXID, IEVTT,
     &                   TTSCR)
C     Add rows to element variable truth table for the SHELL 
C     element blocks of the HEXSHELL.  IEVTT contains the element variable
C     truth table from the exodusII file.  HEXID contains the ID of
C     the HEXSHELL element blocks. TTSCR is a scratch array for the
C     element variable truth table.
C     Design - loop over the number of element blocks. If the element 
C     block id in the IDELB matches the element block id in HEXID, copy the
C     row from the IEVTT array to the TTSCR array and then add a
C     row of zeros into the TTSCR array for the next element block
C     to account for the SHELL element block.  Otherwise, copy the
C     row from the IEVTT array to the TTSCR array.

C ... Current truth table is filled in (NVAREL,NELBLK-ISHEX)
C     and contains valid values only for the non-hexshell-shell blocks
C     (if a block was originally a shell, it has variables, if shell
C     block created from hexshell, then no variables)
C ... The IDELB has both hexshell-shell and hex ids already in it
C     so it's rows don't necessarily correspond to the truth table
C     rows.
      
      INTEGER NELBLK, NVAREL, ISHEX
      INTEGER IDELB(*), HEXID(*) 
      INTEGER IEVTT(NVAREL,NELBLK), TTSCR(NVAREL,NELBLK)

C     Index into HEXID array
      IHEX  = 1
C     Index into truth table array (IEVTT), max value = NELBLK-ISHEX
      IFROM = 1
C     Index into block id array (IDELB)      
      IID   = 1
C     Index into TTSCR array
      ITO   = 1

 200  CONTINUE
C     Exit condition
      IF (IFROM .GT. NELBLK-ISHEX) GOTO 300

C     HEXSHELL element block ID - HEX element block ID
      IF (IHEX .GT. ISHEX) THEN
         IDH = 0
      ELSE
         IDH = HEXID(IHEX)
      ENDIF

      IF (IDELB(IID) .EQ. IDH) THEN
C        Element block id matches hexshell element
         DO 210 J = 1, NVAREL
C           Row for HEX element in element variable truth table
            TTSCR(J,ITO)   = IEVTT(J,IFROM)
C           Row for SHELL element in element variable truth table
            TTSCR(J,ITO+1) = 0
 210     CONTINUE
         ITO  = ITO + 2
         IHEX = IHEX + 1
C ... Skip past hexshell-shell id in idelb
         IID  = IID + 2
      ELSE
C        Element block id doesn't match hexshell element. Copy Element
C        variable truth table values.
         DO 220 J = 1, NVAREL
            TTSCR(J,ITO)   = IEVTT(J,IFROM)
 220     CONTINUE
         ITO = ITO + 1
         IID = IID + 1
      ENDIF
      IFROM = IFROM + 1
      GOTO 200
 300  CONTINUE

C     Copy TTSCR back to IEVTT
      DO 230 I = 1, NELBLK
         DO 240 J = 1, NVAREL
            IEVTT(J,I) = TTSCR(J,I)
 240     CONTINUE
 230  CONTINUE

      RETURN
      END
