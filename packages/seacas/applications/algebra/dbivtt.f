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
      SUBROUTINE DBIVTT (NDB, ISEVOK, ITMP, NELBLK, NVAREL)
C=======================================================================
C   --*** DBIVTT *** Read element variable truth table
C   --   Modified for ExodusII format 8/26/95
C   --*** DBINAM *** (EXOLIB) Read database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBINAM performed a number of different input file read base
C   --on the passed in option argument.  DBINAM was split up
C   --into a number of different subroutins
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   ISEVOK - OUT - the dynamic memory array of the element block variable
C   --                  truth table;variable i,block j exists iff ISEVOK(j,i)
C   --   ITMP   - OUT - temporary scratch array
C   --   NELBLK - IN  - the number of element blocks
C   --   NVAREL - IN  - the number of element variables;

      include 'params.blk'

      INTEGER NDB
      INTEGER NELBLK, NVAREL
      LOGICAL ISEVOK(NELBLK, NVAREL)
      INTEGER ITMP(NVAREL, NELBLK)
      LOGICAL LVALUE

C     Read the element block variable truth table
C     call exgvtt(fileid, num_elem_blks, num_elem_var,
C                 isevok(num_elem_var, num_elem_blks, errorid)
C       isevok - num_elem_var cycles faster
      if (nvarel .gt. 0) then
         CALL EXGVTT(NDB, NELBLK, NVAREL, ITMP, IERR)
         DO 110 I = 1, NVAREL
            DO 100 IELB = 1, NELBLK
               IF (ITMP(I,IELB) .EQ. 0) LVALUE = .FALSE.
               IF (ITMP(I,IELB) .EQ. 1) LVALUE = .TRUE.
               IF ((ITMP(I,IELB) .NE. 0) .AND.
     &             (ITMP(I,IELB) .NE. 1)) THEN
                WRITE(*,*)'Error reading element variable truth table'
               END IF
               ISEVOK(IELB,I) = LVALUE
 100        CONTINUE
 110     CONTINUE
      end if

      RETURN
      END

