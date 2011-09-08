C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C 
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
C 

C $Id: rdnm1.f,v 1.5 2007/10/17 18:47:22 gdsjaar Exp $
C=======================================================================
      SUBROUTINE RDNM1 (NTXT, NDB, NELBLK, NVAREL, ISEVOK, LSEVOK, *)
C=======================================================================

C   --*** RDNM1 *** (TXTEXO) Internal to RDNAME
C   --   Written by Amy Gilkey - revised 02/22/88
C   --
C   --RDNM1 reads the element block variable truth table.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NELBLK - IN - the number of element blocks
C   --   NVAREL - IN - the number of element variables
C   --   ISEVOK - OUT - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   * - return statement if error encountered, including end-of-file;
C   --      NO message is printed
C   --
C   --Database must be positioned in front of truth table upon entry;
C   --upon exit positioned after table.

      INTEGER ISEVOK(NVAREL,*)
      LOGICAL LSEVOK(NVAREL,*)

      CHARACTER*5 STRA

C ... Nothing to read if NVAREL == 0
      if (nvarel .eq. 0) return
      
      IELB = 0
      READ (NTXT, *, END=110, ERR=110)
      DO 100 IELB = 1, NELBLK
        READ (NTXT, *, END=110, ERR=110) (LSEVOK(I,IELB), I=1,NVAREL)
        DO 90 I = 1, NVAREL
          if (lsevok(i,ielb)) then
            isevok(i,ielb) = 1
          else
            isevok(i,ielb) = 0
          end if
 90     continue
 100  CONTINUE
      
      call expvtt(ndb, nelblk, nvarel, isevok, ierr)
      
      RETURN

 110  CONTINUE
      CALL INTSTR (1, 0, IELB, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &  'Reading ELEMENT BLOCK VARIABLE TRUTH TABLE for block '
     &  // STRA(:LSTRA))
      GOTO 120
 120  CONTINUE
      RETURN 1
      END
