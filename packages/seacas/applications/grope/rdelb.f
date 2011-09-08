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
      SUBROUTINE RDELB (NDB, NELBLK, IDELB, NUMELB, NUMLNK, NUMATR,
     &   A, C, KLINK, KATRIB, KATRNM, ISEOF, EBTYPE, EBNAME, NAMLEN)
C=======================================================================

C   --*** RDELB *** (GROPE) Read database element blocks
C   --
C   --RDELB reads the element block information from the database.
C   --Some dynamic dimensioning is done.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NELBLK - IN - the number of element blocks to read
C   --   IDELB - OUT - the element block ID for each block
C   --   NUMELB - OUT - the number of elements for each block
C   --   NUMLNK - OUT - the number of nodes per element for each block
C   --   NUMATR - OUT - the number of attributes for each block
C   --   A - IN - the dynamic memory base array
C   --   KLINK - OUT - the dynamic memory pointer to the connectivity array
C   --      (named 'LINK')
C   --   KATRIB - OUT - the dynamic memory pointer to the attribute array
C   --      (named 'ATRIB')
C   --   ISEOF - IN/OUT - set true if end of file read

      include 'exodusII.inc'

      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      CHARACTER*(MXSTLN) EBTYPE(*)
      CHARACTER*(NAMLEN) EBNAME(*)
      
      DIMENSION A(*)
      CHARACTER*1 C(*)
      LOGICAL ISEOF

      CHARACTER*80 ERRMSG

      CALL INIINT (NELBLK, 0, IDELB)
      CALL INIINT (NELBLK, 0, NUMELB)
      CALL INIINT (NELBLK, 0, NUMLNK)
      CALL INIINT (NELBLK, 0, NUMATR)
      CALL INISTR (NELBLK, ' ', EBTYPE)
      CALL INISTR (NELBLK, ' ', EBNAME)

C ... Get element block ids
      if (nelblk .gt. 0) then
        call exgebi(ndb, idelb, ierr)
        if (ierr .ne. 0) go to 120
      end if

C ... Read element block sizing parameters
      IELNK = 0
      IEATR = 0
      INATR = 0
      DO 100 IELB = 1, NELBLK
        call exgelb(ndb, idelb(ielb), ebtype(ielb), numelb(ielb),
     &       numlnk(ielb), numatr(ielb), ierr)
        if (ierr .ne. 0) go to 120
        
        if (ebtype(ielb) .eq. 'nsided' .or.
     *      ebtype(ielb) .eq. 'NSIDED') THEN
          IELNK = IELNK + NUMLNK(IELB) 
        else
          IELNK = IELNK + NUMLNK(IELB) * NUMELB(IELB) 
      end if
      IEATR = IEATR + NUMATR(IELB) * NUMELB(IELB)
      INATR = INATR + NUMATR(IELB)
 100  CONTINUE

      CALL MDRSRV ('LINK',  KLINK,  IELNK)
      CALL MDRSRV ('ATRIB', KATRIB, IEATR)
      call mcrsrv ('ATRNM', KATRNM, INATR*NAMLEN)
      CALL MDSTAT (NERR, MEM)
      if (nerr .gt. 0) go to 140

C ... Read element block connectivity and attributes
      ielnk = 0
      ieatr = 0
      inatr = 0
      do 110 ielb = 1, nelblk
        islnk = ielnk + 1
        if (ebtype(ielb) .eq. 'nsided' .or.
     *    ebtype(ielb) .eq. 'NSIDED') THEN
          ielnk = islnk + numlnk(ielb) - 1
        else
          ielnk = islnk + numlnk(ielb) * numelb(ielb) - 1
        end if
        isatr = ieatr + 1
        ieatr = isatr + numatr(ielb) * numelb(ielb) - 1

        CALL RDEB1 (NDB,
     &    IDELB(IELB), NUMELB(IELB), NUMLNK(IELB), NUMATR(IELB),
     &    A(KLINK+ISLNK-1), A(KATRIB+ISATR-1), 
     &    C(KATRNM+NAMLEN*INATR), NAMLEN)

        inatr = inatr + numatr(ielb)
 110  CONTINUE

C ... Read element block names (if they exist)
      CALL EXGNAMS(NDB, EXEBLK, nelblk, ebname, ierr) 
      RETURN

  120 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     &   'ELEMENT BLOCK SIZING PARAMETERS for block', IELB
      GOTO 130
  130 CONTINUE
      CALL WDBERR (IERR, ERRMSG)
      ISEOF = .TRUE.
  140 CONTINUE
      RETURN

10000  FORMAT (5 (A, I10))
      END
