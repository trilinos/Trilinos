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
      SUBROUTINE RDESS (NDB, NUMESS, LESSEL, LESSNL,
     &   IDESS, NEESS, NDESS, IXEESS, IXDESS,
     &   LTEESS, LTSESS, FACESS, NAME, ISEOF, NAMLEN)
C=======================================================================

C   --*** RDESS *** (GROPE) Read database element side sets
C   --
C   --RDESS reads the element side set information from the database.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NUMESS - IN  - the number of side sets
C   --   LESSEL - IN  - the length of the element side sets element list
C   --   LESSNL - IN  - the length of the element side sets distribution list
C   --   IDESS  - OUT - the element side set ID for each set
C   --   NEESS  - OUT - the number of elements for each set
C   --   NDESS  - OUT - the number of factors for each set
C   --   IXEESS - OUT - the index of the first element for each set
C   --   IXDESS - OUT - the index of the first factor for each set
C   --   LTEESS - OUT - the elements for all sets
C   --   LTESSS - OUT - the elements for all sets
C   --   LTSESS - OUT - the element sides for all sets
C   --   FACESS - OUT - the distribution factors for all sets
C   --   ISEOF  - IN/OUT - set true if end of file read

      include 'params.blk'
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NDESS(*)
      INTEGER IXEESS(*)
      INTEGER IXDESS(*)
      INTEGER LTEESS(*)
      INTEGER LTSESS(*)
      REAL FACESS(*)
      CHARACTER*(NAMLEN) NAME(*)
      LOGICAL ISEOF

      CHARACTER*80 ERRMSG
      CHARACTER*32 STRA

      if (numess .le. 0) return

      CALL INIINT (NUMESS, 0, IDESS)
      CALL INIINT (NUMESS, 0, NEESS)
      CALL INIINT (NUMESS, 0, NDESS)
      CALL INIINT (NUMESS, 0, IXEESS)
      CALL INIINT (NUMESS, 0, IXDESS)
      CALL INIINT (LESSEL, 0, LTEESS)
      CALL INIINT (LESSEL, 0, LTSESS)
      CALL INIREA (LESSNL, 0.0, FACESS)

C ... Read sideset ids
      call exgssi(ndb, idess, ierr)
      if (ierr .ne. 0) go to 100
      
c ... Check that all ids are unique
      do 80 i = 1, numess
         if (locint (idess(i), i-1, idess) .gt. 0) then
            call intstr (1, 0, idess(i), stra, lstra)
            call prterr ('CMDERR',
     &         'sideset id ' // stra(:lstra) // ' is not unique')
         end if
   80 continue

      
      ies = 1
      ifs = 1
      do 90 i = 1, numess
        ixeess(i) = ies
        ixdess(i) = ifs

C ... Read sideset parameters
        call exgsp  (ndb, idess(i), neess(i), ndess(i), ierr)
        kk = neess(i)
        kkk = ndess(i)
        if (ierr .ne. 0) go to 110
        
C ... Read sideset elements and faces
        if (neess(i) .gt. 0) then
          call exgss  (ndb, idess(i), lteess(ies), ltsess(ies), ierr)
          if (ierr .ne. 0) go to 150
        end if
        
C ... Read sideset distribution factors
        if (ndess(i) .gt. 0) then
           call exgssd (ndb, idess(i), facess(ifs), ierr)
           if (ierr .ne. 0) go to 170
        end if

        ies = ies + neess(i)
        ifs = ifs + ndess(i)
 90   continue      

C ... Read names (if they exist)
      CALL EXGNAMS(NDB, EXSSET, numess, name, ierr) 
      RETURN

 100  CONTINUE
      WRITE (ERRMSG, 10000) 'IDS of element side sets'
      GOTO 180
 110  CONTINUE
      WRITE (ERRMSG, 10000) 'NUMBERS OF ELEMENTS in element side sets'
      GOTO 180
 120  CONTINUE
      WRITE (ERRMSG, 10000) 'NUMBERS OF NODES in element side sets'
      GOTO 180
 130  CONTINUE
      WRITE (ERRMSG, 10000) 'ELEMENT POINTERS to element side sets'
      GOTO 180
 140  CONTINUE
      WRITE (ERRMSG, 10000) 'NODE POINTERS to element side sets'
      GOTO 180
 150  CONTINUE
      WRITE (ERRMSG, 10000) 'Element side set ELEMENTS'
      GOTO 180
 160  CONTINUE
      WRITE (ERRMSG, 10000) 'Element side set NODES'
      GOTO 180
 170  CONTINUE
      WRITE (ERRMSG, 10000) 'Element side set FACTORS'
      GOTO 180
 180  CONTINUE
      CALL WDBERR (IERR, ERRMSG)
      ISEOF = .TRUE.
 190  CONTINUE
      RETURN
      
10000 FORMAT (A)
      END
