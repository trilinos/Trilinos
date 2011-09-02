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
      SUBROUTINE RDNPS (NDB, NUMNPS, LNPSNL, IDNPS, NNNPS, NDNPS,
     $     IXNNPS, IXDNPS, LTNNPS, FACNPS, NAME, ISEOF, NAMLEN)
C=======================================================================

C   --*** RDNPS *** (GROPE) Read database nodal point sets
C   --
C   --RDNPS reads the nodal point set information from the database.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   NUMNPS - IN - the number of nodal points sets
C   --   LNPSNL - IN - the length of the nodal point sets node list
C   --   IDNPS - OUT - the nodal point set ID for each set
C   --   NNNPS - OUT - the number of nodes for each set
C   --   NDNPS - OUT - the number of distribution factors for each set
C   --   IXNNPS - OUT - the index of the first node for each set
C   --   IXDNPS - OUT - the index of the first dist factor for each set
C   --   LTNNPS - OUT - the nodes for all sets
C   --   FACNPS - OUT - the distribution factors for all sets
C   --   ISEOF - IN/OUT - set true if end of file read
C   --
C   --Database must be positioned at start of nodal point set information
C   --upon entry; upon exit at end of nodal point set information.

      include 'params.blk'
      INTEGER IDNPS(*)
      integer NNNPS(*)
      INTEGER NDNPS(*)
      INTEGER IXNNPS(*)
      INTEGER IXDNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)
      CHARACTER*(NAMLEN) NAME(*)
      LOGICAL ISEOF

      CHARACTER*80 ERRMSG
      CHARACTER*32 STRA

      CALL INIINT (NUMNPS, 0, IDNPS)
      CALL INIINT (NUMNPS, 0, NNNPS)
      CALL INIINT (NUMNPS, 0, NDNPS)
      CALL INIINT (NUMNPS, 0, IXNNPS)
      CALL INIINT (NUMNPS, 0, IXDNPS)
      CALL INIINT (LNPSNL, 0, LTNNPS)
      CALL INIREA (LNPSNL, 0.0, FACNPS)

      if (numnps .le. 0) return

C ... Read nodeset ids for all sets
      call exgnsi (ndb, idnps, ierr)
      if (ierr .ne. 0) go to 100

c ... Check that all ids are unique
      do 80 i = 1, numnps
         if (locint (idnps(i), i-1, idnps) .gt. 0) then
            call intstr (1, 0, idnps(i), stra, lstra)
            call prterr ('CMDERR',
     &         'nodeset id ' // stra(:lstra) // ' is not unique')
         end if
   80 continue

C ... Read nodeset parameters
      iens = 1
      ieds = 1
      do 90 i = 1, numnps
        call exgnp  (ndb, idnps(i), nnnps(i), ndnps(i), ierr)
        if (ierr .ne. 0) go to 110
        if (nnnps(i) .ne. ndnps(i) .and. ndnps(i) .ne. 0) then
           WRITE (ERRMSG, 10000)
     &          'Number of nodes does not match number of dist factors',
     &          i
           CALL SQZSTR(ERRMSG, LSTR)
           CALL PRTERR('WARNING', ERRMSG(:LSTR))
        end if
        
        ixnnps(i) = iens
        ixdnps(i) = ieds

        if (nnnps(i) .gt. 0) then
          call exgns  (ndb, idnps(i), ltnnps(iens), ierr)
          if (ierr .ne. 0) go to 130
          if (ndnps(i) .gt. 0) then
            call exgnsd (ndb, idnps(i), facnps(ieds), ierr)
            if (ierr .ne. 0) go to 140
          end if
        end if
        
        iens = iens + nnnps(i)
        ieds = ieds + ndnps(i)
 90   continue

C ... Read names (if they exist)
      CALL EXGNAMS(NDB, EXNSET, numnps, name, ierr) 
      RETURN

  100 CONTINUE
      WRITE (ERRMSG, 10000) 'IDS of nodal point sets'
      GOTO 170
  110 CONTINUE
      WRITE (ERRMSG, 10000) 'nodal point set PARAMETERS'
      GOTO 170
  120 CONTINUE
      WRITE (ERRMSG, 10000) 'POINTERS to nodal point sets'
      GOTO 170
  130 CONTINUE
      WRITE (ERRMSG, 10000) 'Nodal point sets NODES'
      GOTO 170
 140  CONTINUE
      WRITE (ERRMSG, 10000) 'Nodal point sets FACTORS'
      GOTO 170
 170  CONTINUE
      CALL WDBERR (IERR, ERRMSG)
      ISEOF = .TRUE.
 180  CONTINUE


10000  FORMAT (A,' in nodeset ', I10)
      END
