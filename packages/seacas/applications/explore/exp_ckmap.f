C    Copyright(C) 2008-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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
C
C=======================================================================
      SUBROUTINE CKMAP (ICNT, MAP, INDX, TYPE)
C=======================================================================

C   --*** CKMAP *** (EXPLORE) Check database node/element map
C   --
C   --CKMAP checks the node/element order map.
C   --
C   --Parameters:
C   --   ICNT - IN - the number of nodes/elements
C   --   MAP  - IN - the node/element order map
C   --   INDX - SCRATCH - size = ICNT

      include 'exp_errcnt.blk'
      INTEGER MAP(*)
      INTEGER INDX(*)
      CHARACTER*(*) TYPE

      CHARACTER*1024 STRA
C   --Check that each node/element appears once and only once in the map
C     The 'map(i)' values may be larger than icnt, so we can't do a
C     simple check.  Instead, we do an indexed sort and check for no
C     duplicate adjacent values.

      nerr = 0
      if (icnt .le. 1) return

      CALL INDEXX (MAP, INDX, ICNT, .TRUE.)

C ... There has been a request to show min and max ids to help with
C     debugging potential database corruption issues.  Do it here.

      write (stra, 10001) type, map(indx(1)), map(indx(icnt))
10001 FORMAT('INFO: ', A, ' global id range: ',I12, ' to ', I12)
      call sqzstr(stra, lstra)
      CALL PRTERR ('CMDSPEC', STRA(:lstra))

      ILAST = MAP(INDX(1))
      DO 100 IEL = 2, ICNT
        if (map(indx(iel)) .eq. ilast) then
           if (nerr .lt. maxerrs .or. maxerrs .le. 0) then
              write (stra, 10000) type, ilast, type,
     *             indx(iel-1), indx(iel)
10000         FORMAT('MAP ERROR: ',A,' global  id ',I12,
     *             ' assigned to ',A,'s', I12,' and ',I12,'.')
              call sqzstr(stra, lstra)
              CALL PRTERR ('CMDSPEC', STRA(:lstra))
           else if (nerr .eq. maxerrs .and. maxerrs .gt. 0) then
              call prterr('CMDSPEC',
     $             '...skipping additional errors...')
           end if
           nerr = nerr + 1
        end if
        ilast = map(indx(iel))
  100 CONTINUE
      if (nerr .gt. 0) then
         write (stra, 10010) nerr, type
10010    FORMAT('MAP ERROR: Found ',I12,' errors in ',A,' map check.')
         call sqzstr(stra, lstra)
         CALL PRTERR ('CMDSPEC', STRA(:lstra))
      end if
      RETURN
      END
