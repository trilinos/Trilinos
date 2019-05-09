C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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
      SUBROUTINE OPNLOG (LOGU)
C=======================================================================

C   --*** OPNLOG *** (BLOT) Open log file and write header
C   --   Written by Amy Gilkey - revised 01/11/88
C   --
C   --OPNLOG opens the log file and writes the command line as the header
C   --for the log file.
C   --
C   --Parameters:
C   --   NLOG - IN/OUT - the log file number; returned <= if log file
C   --      cannot be opened
C   --
C   --Common Variables:
C   --   Uses NDB of /DBASE/
C   --   Uses QAINFO of /PROGQA/

      include 'params.blk'
      include 'progqa.blk'
      include 'dbase.blk'
      include 'dbname.blk'

      CHARACTER*2048 INLINE, FILNAM, ERRMSG
      CHARACTER*256 STR
      LOGICAL ISON

      NLOG = LOGU
      filnam = basenam(:lenstr(basenam)) // '.blot.log'

      open (unit=nlog, file=filnam(:lenstr(filnam)), form='formatted',
     *  status='unknown', iostat=ierr)
      IF (IERR .NE. 0) THEN
        ERRMSG = 'Log file "'//FILNAM(:LENSTR(FILNAM))//
     *    '" could not be opened.'
        CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
        GOTO 100
      END IF

      INLINE = '$$$ ' // QAINFO(1)
      L = LENSTR (INLINE) + 1

      IF (L .LT. LEN (INLINE)) INLINE(L+1:) = DBNAME
      L = LENSTR (INLINE) + 1

      CALL GRGPAR ('DEVICE', 1, ISON, STR)
      IF (ISON) THEN
         IF (L .LT. LEN (INLINE)) INLINE(L+1:) = STR
      ELSE
         IF (L .LT. LEN (INLINE)) INLINE(L+1:) = '""'
      END IF
      L = LENSTR (INLINE) + 1

      CALL GRGPAR ('DEVICE', 2, ISON, STR)
      IF (ISON) THEN
         IF (L .LT. LEN (INLINE)) INLINE(L+1:) = STR
      ELSE
         IF (L .LT. LEN (INLINE)) INLINE(L+1:) = '""'
      END IF
      L = LENSTR (INLINE) + 1

      WRITE (NLOG, '(A)', IOSTAT=IERR) INLINE(:L-1)
      if (ierr .ne. 0) then
         CALL PRTERR ('WARNING', 'Log file cannot be written')
         NLOG = -1
         GOTO 100
      end if
  100 CONTINUE
      LOGU = NLOG
      RETURN
      END
