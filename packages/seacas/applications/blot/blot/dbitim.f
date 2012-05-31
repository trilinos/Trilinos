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
      SUBROUTINE DBITIM (NDB, OPTION, EXODUS,
     &   NVARNP, NELBLK, NVAREL, ISEVOK,
     &   NSTEPS, NSTEPW, A, KTIMES, KWHOLE, *)
C=======================================================================
C$Log: dbitim.f,v $
CRevision 1.3  2009/03/25 12:36:43  gdsjaar
CAdd copyright and license notice to all files.
CPermission to assert copyright has been granted; blot is now open source, BSD
C
CRevision 1.2  1997/11/11 14:55:54  gdsjaar
CAdded 'external blkdat' to main program to ensure that the block data
Cgets linked into the executable. Wasn't happening on dec alpha
Csystems.
C
CRemoved unreachable lines in several routines
C
CFixed variable name spelling in contor.f
C
CUnsplit strings that were split across lines
C
CRemoved old error variables left over from exodusIIv1
C
CUpped version number
C
CRevision 1.1  1994/04/07 19:58:43  gdsjaar
CInitial checkin of ACCESS/graphics/blotII2
C
C Revision 1.4  1993/10/13  15:03:54  gdsjaar
C Change so does not die if there are variables on database, but no timesteps.
C
C Revision 1.3  1993/07/27  22:39:05  gdsjaar
C Cleaned up some unused variables, fixed some parameter lists.
C
C Revision 1.2  1993/07/27  19:31:38  gdsjaar
C Fixed setting of EXODUS logical
C
c Revision 1.1  1993/07/27  13:30:41  gdsjaar
c Initial checkin of ACCESS/graphics/blotII
c
c Revision 1.2  1990/10/11  16:46:39  gdsjaar
c Fixed problem with singly/doubly dimensioned array
c
c Revision 1.1.1.1  90/08/14  16:13:07  gdsjaar
c Testing
c 
c Revision 1.1  90/08/14  16:13:06  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:13  gdsjaar
c Initial revision
c 

C   --*** DBITIM *** (EXOLIB) Read database time step times
C   --   Written by Amy Gilkey - revised 11/11/87
C   --
C   --DBITIM reads all time steps of the database, storing the time for
C   --each step in dynamic memory.  It also sets the number of time steps.
C   --If there is a dynamic memory problem, the routine returns immediately
C   --without issuing an error message.  If there is a read error, a
C   --warning is printed and the step with the error is ignored.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'W' to store whole times only
C   --      'H' to store history times only
C   --   EXODUS - IN - true iff this is an EXODUS file
C   --   NVARNP - IN - the number of nodal variables
C   --   NELBLK - IN - the number of element blocks
C   --   NVAREL - IN - the number of element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   NSTEPS - OUT - the number of database time steps
C   --   NSTEPW - OUT - the number of whole (versus history) database time steps
C   --   A - IN/OUT - the dynamic memory base array
C   --   KTIMES - OUT - the dynamic memory index of TIMES (if OPTION)
C   --   KWHOLE - OUT - the dynamic memory index of WHOTIM (if OPTION);
C   --      WHOTIM(i) is true iff TIMES(i) is a whole (versus history) time step
C   --   * - return statement if error encountered, message is printed
C   --
C   --Database must be positioned in front of first time step upon entry;
C   --upon exit positioned at end of file.

      include 'exodusII.inc'
      INTEGER NDB
      CHARACTER*(*) OPTION
      LOGICAL EXODUS
      INTEGER NVARNP, NELBLK, NVAREL
C     --NOTE: ISEVOK is doubly-dimensioned array (NELBLK,NVAREL) 
      LOGICAL ISEVOK(*)
      INTEGER NSTEPS, NSTEPW
      REAL A(*)
      INTEGER KTIMES
      INTEGER KWHOLE
      CHARACTER*8 CDUM

      LOGICAL STOWHO, STOHIS

      STOWHO = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'W') .GT. 0))
      STOHIS = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'H') .GT. 0))

      IF (.NOT. EXODUS) GOTO 140

      call exinq(ndb, EXTIMS, nsteps, rdum, cdum, ierr)
      nstepw = nsteps
C      --Reserve memory for step times
      if (nsteps .le. 0) then
        EXODUS = .FALSE.
        goto 140
      end if
      CALL MDRSRV ('TIMES', KTIMES, nsteps)
      IF (STOWHO .AND. STOHIS) THEN
        CALL MDRSRV ('WHOTIM', KWHOLE, nsteps)
      END IF
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 150
      
      WRITE (*, 10000) 'Please wait while the step times are read in'
      
C   --Set NRECST = the number of database records in a whole time step
      
 120  CONTINUE
      IF (.TRUE.) THEN
        
C      --Read step time
        
        call exgatm(ndb, a(ktimes), ierr)
C         --Store step time after time step read correctly
        CALL INILOG (nsteps, .TRUE., A(KWHOLE))
      END IF
      
  140 CONTINUE

  150 CONTINUE
      RETURN

10000  FORMAT (/, 1X, A)
      END
