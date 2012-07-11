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
      SUBROUTINE RNDVAR (A, IVAR, IELBLK, INSTEP, LENVAR, VAR)
C=======================================================================

C   --*** RNDVAR *** (BLOT) Read variable
C   --   Written by Amy Gilkey - revised 07/27/88
C   --
C   --RNDVAR returns the values for the requested variable for the
C   --requested time step.  It either reads the values from the sequential
C   --database file and writes them to a direct access scratch file or it
C   --reads the values from the direct access scratch file.
C   --
C   --This routine uses MDFIND to find the following dynamic memory arrays:
C   --   NUMELB - the number of elements per element block
C   --   ISEVOK - the element block variable truth table;
C   --            variable i of block j exists iff ISEVOK(j,i)
C   --   WHOTIM - true iff whole (versus history) time step
C   --
C   --Parameters:
C   --   A      - IN  - the dynamic memory base array
C   --   IVAR   - IN  - the variable index
C   --                  (0 to initialize random file only)
C   --   IELBLK - IN  - the element block number, <=0 for all
C   --                  (for element blocks only)
C   --   INSTEP - IN  - the time step number
C   --                  = +n to read time step n
C   --                  = -n to transfer time step n to random file only
C   --                  =  0 to transfer all time steps to random file
C   --   LENVAR - IN  - the length of VAR
C   --   VAR    - OUT - the variable values (indeterminate if INSTEP <= 0)
C   --
C   --Common Variables:
C   --   Uses NDB of /DBASE/
C   --   Uses NUMNP, NUMEL, NELBLK, NVARHI, NVARGL, NVARNP, NVAREL,
C   --      NSTEPS, NSTEPW of /DBNUMS/
C   --
C   --Database is rewound upon the first entry of this routine; upon
C   --exit a flag is set to keep track of the database position; the
C   --database should not be moved between calls to this routine.
C   --
C   --A scratch random file is created and read and written in this routine.
C   --It is connected to unit 90.

      include 'exodusII.inc'
      include 'dbase.blk'
      include 'dbnums.blk'

      REAL VAR(*)

      DIMENSION A(*)

      CHARACTER TYP

      INTEGER KNELB, KIEVOK, KIDELB
C      --KNELB - the dynamic memory index of NUMELB - the number of elements
C      --   per element block
C      --KIEVOK - the dynamic memory index of ISEVOK - the element block
C      --   variable truth table; variable i of block j exists iff ISEVOK(j,i)

      IF (IVAR   .EQ. 0) return
      if (instep .le. 0) return

      CALL DBVTYP_BL (IVAR, TYP, IDVAR)

      if (typ .eq. 'G') then
        if (lenvar .lt. nvargl) then
          call prterr ('PROGRAM', 'Invalid Array length in rndvar')
        end if
        call exggv(ndb, instep, nvargl, var, ierr)
      end if

      if (typ .eq. 'N') then
        if (lenvar .lt. numnp) then
          call prterr ('PROGRAM', 'Invalid Array length in rndvar')
        end if
        call exgnv(ndb, instep, idvar, numnp, var, ierr)
      END IF

      if (typ .eq. 'E') then
        if (lenvar .lt. numel) then
          call prterr ('PROGRAM', 'Invalid Array length in rndvar')
        end if
        call inirea(numel, 0.0, var)
        if (ielblk .le. 0) then
          imin = 1
          imax = nelblk
        else
          imin = ielblk
          imax = ielblk
        end if
        CALL MDFIND ('ISEVOK', KIEVOK, IDUM)
        CALL MDFIND ('IDELB',  KIDELB, IDUM)
        CALL MDFIND ('NUMELB', KNELB, IDUM)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 130
        
        call rndelv(ndb, imin, imax, instep, idvar, a(kidelb),
     &              a(knelb), var(1), a(kievok), nelblk)
      END IF
      
 130  CONTINUE
      RETURN
      END

      subroutine rndelv(ndb, imin,imax, instep, idvar,
     &                  idelb, numelb, var, isevok, nelblk)

      real var(*)
      integer idelb(*), numelb(*)
      logical isevok(nelblk, *)

      ibeg = 1
      do 10 iel = imin, imax
C--- Check truth table.
        if (isevok(iel, idvar)) then
           call exgev(ndb, instep, idvar, idelb(iel), numelb(iel),
     &                var(ibeg), ierr)
        end if
        ibeg = ibeg + numelb(iel)
 10   continue
      
      return
      end
C --Converted
