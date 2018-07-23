C Copyright (c) 2007-2017 National Technology & Engineering Solutions of
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
C 
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
C 

C $Id: rwname.f,v 1.4 2007/10/17 18:47:22 gdsjaar Exp $
C=======================================================================
      SUBROUTINE RWNAME (NTXT, NDB, NELBLK, NVARGL, NVARNP, NVAREL,
     &  A, C, KIEVOK, EXODUS, NAMLEN, *)
C=======================================================================

C   --*** RDNAME *** (TXTEXO) Read database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --RDNAME reads the names of the coordinates, the element block types,
C   --and the database variables.  The element block variable truth table
C   --is also read.
C   --
C   --Note that the numbers of variables are read in this routine.
C   --
C   --This routine calls DBVINI and uses DBVIX to get the variable name
C   --indices.
C   --
C   --Dynamic memory is reserved in this routine.  If there is a problem,
C   --the routine returns normally without printing an error message.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NDIM - IN - the number of coordinates per node
C   --   NELBLK - IN - the number of element blocks
C   --   NVARGL - OUT - the number of global variables
C   --   NVARNP - OUT - the number of nodal variables
C   --   NVAREL - OUT - the number of element variables
C   --   A - In - the dynamic memory base array
C   --   C - IN - the dynamic character memory base array
C   --   KIEVOK - OUT - the dynamic memory index of the element block variable
C   --      truth table; variable i of block j exists iff ISEVOK(j,i)
C   --   EXODUS - OUT - set false if GENESIS format, true if EXODUS so far
C   --   * - return statement if error encountered, including end-of-file;
C   --      NOT used if valid GENESIS file; message is printed
C   --
C   --Database must be positioned in front of coordinate names upon entry;
C   --upon exit positioned after names.

      include 'exodusII.inc'
      character*1 c(*)
      DIMENSION A(*)
      logical exodus

      NVARGL = -999
      NVARNP = -999
      NVAREL = -999

C   --Read variable names
      READ (NTXT, *, END=130, ERR=130)
      READ (NTXT, *, END=130, ERR=130) NVARGL, NVARNP, NVAREL

C ... Allocate array for names
      MAXVAR = max(nvargl, nvarnp, nvarel)
      if (maxvar .gt. 0) then
        exodus = .true.
      else
        exodus = .false.
        go to 100
      end if

      call mcrsrv('NAMES', kname, namlen*maxvar)
      call mcstat(nerr, mem)
      if (nerr .gt. 0) then
        return 1
      end if
      
      call rwnam1(ntxt, ndb, 'G', nvargl, c(kname), namlen, *180)
      call rwnam1(ntxt, ndb, 'N', nvarnp, c(kname), namlen, *180)
      call rwnam1(ntxt, ndb, 'E', nvarel, c(kname), namlen, *180)
      
      call mcdel('NAMES')
      call mcstat(nerr, mem)
      if (nerr .gt. 0) then
        return 1
      end if


C   --Read the element block variable truth table
      CALL MDRSRV ('ISEVOK', KIEVOK, NELBLK * NVAREL)
C ... Temporary logical array
      CALL MDRSRV ('LSEVOK', KLEVOK, NELBLK * NVAREL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

      if (nvarel .gt. 0) then
        CALL RDNM1 (NTXT, NDB, NELBLK, NVAREL,
     &    A(KIEVOK), A(KLEVOK), *180)
      end if

      CALL MDDEL ('LSEVOK')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100
  100 CONTINUE
      RETURN

  130 CONTINUE
      CALL PRTERR ('FATAL', 'Reading NUMBER OF VARIABLES')
      GOTO 180
  180 CONTINUE
      RETURN 1
      END

      subroutine rwnam1(ntxt, ndb, flag, nvar, names, namlen, *)
      include 'exodusII.inc'
      character*1 flag
      character*(namlen) names(*)
      
      if (nvar .eq. 0) then
        read (ntxt,*,end=150, err=150)
        return
      else
        READ (NTXT, '(2 (A, 1X))', END=150, ERR=150)
     &    (NAMES(I), I=1, nvar)
        
        call expvp(ndb, flag, nvar, ierr)
        call expvan(ndb, flag, nvar, names, ierr)
        if (ierr .ne. 0) go to 160
      end if
      return
      
  150 CONTINUE
      CALL PRTERR ('FATAL', 'Reading VARIABLE NAMES')
      return 1
  160 CONTINUE
      CALL PRTERR ('FATAL', 'Writing VARIABLE NAMES')
      return 1
      end
      
