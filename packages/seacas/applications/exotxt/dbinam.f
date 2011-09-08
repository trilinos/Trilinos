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

C     $Id: dbinam.f,v 1.7 2007/10/17 18:46:09 gdsjaar Exp $
C=======================================================================
      SUBROUTINE DBINAM (NDB, C, KNAMES, NVARGL, NVARNP, NVAREL,
     &                   IXGV, IXNV, IXEV, IOERR, NAMLEN)
C=======================================================================

C   --*** DBINAM *** Read the names of the of global, node, and
C   --               element result variables
C   --   Modified for ExodusIIV2 database 8/26/95
C   --This code originated from the subroutine
C   --*** DBINAM *** (EXOLIB) Read database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBINAM performed a number of different input file read base
C   --on the passed in option argument.  DBINAM was split up
C   --into a number of different subroutins
C   --
C   --This routine calls DBVINI and uses DBVIX to get the variable name
C   --indices.  This subroutine includes and call DBINM2.
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   C      - OUT - the dynamic memory base array (character)
C   --   NVARGL - IN  - the number of global variables; <0 if end-of-file
C   --   NVARNP - IN  - the number of nodal variables; <0 if end-of-file
C   --   NVAREL - IN  - the number of element variables; <0 if end-of-file
C   --   IXGV   - OUT - the VNAMES index of the global variable names
C   --   IXNV   - OUT - the VNAMES index of the nodal variable names
C   --   IXEV   - OUT - the VNAMES index of the element variable names
C   --   IOERR  - OUT - error flag
C   --

C   --Routines Called:
C   --   EXUPCS - (SUPES) Convert to uppercase and blank non-standard
C   --   MDRSRV - (SUPES) Reserve dynamic memory
C   --   PCKSTR - (STRLIB) Remove embedded blanks

      PARAMETER (MAXDIM=6)

      INTEGER NDB
      CHARACTER*1 C(*)
      CHARACTER*(NAMLEN) KNAMES(*)
      INTEGER NVARGL, NVARNP, NVAREL
      INTEGER IXGV, IXNV, IXEV
      INTEGER IOERR
      IOERR = 0

C     There are no global, nodal, or element variables
      IF ((nvargl + nvarnp + nvarel) .eq. 0) RETURN

C     DBVINI Initializes the indices for DBVTYP and DBVIX
C     DBVINI(num_global_var, num_nodal_var, num_elem_var)
      CALL DBVINI (NVARGL, NVARNP, NVAREL)

C     call dbvix(var_index, ret_var_type, ret_var_num)
      CALL DBVIX ('G', 1, IXGV)
      CALL DBVIX ('G', NVARGL, IXGVE)
      CALL DBVIX ('N', 1, IXNV)
      CALL DBVIX ('N', NVARNP, IXNVE)
      CALL DBVIX ('E', 1, IXEV)
      CALL DBVIX ('E', NVAREL, IXEVE)

C      --Read and pack variable names
      call dbinm2 (ndb, nvargl, nvarnp, nvarel, KNAMES(IXGV),
     &             KNAMES(IXNV), KNAMES(IXEV), IOERR, namlen)
      IF (IOERR .EQ. 1) RETURN

      RETURN
      END

C     *********************************************************
      SUBROUTINE DBINM2 (ndb, nvargl, nvarnp, nvarel,
     &  namgv, namnv, namev, ioerr, namlen)
C     *********************************************************
C     NDB    - IN  - File ID
C     NVARGL - IN  - Number of global variables read
C     NVARNP - IN  - Number of nodal variables read
C     NVAREL - IN  - Number of element variables rea
C     NAMGV  - OUT - Returned character array with global variable names
C     NAMNV  - OUT - Returned character array with nodal variable names
C     NAMEV  - OUT - Returned character array with element variable names
C     IOERR   - OUT - error flag  
   
      character*(namlen) namgv(*), namnv(*), namev(*)
      integer ioerr
      
      IOERR = 0

C     Reads the names of the results variables from the database
C     call exgvan(fileid, vartyp, num_var, array_var_names(num_var), errorid)
      if (nvargl .gt. 0) then
         call exgvan(ndb, 'G', nvargl, namgv, ierr)
      end if
      if (nvarnp .gt. 0) then
         call exgvan(ndb, 'N', nvarnp, namnv, ierr)
      end if
      if (nvarel .gt. 0) then
         call exgvan(ndb, 'E', nvarel, namev, ierr)
      end if
            
      DO 130 I = 1, nvargl
         CALL EXUPCS (namgv(i))
 130  CONTINUE
      CALL PCKSTR (NVARGL, namgv)
      
      DO 140 I = 1, nvarnp
         CALL EXUPCS (namnv(i))
 140  CONTINUE
      CALL PCKSTR (NVARNP, namnv)
      
      DO 150 I = 1, nvarel
         CALL EXUPCS (namev(i))
 150  CONTINUE
      CALL PCKSTR (NVAREL, namev)
      
      return
      end
