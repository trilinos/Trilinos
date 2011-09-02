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

C=======================================================================
      SUBROUTINE DBIVTT (NDB, A, ISEVOK, ITMP, NELBLK, NVAREL,
     &                   NELBDM)
C=======================================================================
C$Id: dbivtt.f,v 1.4 2007/10/17 18:46:09 gdsjaar Exp $
C$Log: dbivtt.f,v $
CRevision 1.4  2007/10/17 18:46:09  gdsjaar
CAdded copyright notice to all files.
C
Cexotxt2 is licensed under the BSD license
C
CRevision 1.3  1996/05/21 16:52:21  caforsy
CAdded read/write for property data.  Cleaned up exodusII error checks
C
CRevision 1.2  1995/11/07 15:01:42  gdsjaar
CInitial checkin of ACCESS/translate/exotxt2
C
C   --*** DBIVTT *** Read element variable truth table
C   --   Modified for ExodusII format 8/26/95
C   --*** DBINAM *** (EXOLIB) Read database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --
C   --DBINAM performed a number of different input file read base
C   --on the passed in option argument.  DBINAM was split up
C   --into a number of different subroutins

C   --DBINAM reads the names of the coordinates, the element block types,
C   --and the database variables from the database.  All names are converted
C   --to uppercase and all embedded blanks within a name are removed.
C   --The element block variable truth table is also read.
C   --
C   --Note that the numbers of variables are read in this routine.
C   --
C   --This routine calls DBVINI and uses DBVIX to get the variable name
C   --indices.
C   --
C   --Dynamic memory is reserved in this routine.  If there is a problem,
C   --MEMERR is called and then this routine returns with and error IOERR=1
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   A      - OUT - the dynamic memory base array
C   --   NELBLK - IN  - the number of element blocks
C   --   NVAREL - IN  - the number of element variables; <0 if end-of-file
C   --   ISEVOK - OUT - the dynamic memory array of the element block variable
C   --                  truth table;variable i,block j exists iff ISEVOK(j,i)

C   --Routines Called:
C   --   MDRSRV - (SUPES) Reserve dynamic memory

      INTEGER NDB
      DIMENSION A(*)
      INTEGER NELBLK, NVAREL
      LOGICAL ISEVOK(NELBDM,*)
      INTEGER ITMP(NVAREL, NELBDM)

C     Read the element block variable truth table
C     call exgvtt(fileid, num_elem_blks, num_elem_var,
C                   isevok(num_elem_var, num_elem_blks, errorid)
C       isevok - num_elem_var cycles faster
      CALL EXGVTT(NDB, NELBLK, NVAREL, ITMP, IERR)
      DO 110 I = 1, NVAREL
         DO 100 IELB = 1, NELBLK
            ISEVOK(IELB,I) = (ITMP(I,IELB) .NE. 0)
 100     CONTINUE
 110  CONTINUE

      RETURN
      END

