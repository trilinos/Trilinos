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
      SUBROUTINE DBINAM (NDB, OPTION, NDIM, NELBLK,
     &   NNDIM, NNELB, NVARHI, NVARGL, NVARNP, NVAREL,
     &   NAMECO, NAMELB, VNAMES, IXHV, IXGV, IXNV, IXEV,
     &   IA, KIEVOK, EXODUS, *)
C=======================================================================
C$Id: dbinam.f,v 1.5 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbinam.f,v $
CRevision 1.5  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.4  1998/03/24 15:24:55  gdsjaar
CRemoved several warnings, converted to ANSI C, cleaned up quite a bit
CCompiles without warnings on Sun, Dec, HP, gcc, g77.
C
CRevision 1.3  1992/07/16 22:38:08  gdsjaar
CCheck for zero number of variables to determine if exodus.  exo2exo1
Cwrites zeroes for all variables which screwed up previous method of
Cchecking.
C
c Revision 1.2  1990/10/11  17:05:03  gdsjaar
c Fixed first call to DBINM1 - added MAX()
c
c Revision 1.1.1.1  90/08/14  16:12:48  gdsjaar
c Testing
c 
c Revision 1.1  90/08/14  16:12:47  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:10  gdsjaar
c Initial revision
c 

C   --*** DBINAM *** (EXOLIB) Read database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --
C   --DBINAM reads the names of the coordinates, the element block types,
C   --and the database variables from the database.  All names are converted
C   --to uppercase and all embedded blanks within a name are removed.
C   --The element block variable truth table is also read.
C   --
C   --Note that the numbers of variables are read in this routine.
C   --
C   --This routine calls DBIV1 if the truth table is read.
C   --
C   --This routine calls DBVINI and uses DBVIX to get the variable name
C   --indices.
C   --
C   --Dynamic memory is reserved in this routine.  If there is a problem,
C   --the routine returns normally without printing an error message.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'C' to store coordinate names
C   --      'B' to store element block names
C   --      'V' to store variables names
C   --      'T' to store element block variable truth table
C   --   NDIM - IN - the number of coordinates per node
C   --   NELBLK - IN - the number of element blocks
C   --   NNDIM - OUT - the number of coordinates per node; <0 if end-of-file
C   --   NNELB - OUT - the number of element blocks; <0 if end-of-file
C   --   NVARHI - OUT - the number of history variables; <0 if end-of-file
C   --   NVARGL - OUT - the number of global variables; <0 if end-of-file
C   --   NVARNP - OUT - the number of nodal variables; <0 if end-of-file
C   --   NVAREL - OUT - the number of element variables; <0 if end-of-file
C   --   NAMECO - OUT - the names of the coordinates; max size = 6 (if OPTION)
C   --   NAMELB - OUT - the names of the element block types; max size = 256
C   --      (if OPTION)
C   --   VNAMES - OUT - the names of the global, nodal, and element variables;
C   --      max size = 256 (if OPTION)
C   --   IXHV - OUT - the VNAMES index of the history variable names (if OPTION)
C   --   IXGV - OUT - the VNAMES index of the global variable names (if OPTION)
C   --   IXNV - OUT - the VNAMES index of the nodal variable names (if OPTION)
C   --   IXEV - OUT - the VNAMES index of the element variable names (if OPTION)
C   --   A - OUT - the dynamic memory base array
C   --   KIEVOK - OUT - the dynamic memory index of the element block variable
C   --      truth table; (if OPTION)
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   EXODUS - OUT - false if GENESIS file, true if EXODUS file so far
C   --   * - return statement if error encountered, including end-of-file;
C   --      NOT used if valid GENESIS file; message is printed
C   --
C   --Database must be positioned in front of coordinate names upon entry;
C   --upon exit positioned after names.

C   --Routines Called:
C   --   EXUPCS - (SUPES) Convert to uppercase and blank non-standard
C   --   MDRSRV - (SUPES) Reserve dynamic memory
C   --   PCKSTR - (STRLIB) Remove embedded blanks

      PARAMETER (MAXDIM=6, MAXELB=256, MAXVAR=256)

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NDIM, NELBLK
      INTEGER NNDIM, NNELB
      CHARACTER*8 NAMECO(*)
      CHARACTER*8 NAMELB(*)
      CHARACTER*8 VNAMES(*)
      INTEGER NVARHI, NVARGL, NVARNP, NVAREL
      INTEGER IXHV, IXGV, IXNV, IXEV
      INTEGER IA(*)
      INTEGER KIEVOK
      LOGICAL EXODUS

      CHARACTER*80 ERRMSG
      INTEGER LDUM

      EXODUS = .FALSE.
      NNDIM = -999
      NNELB = -999
      NVARHI = -999
      NVARGL = -999
      NVARNP = -999
      NVAREL = -999

C   --Read and pack coordinate names

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
         IF (NDIM .GT. MAXDIM) CALL PRTERR ('WARNING',
     &      'Too many coordinate names in the database')

         READ (NDB, END=160, ERR=170, IOSTAT=IERR)
     &      (NAMECO(I), I=1,MIN(NDIM,MAXDIM))

         DO 100 I = 1, MIN(NDIM,MAXDIM)
            CALL EXUPCS (NAMECO(I))
  100    CONTINUE
         CALL PCKSTR (MIN(NDIM,MAXDIM), NAMECO)
      ELSE
         READ (NDB, END=160, ERR=160, IOSTAT=IERR)
      END IF
      NNDIM = NDIM

C   --Read and pack element block type names

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'B') .GT. 0)) THEN
         IF (NELBLK .GT. MAXELB) CALL PRTERR ('WARNING',
     &      'Too many element block names in the database')

         READ (NDB, END=160, ERR=180, IOSTAT=IERR)
     &      (NAMELB(I), I=1,MIN(NELBLK,MAXELB))

         DO 110 I = 1, MIN(NELBLK,MAXELB)
            CALL EXUPCS (NAMELB(I))
  110    CONTINUE
         CALL PCKSTR (MIN(NELBLK,MAXELB), NAMELB)
      ELSE
         READ (NDB, END=160, ERR=180, IOSTAT=IERR)
      END IF
      NNELB = NELBLK

C   --Read the number of variables

      READ (NDB, END=160, ERR=190, IOSTAT=IERR)
     &   NVARHI, NVARGL, NVARNP, NVAREL
      if ((nvarhi + nvargl + nvarnp + nvarel) .eq. 0) go to 160

      EXODUS = .TRUE.

C   --Initialize for DBVTYP and DBVIX

      CALL DBVINI (NVARHI, NVARGL, NVARNP, NVAREL)

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0)) THEN
         IF (NVARHI + NVARGL + NVARNP + NVAREL .GT. MAXVAR)
     &      CALL PRTERR ('WARNING',
     &      'Too many variable names in the database')

C      --Get the name indices

         CALL DBVIX ('H', 1, IXHV)
         CALL DBVIX ('H', NVARHI, IXHVE)
         IXHVE = MIN(IXHVE,MAXVAR)
         CALL DBVIX ('G', 1, IXGV)
         CALL DBVIX ('G', NVARGL, IXGVE)
         IXGVE = MIN(IXGVE,MAXVAR)
         CALL DBVIX ('N', 1, IXNV)
         CALL DBVIX ('N', NVARNP, IXNVE)
         IXNVE = MIN(IXNVE,MAXVAR)
         CALL DBVIX ('E', 1, IXEV)
         CALL DBVIX ('E', NVAREL, IXEVE)
         IXEVE = MIN(IXEVE,MAXVAR)

C      --Read and pack variable names

         READ (NDB, END=200, ERR=200, IOSTAT=IERR)
     &      (VNAMES(I), I=IXHV,IXHVE),
     &      (VNAMES(I), I=IXGV,IXGVE),
     &      (VNAMES(I), I=IXNV,IXNVE),
     &      (VNAMES(I), I=IXEV,IXEVE)

         DO 120 I = IXHV, IXHVE
            CALL EXUPCS (VNAMES(I))
  120    CONTINUE
         CALL PCKSTR (NVARHI, VNAMES(IXHV))

         DO 130 I = IXGV, IXGVE
            CALL EXUPCS (VNAMES(I))
  130    CONTINUE
         CALL PCKSTR (NVARGL, VNAMES(IXGV))

         DO 140 I = IXNV, IXNVE
            CALL EXUPCS (VNAMES(I))
  140    CONTINUE
         CALL PCKSTR (NVARNP, VNAMES(IXNV))

         DO 150 I = IXEV, IXEVE
            CALL EXUPCS (VNAMES(I))
  150    CONTINUE
         CALL PCKSTR (NVAREL, VNAMES(IXEV))

      ELSE
         READ (NDB, END=200, ERR=200, IOSTAT=IERR)
      END IF

C   --Read the element block variable truth table

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
         CALL MDRSRV ('ISEVOK', KIEVOK, NELBLK * NVAREL)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 160

         CALL DBINM1 (NDB, OPTION, NELBLK, NVAREL, IA(KIEVOK),
     *     IA(KIEVOK), IERR, MAX(NELBLK,1), *210)

C      --Call DBIV1

         CALL DBIV1 (NELBLK,
     &      NVARHI, NVARGL, NVARNP, NVAREL, IA(KIEVOK))

      ELSE
         CALL DBINM1 (NDB, OPTION, NELBLK, NVAREL, LDUM, IDUM,
     &      IERR, MAX(NELBLK,1), *210)
      END IF

  160 CONTINUE
      RETURN

  170 CONTINUE
      ERRMSG = 'COORDINATE NAMES'
      GOTO 220
  180 CONTINUE
      ERRMSG = 'ELEMENT BLOCK NAMES'
      GOTO 220
  190 CONTINUE
      ERRMSG = 'NUMBER OF VARIABLES'
      GOTO 220
  200 CONTINUE
      ERRMSG = 'VARIABLE NAMES'
      GOTO 220
  210 CONTINUE
      ERRMSG = 'ELEMENT BLOCK VARIABLE TRUTH TABLE'
      GOTO 220
  220 CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
