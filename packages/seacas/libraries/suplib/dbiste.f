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
      SUBROUTINE DBISTE (NDB, OPTION, ISTEP,
     &   NVARHI, NVARGL, NVARNP, NUMNP, NVAREL, NELBLK, NUMELB, ISEVOK,
     &   TIME, WHOTIM, VARHI, VARGL, VARNP, VAREL, *)
C=======================================================================
C$Id: dbiste.f,v 1.6 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbiste.f,v $
CRevision 1.6  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.5  1997/03/20 19:40:12  caforsy
CUpdated Imakefile for Imake 6.1.  Changed printing routines to handle
Clarger problems.
C
CRevision 1.4  1992/04/08 21:13:20  gdsjaar
CFixed problem with singly accessing doubly dimensioned array
CAdded params to dbist2 and dbist1 so error messages would print
C
c Revision 1.3  1990/11/30  09:50:52  gdsjaar
c Modified to work on Unicos
c
c Revision 1.1.1.1  90/08/14  16:13:04  gdsjaar
c Testing
c 
c Revision 1.1  90/08/14  16:13:03  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:12  gdsjaar
c Initial revision
c 

C   --*** DBISTE *** (EXOLIB) Read database variables for one time step
C   --   Written by Amy Gilkey - revised 10/14/87
C   --
C   --DBISTE reads the database global, nodal, and element variables
C   --for one time step.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'H' to store history variables
C   --      'G' to store global variables
C   --      'E' to store element variables
C   --      'N' to store nodal variables
C   --   ISTEP - IN - the time step number
C   --   NVARHI - IN - the number of history variables
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NUMNP - IN - the number of nodes
C   --   NVAREL - IN - the number of element variables
C   --   NELBLK - IN - the number of element blocks
C   --   NUMELB - IN - the number of elements per block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   TIME - OUT - the time step time
C   --   WHOTIM - OUT - true iff whole (versus history) time step
C   --   VARHI - OUT - the history variables for the time step (if OPTION)
C   --   VARGL - OUT - the global variables for the time step (if OPTION)
C   --   VARNP - OUT - the nodal variables for the time step (if OPTION)
C   --   VAREL - OUT - the element variables for the time step (if OPTION)
C   --   * - return statement if error encountered, including end-of-file;
C   --      message is printed
C   --
C   --Database must be positioned in front of time step upon entry;
C   --upon exit positioned after time step.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER ISTEP
      INTEGER NVARHI, NVARGL, NVARNP, NUMNP, NVAREL, NELBLK
      INTEGER NUMELB(*)
      LOGICAL ISEVOK(*)
      REAL TIME
      LOGICAL WHOTIM
      REAL VARHI(*)
      REAL VARGL(*)
C     --NOTE: VARNP and VAREL are passed into DBIST? as doubly-dimensioned
      REAL VARNP(*)
      REAL VAREL(*)

      CHARACTER*80 ERRMSG

C   --Read step time

      READ (NDB, END=220, ERR=160, IOSTAT=IERR) TIME, HISTFL
      WHOTIM = (HISTFL .EQ. 0.0)

C   --Read history variables

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'H') .GT. 0)) THEN
         READ (NDB, END=170, ERR=170, IOSTAT=IERR)
     &      (VARHI(IVAR), IVAR=1,NVARHI)
      ELSE
         READ (NDB, END=170, ERR=170, IOSTAT=IERR)
      END IF

      IF (WHOTIM) THEN

C      --Read global variables

         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'G') .GT. 0)) THEN
            READ (NDB, END=180, ERR=180, IOSTAT=IERR)
     &         (VARGL(IVAR), IVAR=1,NVARGL)
         ELSE
            READ (NDB, END=180, ERR=180, IOSTAT=IERR)
         END IF

C      --Read nodal variables

         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
            CALL DBIST1 (NDB, NVARNP, MAX(NVARNP,1), NUMNP, VARNP, IVAR,
     &           *190)
         ELSE
            DO 110 IVAR = 1, NVARNP
               READ (NDB, END=190, ERR=190, IOSTAT=IERR)
  110       CONTINUE
         END IF

C      --Read element variables

         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'E') .GT. 0)) THEN
            CALL DBIST2 (NDB, NVAREL, MAX(NVAREL,1), NELBLK,
     $           MAX(NELBLK,1), ISEVOK, VAREL, NUMELB,
     &           IVAR, IELB, *200)
         ELSE
            DO 150 IELB = 1, NELBLK
               DO 140 IVAR = 1, NVAREL
                  IF (ISEVOK( (IVAR-1)*NELBLK+IELB )) THEN
                     READ (NDB, END=200, ERR=200, IOSTAT=IERR)
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      END IF

      RETURN

  160 CONTINUE
      WRITE (ERRMSG, '(A, I8)', IOSTAT=IDUM)
     &   'TIME for TIME STEP', ISTEP
      GOTO 210
  170 CONTINUE
      WRITE (ERRMSG, '(A, I8)', IOSTAT=IDUM)
     &   'HISTORY VARIABLES for TIME STEP', ISTEP
      GOTO 210
  180 CONTINUE
      WRITE (ERRMSG, '(A, I8)', IOSTAT=IDUM)
     &   'GLOBAL VARIABLES for TIME STEP', ISTEP
      GOTO 210
  190 CONTINUE
      WRITE (ERRMSG, '(A, I8, A, I8)', IOSTAT=IDUM)
     &   'NODAL VARIABLE', IVAR, ' for TIME STEP', ISTEP
      GOTO 210
  200 CONTINUE
      WRITE (ERRMSG, '(A, I8, A, I8, A, I8)', IOSTAT=IDUM)
     &   'ELEMENT VARIABLE', IVAR, ' of BLOCK', IELB,
     &   ' for TIME STEP', ISTEP
      GOTO 210
  210 CONTINUE
      CALL DBERR (IERR, ERRMSG)
  220 CONTINUE
      RETURN 1
      END
