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
      SUBROUTINE DBISTE (NDB, OPTION, ISTEP, NELBLK, TIME,
     &                   NVARGL, NVARNP, NVAREL, NUMNP,
     &                   IDELB, NUMELB, ISEVOK,
     &                   VARGL, VARNP, VAREL, IOERR)
C=======================================================================
C$Id: dbiste.f,v 1.4 2007/10/17 18:46:09 gdsjaar Exp $
C$Log: dbiste.f,v $
CRevision 1.4  2007/10/17 18:46:09  gdsjaar
CAdded copyright notice to all files.
C
Cexotxt2 is licensed under the BSD license
C
CRevision 1.3  1996/05/21 16:52:20  caforsy
CAdded read/write for property data.  Cleaned up exodusII error checks
C
CRevision 1.2  1995/11/07 15:01:39  gdsjaar
CInitial checkin of ACCESS/translate/exotxt2
C

C   --*** DBISTE *** (EXOLIB) Read database variables for one time step
C   --   Written by Amy Gilkey - revised 10/14/87
C   --   Modified for ExodusIIv2 database format 10/16/95
C   --
C   --DBISTE reads the database global, nodal, and element variables
C   --for one time step.
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   OPTION - IN  - ' ' to not store
C   --                  '*' to store all, else store options:
C   --                  'G' to store global variables
C   --                  'E' to store element variables
C   --                  'N' to store nodal variables
C   --   ISTEP  - IN  - the time step number
C   --   NELBLK - IN  - the number of element blocks
C   --   TIME   - OUT - the time step time
C   --   NVARGL - IN  - the number of global variables
C   --   NVARNP - IN  - the number of nodal variables
C   --   NVAREL - IN  - the number of element variables
C   --   NUMNP  - IN  - the number of nodes
C   --   IDELB  - IN  - the element block ID's
C   --   NUMELB - IN  - the number of elements per block
C   --   ISEVOK - IN  - the element block variable truth table;
C   --                  variable i of block j exists iff ISEVOK(j,i)
C   --   VARGL  - OUT - the global variables for the time step (if OPTION)
C   --   VARNP  - OUT - the nodal variables for the time step (if OPTION)
C   --   VAREL  - OUT - the element variables for the time step (if OPTION)
C   --   IOERR  - OUT - I/O error flag


      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER ISTEP
      INTEGER NELBLK 
      REAL    TIME
      INTEGER NVARGL, NVARNP, NVAREL, NUMNP
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      LOGICAL ISEVOK(*)
      REAL    VARGL(*)
      REAL    VARNP(*)
      REAL    VAREL(*)
      INTEGER IOERR
      INTEGER EBID, NEPB

C   --Read step time

      CALL EXGTIM(NDB, ISTEP, TIME, IERR)

C     --Read global variables
      IF (NVARGL .GT. 0) THEN
         CALL EXGGV (NDB, ISTEP, NVARGL, VARGL, IERR)
      END IF


C     --Read nodal variables
      IF (NVARNP .GT. 0) THEN
         INP = 1
         DO 10 I = 1, NVARNP
            CALL EXGNV (NDB, ISTEP, I, NUMNP, VARNP(INP), IERR)
            INP = INP + NUMNP
 10      CONTINUE
      END IF

C      --Read element variables

      IF (NVAREL .GT. 0) THEN
         IEL = 1
         DO 20 I = 1, NELBLK
            EBID = IDELB(I)
            NEPB = NUMELB(I)
            DO 30 J = 1, NVAREL
               CALL EXGEV(NDB,ISTEP,J,EBID,NEPB,VAREL(IEL),IERR)
               IEL = IEL + NEPB
 30         CONTINUE
 20      CONTINUE
      END IF

      RETURN
      END
