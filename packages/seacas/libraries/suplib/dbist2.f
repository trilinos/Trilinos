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
      SUBROUTINE DBIST2 (NDB, NVAREL, NVARDM, NELBLK, NELBDM, ISEVOK,
     $     VAREL, NUMELB, IVAR, IELB, *)
C=======================================================================
C$Id: dbist2.f,v 1.5 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbist2.f,v $
CRevision 1.5  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.4  1992/04/08 21:13:22  gdsjaar
CFixed problem with singly accessing doubly dimensioned array
CAdded params to dbist2 and dbist1 so error messages would print
C
c Revision 1.3  1990/11/30  09:50:50  gdsjaar
c Modified to work on Unicos
c
c Revision 1.1.1.1  90/08/14  16:13:01  gdsjaar
c Testing
c 
c Revision 1.1  90/08/14  16:13:00  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:12  gdsjaar
c Initial revision
c 

C   --*** DBIST2 *** (EXOLIB) Internal to DBISTE, Read element variables 
C   --   Written by Greg Sjaardema 8/8/90, to remove MAX from dimensions
C   --
C   --DBIST2 reads the database element variables for one time step.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NVAREL - IN - the number of element variables
C   --   NVARDM - IN - the row dimension of VAREL
C   --   NELBLK - IN - the number of element blocks
C   --   NELBDM - IN - the row dimension of ISEVOK
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   VAREL - OUT - the element variables for the time step (if OPTION)
C   --   IVAR  - OUT - the nodal variable number if error on read.
C   --   IELB  - OUT - the element block number if error on read.
C   --   * - return statement if error encountered, including end-of-file;
C   --      message is printed
C   --
      INTEGER NDB
      INTEGER NVAREL, NVARDM, NELBLK, NELBDM
      INTEGER NUMELB(*)
      LOGICAL ISEVOK(NELBDM,*)
      REAL VAREL(NVARDM,*)

      IEL0 = 0
      DO 130 IELB = 1, NELBLK
         DO 120 IVAR = 1, NVAREL
            IF (ISEVOK(IELB,IVAR)) THEN
               READ (NDB, END=200, ERR=200, IOSTAT=IERR)
     &              (VAREL(IVAR,IEL0+N), N=1,NUMELB(IELB))
            END IF
  120    CONTINUE
         IEL0 = IEL0 + NUMELB(IELB)
  130 CONTINUE
      RETURN
  200 CONTINUE
      RETURN 1
      END
