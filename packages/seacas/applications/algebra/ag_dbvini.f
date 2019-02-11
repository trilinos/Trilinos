C    Copyright(C) 2008-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C
C    * Neither the name of NTESS nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C
C=======================================================================
      SUBROUTINE DBVINI (NVARGL, NVARNP, NVAREL)
C=======================================================================

C   --*** DBVINI *** (EXOLIB) Initialize for DBVTYP and DBVIX
C   --   Written by Amy Gilkey - revised 10/14/87
C   --
C   --DBVINI initializes the indices for DBVTYP and DBVID.  It must be
C   --called before either of the other two routines are called.
C   --
C   --Note that the indices are shared because the other two routines
C   --are ENTRY routines of DBVINI.
C   --
C   --Parameters:
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables

      INTEGER NVARGL, NVARNP, NVAREL

      INTEGER IXGV, IXNV, IXEV, IXGVE, IXNVE, IXEVE
      common /dbv/ IXGV, IXNV, IXEV, IXGVE, IXNVE, IXEVE

      IXGV  = 1
      IXGVE = IXGV + NVARGL - 1
      IXNV  = IXGVE + 1
      IXNVE = IXNV + NVARNP - 1
      IXEV  = IXNVE + 1
      IXEVE = IXEV + NVAREL - 1

      RETURN
      END

C=======================================================================
      SUBROUTINE DBVTYP (IIX, TYP, ID)
C=======================================================================

C   --*** DBVTYP *** (EXOLIB) Return the variable type and number
C   --   Written by Amy Gilkey - revised 03/18/88
C   --
C   --DBVTYP is passed a variable index.  It returns the variable type
C   --and variable number.
C   --
C   --Note that DBVINI must be called before this routine to initialize
C   --the indices.  The indices are shared because this routine is an
C   --ENTRY routine of DBVINI.
C   --
C   --Parameters:
C   --   IIX - IN  - the variable index
C   --   TYP - OUT - the variable type: 'G'lobal, 'N'odal, 'E'lement
C   --   ID  - OUT - the variable number within the type

      INTEGER IXGV, IXNV, IXEV, IXGVE, IXNVE, IXEVE
      common /dbv/ IXGV, IXNV, IXEV, IXGVE, IXNVE, IXEVE

      INTEGER IIX, ID
      CHARACTER TYP

      IF ((IXGV .LE. 0) .AND. (IXNV .LE. 0)
     &    .AND. (IXEV .LE. 0)) RETURN

      IF ((IIX .GE. IXGV) .AND. (IIX .LE. IXGVE)) THEN
         TYP = 'G'
         ID = IIX - IXGV + 1
      ELSE IF ((IIX .GE. IXNV) .AND. (IIX .LE. IXNVE)) THEN
         TYP = 'N'
         ID = IIX - IXNV + 1
      ELSE IF ((IIX .GE. IXEV) .AND. (IIX .LE. IXEVE)) THEN
         TYP = 'E'
         ID = IIX - IXEV + 1
      ELSE
         TYP = ' '
         ID = 0
      END IF

      RETURN
      END

C=======================================================================
      SUBROUTINE DBVIX (ITYP, IID, IX)
C=======================================================================

C   --*** DBVIX *** (EXOLIB) Return the variable index
C   --   Written by Amy Gilkey - revised 10/14/87
C   --
C   --DBVIX is passed a variable type and number.  It returns the variable
C   --index.
C   --
C   --Note that DBVINI must be called before this routine to initialize
C   --the indices.  The indices are shared because this routine is an
C   --ENTRY routine of DBVINI.
C   --
C   --Parameters:
C   --   ITYP - IN  - the variable type: 'G'lobal, 'N'odal, 'E'lement
C   --   IID  - IN  - the variable number within the type
C   --   IX   - OUT - the variable index

      INTEGER IXGV, IXNV, IXEV, IXGVE, IXNVE, IXEVE
      common /dbv/ IXGV, IXNV, IXEV, IXGVE, IXNVE, IXEVE

      CHARACTER ITYP
      INTEGER IID, IX

      IF ((IXGV .LE. 0) .AND. (IXNV .LE. 0)
     &    .AND. (IXEV .LE. 0)) RETURN

      IF (ITYP .EQ. 'G') THEN
         IX = IID + IXGV - 1
      ELSE IF (ITYP .EQ. 'N') THEN
         IX = IID + IXNV - 1
      ELSE IF (ITYP .EQ. 'E') THEN
         IX = IID + IXEV - 1
      END IF

      RETURN
      END
