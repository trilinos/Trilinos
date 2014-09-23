C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

C $Id: skinit.f,v 1.1 1990/11/30 11:15:45 gdsjaar Exp $
C $Log: skinit.f,v $
C Revision 1.1  1990/11/30 11:15:45  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]SKINIT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SKINIT (STACK, NDIM, LENGTH, IERROR)
C***********************************************************************
C
C  SUBROUTINE SKINIT = STACK MANAGEMENT ROUTINE
C
C***********************************************************************
C
C** PARAMETERS
C      STACK  = STACK ARRAY
C      NDIM   = DIMENSIONED SIZE OF STACK IN CALLING PROGRAM
C      LENGTH = LENGTH OF STACK .LE. NDIM - 2
C      IERROR = 0 - NO ERROR
C               1 - STACK TOO SHORT (I.E. LENGTH > NDIM - 2
C               2 - STACK EMPTY
C               3 - STACK FULL
C               4 - INVALID STACK TYPE
C
C**********************************************************************
C
      PARAMETER (LSYOUT = 6)
      CHARACTER*(*) TYPE
      INTEGER STACK(NDIM)
C
      IF (NDIM .LT. LENGTH + 2) THEN
         IERROR = 1
      ELSE
         STACK(1) = 0
         STACK(2) = LENGTH
         IERROR = 0
      END IF
C
      RETURN
C
C=======================================================================
      ENTRY SKPOP (STACK, NDIM, IVALUE, IERROR)
C
      IF (STACK(1) .EQ. 0) THEN
         IERROR = 2
      ELSE
         IVALUE = STACK(STACK(1) + 2)
         STACK(1) = STACK(1) - 1
         IERROR = 0
      END IF
C
      RETURN
C
C=======================================================================
      ENTRY SKPUSH (STACK, NDIM, IVALUE, IERROR)
C
      IF (STACK(1) .EQ. STACK(2)) THEN
         IERROR = 3
      ELSE
         STACK(1) = STACK(1) + 1
         STACK(STACK(1) + 2) = IVALUE
         IERROR = 0
      END IF
C
      RETURN
C
C=======================================================================
      ENTRY SKEROR (LOUT, IERROR)
C
      IF (LOUT .EQ. 0) THEN
         LUNIT = LSYOUT
      ELSE
         LUNIT = LOUT
      END IF
C
      IF (IERROR .EQ. 0) THEN
      ELSE IF (IERROR .EQ. 1) THEN
         WRITE (LUNIT, '(A)') ' STACK ERROR:  ARRAY TOO SHORT'
      ELSE IF (IERROR .EQ. 2) THEN
         WRITE (LUNIT, '(A)') ' STACK ERROR:  STACK EMPTY'
      ELSE IF (IERROR .EQ. 3) THEN
         WRITE (LUNIT, '(A)') ' STACK ERROR:  STACK FULL'
      ELSE IF (IERROR .EQ. 4) THEN
         WRITE (LUNIT, '(A)') ' STACK ERROR:  INVALID TYPE'
      ELSE
         WRITE (LUNIT, '(A)') ' STACK ERROR:  UNKNOWN ERROR'
      END IF
C
      IERROR = 0
      RETURN
C
C=======================================================================
      ENTRY SKPRIN (LOUT, STACK, NDIM, TYPE, IERROR)
C
      IF (LOUT .EQ. 0) THEN
         LUNIT = LSYOUT
      ELSE
         LUNIT = LOUT
      END IF
C
      IF (STACK(1) .EQ. 0) THEN
         IERROR = 2
      ELSE IF (TYPE .EQ. 'I') THEN
         WRITE (LUNIT, '(2I8)') (I, STACK(I + 2),
     &      I = STACK(1), 1, -1)
      ELSE IF (TYPE .EQ. 'R') THEN
         CALL SKPRNT (LUNIT, STACK(1), STACK(1), NDIM)
      ELSE
         IERROR = 4
      END IF
C
      RETURN
C
      END
