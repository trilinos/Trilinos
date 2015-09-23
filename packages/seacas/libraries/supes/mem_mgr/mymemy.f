C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
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
C    * Neither the name of Sandia Corporation nor the names of its
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
      SUBROUTINE MYMEMY( MEMREQ, LOCBLK, MEMRTN, MAXSIZ )
      SAVE NUSED
      DATA NUSED /0/
C
************************************************************************
C
C     FORTRAN EXTENSION LIBRARY - ANSI FORTRAN - USER INTERFACE ROUTINE
C
C     DESCRIPTION:
C     This routine requests the operating system to allocate or release
C     numeric storage. A postive MEMREQ indicates a request for memory,
C     while a negative MEMREQ indicates a release. All locations and
C     sizes are measured in numeric storage units.
C
C     In memory request mode, MEMRTN .LT. MEMREQ indicates an error.
C
C     In memory release mode, MEMRTN .LE. -MEMREQ. Furthermore, memory
C     must be released from the top down, i.e., LOCBLK must not change.
C
C     This version actually allocates storage from a static pool, whose
C     size is defined by the parameter MAXSIZ. If system dependent
C     support for the function IXLNUM is not implemented, the PARAMETER
C     and COMMON statements above must be duplicated in the caller.
C
C     FORMAL PARAMETERS:
C     MEMREQ    INTEGER         Number of numeric units
C     LOCBLK    INTEGER         Location of memory block
C     MEMRTN    INTEGER         Size of memory block at routine completion
C     MAXSIZ    INTEGER         Size of character memory - dimension in
C                               MDINIT.
C
C     SAVED VARIABLES:
C     NUSED     INTEGER         Number of units dynamically allocated
C
************************************************************************
C
      IF ( MEMREQ .GE. 0 ) THEN
C
C Allocate storage -
         LOCBLK = 1 + NUSED
         MEMRTN = MIN( MAXSIZ-NUSED , MEMREQ )
         NUSED = NUSED + MEMRTN
      ELSE
C
         MEMRTN = -MEMREQ
      END IF
C
      RETURN
      END
