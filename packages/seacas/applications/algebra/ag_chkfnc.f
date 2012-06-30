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

C=======================================================================
      SUBROUTINE CHKFNC (FNAM, ENAM)
C=======================================================================

C   --*** CHKFNC *** (ALGEBRA) Check that the function names match
C   --   Written by Amy Gilkey - revised 07/24/87
C   --
C   --CHKFNC prints a program error message if the two function names
C   --do not match.  This is a debugging tool to check that the functions
C   --in EVAL are ordered correctly.
C   --
C   --Parameters:
C   --   FNAM - IN - the name of the function from the EVAL code
C   --   ENAM - IN - the name of the function from the equation entry

      CHARACTER*(*) FNAM, ENAM

      CHARACTER*80 ERRMSG

      IF (FNAM .NE. ENAM) THEN
         ERRMSG = 'Function ' // FNAM(:LENSTR(FNAM))
     &      // ' should be ' // ENAM(:LENSTR(ENAM))
         CALL PRTERR ('PROGRAM', ERRMSG(:LENSTR(ERRMSG)))
      END IF

      RETURN
      END
