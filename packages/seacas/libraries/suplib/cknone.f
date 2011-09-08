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
      SUBROUTINE CKNONE (NVAL, ISSEL, VALNAM, *)
C=======================================================================

C   --*** CKNONE *** (ETCLIB) Check number of values is zero
C   --   Written by Amy Gilkey - revised 12/23/87
C   --
C   --CKNONE prints an error message if the number of values is zero.
C   --
C   --Parameters:
C   --   NVAL - IN - the value being checked
C   --   ISSEL - IN - print none selected error message iff true
C   --   VALNAM - IN - the name of the value being checked (plural)
C   --   * - return statement if error

      LOGICAL ISSEL
      CHARACTER*(*) VALNAM

      CHARACTER*1024 ERRMSG

      IF (NVAL .LE. 0) THEN
         IF (ISSEL) THEN
            ERRMSG = 'No ' // VALNAM // ' are selected'
            CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
         ELSE
            ERRMSG = 'There are no ' // VALNAM
            CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
         END IF
         RETURN 1
      END IF

      RETURN
      END
