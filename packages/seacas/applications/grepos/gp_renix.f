C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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
      SUBROUTINE RENIX (LLIST, IOFFIX, IX, LIST, ZEROIT)
C=======================================================================

C   --*** RENIX *** (GJOIN) Renumbers items according to an index list
C   --   Written by Amy Gilkey - revised 09/29/87
C   --
C   --RENIX renumbers each item in the list to the index given in another
C   --list or adds an offset to the item.
C   --
C   --Parameters:
C   --   LLIST - IN - the length of LIST
C   --   IOFFIX - IN - the offset to be added; if <0, use IX
C   --   IX - IN - the new index of the specified item
C   --   LIST - IN/OUT - the list of integers; renumbered
C   --   ZEROIT - IN - logical, if true and IX(LIST(I)) <= 0,
C               set LIST(I) = 0, else set to ABS value.
      INTEGER IX(*)
      INTEGER LIST(*)
      LOGICAL ZEROIT

      IF (IOFFIX .GE. 0) THEN
         DO 100 I = 1, LLIST
            LIST(I) = LIST(I) + IOFFIX
  100    CONTINUE
      ELSE
         DO 110 I = 1, LLIST
           IF (IX(LIST(I)) .LE. 0 .AND. ZEROIT) THEN
              LIST(I) = 0
           ELSE
              LIST(I) = ABS(IX(LIST(I)))
           END IF
  110    CONTINUE
      END IF

      RETURN
      END
