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
      SUBROUTINE SWPESS (NUMESS, IDESS, NEESS, NEDSS,
     *  IXEESS, IXEDSS, LTEESS, LTSSS, LTSNC, FAC, ALLONE, NONQUD,
     *  COMTOP)
C=======================================================================

C   --Parameters:
C   --
C   --   NUMESS - IN/OUT - the number of element side sets
C   --   LESSEL - IN/OUT - the length of the element side sets element list
C   --   IDESS - IN/OUT - the element side set ID for each set
C   --   NEESS - IN/OUT - the number of elements for each set
C   --   NEDSS - IN/OUT - the number of dist-fac for each set
C   --   IXEESS - IN/OUT - the index of the first element for each set
C   --   IXEDSS - IN/OUT - the index of the first dist-fac for each set
C   --   LTEESS - IN/OUT - the elements for all sets
C   --   LTSSS - IN/OUT - the sides for all sets
C   --   LTSNC - IN/OUT - the face count for each element/side in the list
C   --   FACESS - IN/OUT - the distribution factors for all sets????????????

C   --   USESDF - IN - true if df are non-unity, false if all unity
C   --   NONQUAD - IN - true if model contains non-hex/non-quad elements
      
      INTEGER IDESS(*)   ! NUMESS
      INTEGER NEESS(*)   ! NUMESS
      INTEGER NEDSS(*)   ! NUMESS
      INTEGER IXEESS(*)  ! NUMESS
      INTEGER IXEDSS(*)  ! NUMESS
      INTEGER LTEESS(*)  ! LESSEL
      INTEGER LTSSS(*)   ! LESSEL
      INTEGER LTSNC(*)   ! LESSEL
      REAL    FAC(*)     ! LESSDL
      LOGICAL ALLONE, NONQUD
      CHARACTER*(*)   COMTOP

      IF (COMTOP(:3) .NE. 'SHE') THEN
        CALL PRTERR ('PROGRAM',
     *    'Swapping of sidesets on non-shell elements not supported')
        return
      END IF
      DO 30 I=1, NUMESS
         IF (IDESS(I) .LT. 0) THEN
            CALL SWPS1(NEESS(I), LTSSS(IXEESS(I)))
            IDESS(I) = -IDESS(I)
         END IF
   30 CONTINUE
      RETURN
      END

      subroutine swps1(numsid, side)
      
      integer numsid
      integer side(*)

      do 10 i = 1, numsid
        side(i) = 3 - side(i)
 10   continue
      return
      end
