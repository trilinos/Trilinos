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
      SUBROUTINE CKESS (NUMESS, LESSEL, LESSNL, NUMEL, NUMNP,
     &   IDESS, NEESS, NNESS, IXEESS, IXNESS,
     &   LTEESS, LTSESS, FACESS, ICHECK, NDIM)
C=======================================================================

C   --*** CKESS *** (GROPE) Check database element side sets
C   --
C   --CKESS checks the element side set information.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NUMESS - IN - the number of element side sets
C   --   LESSEL - IN - the number of elements for all sets
C   --   LESSNL - IN - the number of nodes for all sets
C   --   NUMEL - IN - the number of elements
C   --   NUMNP - IN - the number of nodes
C   --   IDESS - IN - the element side set ID for each set
C   --   NEESS - IN - the number of elements for each set
C   --   NNESS - IN - the number of nodes for each set
C   --   IXEESS - IN - the index of the first element for each set
C   --   IXNESS - IN - the index of the first node for each set
C   --   LTEESS - IN - the elements for all sets
C   --   LTSESS - IN - the element faces for all sets
C   --   FACESS - IN - the distribution factors for all sets
C   --   ICHECK - SCRATCH - size = MAX (NUMESS, LESSEL, LESSNL)

      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER IXEESS(*)
      INTEGER IXNESS(*)
      INTEGER LTEESS(*)
      INTEGER LTSESS(*)
      REAL FACESS(*)
      INTEGER ICHECK(*)

      CHARACTER*128 STRA

C   --Check for unique identifier

      DO 100 IESS = 1, NUMESS
         IF (LOCINT (IDESS(IESS), IESS-1, IDESS) .GT. 0) THEN
            CALL INTSTR (1, 0, IDESS(IESS), STRA, LSTRA)
            CALL PRTERR ('CMDSPEC', 'Element side set ID '
     &         // STRA(:LSTRA) // ' is not unique')
         END IF
  100 CONTINUE

C   --Check number of elements in element side sets

      NESS = 0
      DO 110 IESS = 1, NUMESS
         NESS = MAX (NESS, IXEESS(IESS) + NEESS(IESS) - 1)
  110 CONTINUE

      IF (NESS .NE. LESSEL) THEN
         CALL PRTERR ('CMDSPEC', 'Maximum element index'
     &      // ' in all element side sets does not match total')
      END IF

C   --Check number of nodes in element side sets

      NESS = 0
      DO 120 IESS = 1, NUMESS
         NESS = MAX (NESS, IXNESS(IESS) + NNESS(IESS) - 1)
  120 CONTINUE

      IF (NESS .NE. LESSNL .AND. LESSNL .GT. 0) THEN
         CALL PRTERR ('CMDSPEC', 'Maximum node index'
     &      // ' in all element side sets does not match total')
      END IF

C   --Check all elements in element side sets are within element range

      CALL CHKRNG (LTEESS, LESSEL, NUMEL, NZERO, NERR)
      IF (NERR .GT. 0) THEN
         CALL PRTERR ('CMDSPEC',
     &      'Element side set element ids are out of range')
      END IF
      IF (NZERO .GT. 0) THEN
         CALL PRTERR ('CMDSPEC',
     &      'Element side set element ids are zero')
      END IF

C   --Check all element faces in element side sets are within range
C ... Since we don't know (or don't want to expend the effort...) the
C     the number of faces for each element, we assume that the maximum
C     number of faces is 4 for 2D and 6 for 3D      
      if (ndim .eq. 2) then
        CALL CHKRNG (LTSESS, LESSEL, 4, NZERO, NERR)
      else if (ndim .eq. 3) then
        CALL CHKRNG (LTSESS, LESSEL, 6, NZERO, NERR)
      end if
      
C ... Check for duplicate element/sides in a sideset. This causes
C     problems with some analysis codes
      do iess = 1, numess
        call iniint(numel, 0, icheck)
        nel = neess(iess)
        indx = ixeess(iess)
        do j = 0, nel-1
          iel = lteess(indx+j)
          ifa = ltsess(indx+j)
          if (btest(icheck(iel), ifa)) then
            write (stra, 10000) iel, ifa, idess(iess)
10000       FORMAT('SIDESET ERROR: The element face pair ',I10,'.',I1,
     $        ' is duplicated in sideset ', I10,'.')
            call sqzstr(stra, lstra)
            CALL PRTERR ('CMDSPEC', STRA(:lstra))
          else
            icheck(iel) = ibset(icheck(iel), ifa)
          end if
        end do
      end do

      IF (NERR .GT. 0) THEN
        CALL PRTERR ('CMDSPEC',
     &    'Element side set faces are out of range')
      END IF
      IF (NZERO .GT. 0) THEN
        CALL PRTERR ('CMDSPEC',
     &    'Element side set faces are zero')
      END IF
      RETURN
      END
