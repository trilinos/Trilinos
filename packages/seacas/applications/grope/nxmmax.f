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
      SUBROUTINE NXMMAX (ISTEP, MMNUM, IVAR, NVAR, NUM, VAR,
     &   OLDMIN, OLDMAX, XMIN, XMAX,
     &   MINSTE, MAXSTE, MINNE, MAXNE, NMIN, NMAX)
C=======================================================================

C   --*** NXMMAX *** (GROPE) Find the min/max for this time step
C   --
C   --NXMMAX finds an incremental minimum and maximum of a variable
C   --for this time step.  It does not reset the min/max if not on the
C   --first step.  It may find the min/max that are greater/less than
C   --the previous values.
C   --
C   --Parameters:
C   --   ISTEP - IN - the current step number, previous XMIN and XMAX
C   --      are overwritten if <=1
C   --   MMNUM - IN - number of sequential min/max requests for this
C   --      variable, if >1 then find next min/max
C   --   IVAR - IN - min/max variable number
C   --   NUM - IN - the number of nodes or elements or 1 for global
C   --   VAR - IN - the variables for current time step
C   --   OLDMIN, OLDMAX - IN/OUT - the last minimum and maximum values,
C   --      only if MMNUM > 1
C   --   XMIN, XMAX - IN/OUT - the minimum and maximum values
C   --   MINSTE, MAXSTE - IN/OUT - the step number where the minimum and
C   --      maximum were found
C   --   MINNE, MAXNE - IN/OUT - the node or element number where the
C   --      minimum and maximum were found
C   --   NMIN, NMAX - IN/OUT - the number of values equal to the minimum
C   --      and maximum were found

      REAL VAR(NUM, NVAR)

      IF (ISTEP .LE. 1) THEN
         XMIN = 1.0E36
         XMAX = - 1.0E36
         IF (MMNUM .LE. 1) THEN
            OLDMIN = - 1.0E36
            OLDMAX = 1.0E36
         ELSE
            OLDMIN = XMAX
            OLDMAX = XMIN
         END IF
      END IF

      DO 100 I = 1, NUM
         IF ((XMIN .GT. VAR(I,IVAR))
     &      .AND. (VAR(I,IVAR) .GT. OLDMIN)) THEN
            XMIN = VAR(I,IVAR)
            MINSTE = ISTEP
            MINNE = I
            NMIN = 1
         ELSE IF (XMIN .EQ. VAR(I,IVAR)) THEN
            NMIN = NMIN + 1
         END IF

         IF ((XMAX .LT. VAR(I,IVAR))
     &      .AND. (VAR(I,IVAR) .LT. OLDMAX)) THEN
            XMAX = VAR(I,IVAR)
            MAXSTE = ISTEP
            MAXNE = I
            NMAX = 1
         ELSE IF (XMAX .EQ. VAR(I,IVAR)) THEN
            NMAX = NMAX + 1
         END IF
  100 CONTINUE

      RETURN
      END
