C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
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

C $Id: unique.f,v 1.1 1991/02/21 15:46:09 gdsjaar Exp $
C $Log: unique.f,v $
C Revision 1.1  1991/02/21 15:46:09  gdsjaar
C Initial revision
C
      SUBROUTINE UNIQUE (LSTSN, NSEG, MAP, ITMP, NUMNIQ, NUMNP)
C
C***********************************************************************
C
C     DESCRIPTION:
C       This routine determines the number of unique node numbers in
C       a side set.
C
C     FORMAL PARAMETERS:
C       LSTSN   INTEGER   List of nodes on this boundary
C       NSEG    INTEGER   Number of nodes in side set
C       MAP     INTEGER   Relates node in side set to list of unique nodes
C       ITMP    INTEGER   Temporary array for sorting nodes
C       NUMNIQ  INTEGER   Number of unique nodes
C       NDIM    INTEGER   Number of spatial dimensions
C
C     CALLED BY:
C
C***********************************************************************
C
      DIMENSION LSTSN(*), MAP(*), ITMP(*)
C
      CALL INIINT (NSEG, 0, MAP)
      CALL INIINT (NUMNP, 0, ITMP)
C
      NUMNIQ = 0
      DO 30 I = 1 , NSEG
          IF ( ITMP(LSTSN(I)) .EQ. 0 ) THEN
              NUMNIQ = NUMNIQ + 1
              ITMP(LSTSN(I)) = NUMNIQ
              MAP(I) = NUMNIQ
          ELSE
              MAP(I) = ITMP(LSTSN(I))
          END IF
   30 CONTINUE
C
      RETURN
      END
