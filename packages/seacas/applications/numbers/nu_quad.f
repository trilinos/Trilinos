C    Copyright(C) 1988-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

C $Id: quad.f,v 1.1 1991/02/21 15:45:19 gdsjaar Exp $
C $Log: quad.f,v $
C Revision 1.1  1991/02/21 15:45:19  gdsjaar
C Initial revision
C
      SUBROUTINE QUAD(XXX, XI, XG, NDIM, NNODES, NQUAD, WT)
      DIMENSION XXX(NDIM+1,NNODES,NQUAD), XI(NDIM,*), XG(NDIM,*)
C
      IF (NQUAD .EQ. 1) THEN
          QUADL = 0.0
      ELSE
          QUADL = 1./SQRT(3.)
      END IF
C
      WT = 2.**NDIM / FLOAT(NQUAD)
      IF (NQUAD .EQ. 1) THEN
          XG(1,1) = 0.0
          XG(2,1) = 0.0
          XG(3,1) = 0.0
      ELSE
          DO 20 I=1, NNODES
              DO 10 J=1, NDIM
                  XG(J,I) = XI(J,I) * QUADL
   10         CONTINUE
   20     CONTINUE
      END IF
C
      IF (NDIM .EQ. 3) THEN
          DO 40 I=1, NQUAD
              DO 30 J=1, NNODES
                  TMP1 = (1. + XI(1,J) * XG(1,I))
                  TMP2 = (1. + XI(2,J) * XG(2,I))
                  TMP3 = (1. + XI(3,J) * XG(3,I))
C
                  XXX(1,J,I) = TMP1    * TMP2 * TMP3 / 8.0
                  XXX(2,J,I) = XI(1,J) * TMP2 * TMP3 / 8.0
                  XXX(3,J,I) = XI(2,J) * TMP1 * TMP3 / 8.0
                  XXX(4,J,I) = XI(3,J) * TMP1 * TMP2 / 8.0
   30         CONTINUE
   40     CONTINUE
      ELSE
          DO 60 I=1, NQUAD
              DO 50 J=1, NNODES
                  TMP1 = (1. + XI(1,J) * XG(1,I))
                  TMP2 = (1. + XI(2,J) * XG(2,I))
C
                  XXX(1,J,I) = TMP1    * TMP2 / 4.0
                  XXX(2,J,I) = XI(1,J) * TMP2 / 4.0
                  XXX(3,J,I) = XI(2,J) * TMP1 / 4.0
   50         CONTINUE
   60     CONTINUE
      END IF
      RETURN
      END
