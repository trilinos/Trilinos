C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C 
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
C 

C $Id: srchge.f,v 1.3 2007/10/17 18:40:36 gdsjaar Exp $
C $Log: srchge.f,v $
C Revision 1.3  2007/10/17 18:40:36  gdsjaar
C Added copyright notice to all files.
C
C Mapvar is licensed under the BSD license
C
C Revision 1.2  2000/11/14 17:30:40  gdsjaar
C Removed old cray compiler directives
C
C Revision 1.1  1998/03/13 18:12:28  gdsjaar
C New code -- mapvar. Interpolates results form an exodusII results file
C to a differently mesh geometry.  Written by Gerry Wellman,
C 9117. Loosely based on MERLIN. Provides a superset of merlin
C functionality.
C
C
      SUBROUTINE SRCHGE( LBLK,NE,X,IND,XV,IMIN,IMAX,NDIM,I,
     *                  IL,IU,IT,INDX1,INDX2,XTST )
C
C-----------------------------------------------------------------------
C     
C DESCRIPTION:
C
C  PERFORM A BINARY SEARCH TO FIND THE ELEMENT NUMBER I
C  OF A SORTED ARRAY FOR WHICH ALL ELEMENTS AT I OR ABOVE ARE
C  GREATER OR EQUAL TO THAN SOME VALUE XV,
C  WHILE ALL ELEMENTS BELOW I ARE LESS THAN THE VALUE.
C
C       X(I-2)     X(I-1)    X(I)    X(I+1)   X(I+2)
C                      XV                             X>
C
C  ASSUMED THAT THE ARRAY X HAS BEEN SORTED IN INCREASING ORDER,
C  BUT THE ELEMENTS HAVE NOT BEEN MOVED.
C  THE SORTED LIST IS DETERMINED BY THE ARRAY INDX,
C  WHICH POSITIONS THE ORIGINAL UNSORTED X ARRAY ELEMENTS
C  IN THE SORTED LIST. THUS, THE 5TH ELEMENT IN THE SORTED LIST IS
C    X(IND(5))
C
C-----------------------------------------------------------------------
C
C  INPUT:
C
C  X      -  ARRAY IN UNSORTED ORDER
C  IND    -  INDEX ARRAY GIVING THE ELEMENT ORDER IN THE SORTED LIST
C  XV     -  X VALUE TO TEST AGAINST
C  IMIN   -  THE LOWEST NUMBERED POSITION IN THE SORTED LIST TO TEST
C  IMAX   -  THE HIGHEST NUMBERED POSITION IN THE SORTED LIST TO TEST
C  NDIM   -  THE DIMENSION OF THE ARRAYS
C
C  OUTPUT:
C
C  I      -  THE FIRST POSITION IN THE SORTED LIST .GT. XV
C
C  SCRATCH:
C
C  IL
C  IU
C  IT
C  INDX1
C  INDX2
C  XTST
C
C-----------------------------------------------------------------------
C
      DIMENSION
     *  X(NDIM), IND(NDIM), XV(LBLK), I(LBLK)
      DIMENSION
     *  IL(LBLK),IU(LBLK),IT(LBLK),INDX1(LBLK),INDX2(LBLK),
     *  XTST(LBLK)
C
      IF (IMAX.EQ.0.OR.NDIM.EQ.0) THEN
         DO J = 1, NE
            I(J) = 0
         ENDDO
         RETURN
      ENDIF
      DO 25 J = 1, NE
        IL(J) = IMIN
        IU(J) = IMAX
        INDX1(J) = J
 25   CONTINUE
        ILOOP = NE
 1000 CONTINUE
C
      DO 50 JJ = 1, ILOOP
        INDX2(JJ) = INDX1(JJ)
        J = INDX1(JJ)      
        IT(J) =  (IU(J) + IL(J)) / 2 
 50   CONTINUE
      DO 35 J = 1, NE
       XTST(J) = X( IND(IT(J)) )
 35   CONTINUE
C
      IF ( ILOOP .GT. 64) THEN

      ILP = 0
      DO 60 JJ = 1, ILOOP
        J = INDX2(JJ)
        IF( XTST(J) .LT. XV(J) )THEN
          IL(J) =IT(J) + 1
        ELSE
          IU(J) =IT(J) - 1
        ENDIF
        IF( IL(J) .LE. IU(J)) THEN
          ILP = ILP + 1
          INDX1(ILP) = J 
        ENDIF
 60   CONTINUE
      ELSE

      ILP = 0
      DO 51 JJ = 1, ILOOP
        J = INDX2(JJ)      
        IF( XTST(J) .LT. XV(J) )THEN
          IL(J) =IT(J) + 1
        ELSE
          IU(J) =IT(J) - 1
        ENDIF
        IF( IL(J) .LE. IU(J)) THEN
          ILP = ILP + 1
          INDX1(ILP) = J 
        ENDIF
 51   CONTINUE
      ENDIF
c
      ILOOP = ILP
c
      IF(ILOOP .NE. 0 )GO TO 1000
C  RANGE HAD NARROWED TO 1 LOCATION. HOWEVER, THE POINT LAST TESTED
C  COULD BE ABOVE, BELOW, OR ON THE SEARCH POINT. CHECK FOR PROPER CASE
      DO 200 J = 1, NE
        IF( XTST(J) .LT. XV(J) )THEN
          I(J) = IT(J) + 1
        ELSE
          I(J) = IT(J) 
        ENDIF
 200  CONTINUE
C
      RETURN
      END
C
