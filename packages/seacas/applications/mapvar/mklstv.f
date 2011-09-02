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

C $Id: mklstv.f,v 1.2 2007/10/17 18:40:35 gdsjaar Exp $
C $Log: mklstv.f,v $
C Revision 1.2  2007/10/17 18:40:35  gdsjaar
C Added copyright notice to all files.
C
C Mapvar is licensed under the BSD license
C
C Revision 1.1  1998/03/13 18:12:24  gdsjaar
C New code -- mapvar. Interpolates results form an exodusII results file
C to a differently mesh geometry.  Written by Gerry Wellman,
C 9117. Loosely based on MERLIN. Provides a superset of merlin
C functionality.
C
C
      SUBROUTINE MKLSTV( NUMPTS,IND,IRNK2,IUP,ILO,INDX,
     *                   IE,LIST,NLIST,NBLKSZ,NSPC)
C     
C-----------------------------------------------------------------------
C     
C DESCRIPTION:
C
C VECTOR MAKE LIST (3D) 
C GIVEN A LIST OF PARTICLES (IE  THEIR INDEX AND RANK) FIND
C THE LIST OF PARTICLES WITHIN THE BOUNDS SET BY XMIN AND XMAX
C FOR THE IE'TH PARTICLE IN THE VECTOR BLOCK
C     
C-----------------------------------------------------------------------
C
C  CALLING ARGUMENTS:
C
C     NUMPTS INTEGER   NUMBER OF POINTS TO BE SEACHED
C     IND    INTEGER   ORDER INDEX
C     IRNK2  INTEGER   RANK 
C     IUP    INTEGER   SCRACTH (NBLKSZ LONG)
C     ILO    INTEGER   SCRACTH (NBLKSZ LONG)
C     INDX   INTEGER   SCRATCH (NBLKSZ LONG)
C     IE     INTEGER   PARTICLE NUMBER
C     LIST   INTEGER   LIST OF FOUND PARTICLES
C     NLIST  INTEGER   NUMBER OF PARTICLES FOUND
C     NBLKSZ INTEGER   BLOCK SIZE OF IUP AND ILO BLOCKS
C     NSPC   INTEGER   NUMBER OF SPACIAL COORD. (NUMBER OF DIMENSIONS)
C
C-----------------------------------------------------------------------
C     
      DIMENSION
     *  IND(NUMPTS,NSPC),IRNK2(NUMPTS,NSPC,*),
     *  IUP(NBLKSZ,NSPC),INDX(NUMPTS),
     *  ILO(NBLKSZ,NSPC), LIST(NUMPTS)
C
C BUILD A LIST OF POINTS THAT ARE CLOSE TO SURFACE IE
      J = IE
      NLIST = 0
      IF( NSPC .EQ. 1)THEN
C============================== o n e   -  d ======================
        NUM1 = IUP(J,1) - ILO(J,1) + 1
        ILOW = ILO(J,1)
        IUPR = IUP(J,1)      
        DO 101 I1 = ILOW, IUPR
          NLIST = NLIST +1
          LIST(NLIST) = IND(I1,1)
 101    CONTINUE
C
      ELSE IF( NSPC .EQ. 2 )THEN
C============================== t w o   -  d ======================
        NUM1 = IUP(J,1) - ILO(J,1) + 1
        NUM2 = IUP(J,2) - ILO(J,2) + 1
C DO WE HAVE A LIST ? 
        IF( NUM2.LE.0 .OR. NUM1.LE.0 ) RETURN
C SELECT WHICH LIST IS THE SMALLEST
        IF( NUM1 .LE. NUM2 )THEN
          IXYZ = 1
          IY = 2
          NUM = NUM1
        ELSE
          IXYZ = 2
          IY = 1
          NUM = NUM2
        ENDIF
C
        ILOW = ILO(J,IXYZ)
        IUPR = IUP(J,IXYZ)      
C FIRST TEST
        IF( NUM .GT. 64 ) THEN
          DO 201 I1 = ILOW, IUPR
            IF( IRNK2(I1,IXYZ,1) .GE. ILO(J,IY) .AND. 
     *        IRNK2(I1,IXYZ,1) .LE. IUP(J,IY) )THEN
              NLIST = NLIST +1
              LIST(NLIST) = IND(I1,IXYZ)
            ENDIF
 201      CONTINUE
        ELSE
CDIR$ SHORT LOOP
          DO 202 I1 = ILOW, IUPR
            IF( IRNK2(I1,IXYZ,1) .GE. ILO(J,IY) .AND. 
     *        IRNK2(I1,IXYZ,1) .LE. IUP(J,IY) )THEN
              NLIST = NLIST +1
              LIST(NLIST) = IND(I1,IXYZ)
            ENDIF
 202      CONTINUE
        ENDIF
C
      ELSE IF( NSPC .EQ. 3 )THEN
C============================== t h r e e   -  d ======================
        NUM1 = IUP(J,1) - ILO(J,1) + 1
        NUM2 = IUP(J,2) - ILO(J,2) + 1
        NUM3 = IUP(J,3) - ILO(J,3) + 1
C DO WE HAVE A LIST ? 
        IF( NUM3.LE.0 .OR. NUM2.LE.0 .OR. NUM1.LE.0 ) RETURN
C SELECT WHICH LIST IS THE SMALLEST
        IF( NUM1 .LE. NUM2 .AND. NUM1 .LE. NUM3 )THEN
          IXYZ = 1
          IY = 2
          IZ = 3
          NUM = NUM1
        ELSEIF( NUM2 .LE. NUM1 .AND. NUM2 .LE. NUM3 )THEN
          IXYZ = 2
          IY = 1
          IZ = 3
          NUM = NUM2
        ELSE
          IXYZ = 3
          IY = 1
          IZ = 2
          NUM = NUM3
        ENDIF
C
        ILOW = ILO(J,IXYZ)
        IUPR = IUP(J,IXYZ)      
        IF (ILOW.EQ.0) THEN
          NLIST = 0
          RETURN
        ENDIF
        ILP = 0
C FIRST TEST
        IF( NUM .GT. 64 ) THEN
          DO 301 I1 = ILOW, IUPR
            IF( IRNK2(I1,IXYZ,1) .GE. ILO(J,IY) .AND. 
     *        IRNK2(I1,IXYZ,1) .LE. IUP(J,IY) )THEN
              ILP = ILP +1
              INDX(ILP) = I1
            ENDIF
 301      CONTINUE
        ELSE
CDIR$ SHORT LOOP
          DO 302 I1 = ILOW, IUPR
            IF( IRNK2(I1,IXYZ,1) .GE. ILO(J,IY) .AND. 
     *        IRNK2(I1,IXYZ,1) .LE. IUP(J,IY) )THEN
              ILP = ILP +1
              INDX(ILP) = I1
            ENDIF
 302      CONTINUE
        ENDIF
C SECOND TEST
        IF( ILP .GT. 64 ) THEN
CDIR$ IVDEP
          DO 311 I1 = 1, ILP
            IF( IRNK2(INDX(I1),IXYZ,2) .GE. ILO(J,IZ) .AND.
     *        IRNK2(INDX(I1),IXYZ,2) .LE. IUP(J,IZ) )THEN
              NLIST = NLIST + 1
              LIST(NLIST) = IND(INDX(I1),IXYZ)
            ENDIF
 311      CONTINUE
        ELSE
CDIR$ SHORT LOOP
CDIR$ IVDEP
          DO 313 I1 = 1, ILP
            IF( IRNK2(INDX(I1),IXYZ,2) .GE. ILO(J,IZ) .AND.
     *        IRNK2(INDX(I1),IXYZ,2) .LE. IUP(J,IZ) )THEN
              NLIST = NLIST + 1
              LIST(NLIST) = IND(INDX(I1),IXYZ)
            ENDIF
 313      CONTINUE
        ENDIF
      ENDIF
C
      RETURN
      END
C
