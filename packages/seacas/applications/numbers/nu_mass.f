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

C $Id: mass.f,v 1.7 2007/03/21 20:12:37 gdsjaar Exp $
C $Log: mass.f,v $
C Revision 1.7  2007/03/21 20:12:37  gdsjaar
C Several commands which can work on the deformed geometry were only
C checking whether the file was an exodus file (had timesteps) when
C requesting deformed coordinates.  Changed to also check whether the
C file had valid displacments also.
C
C Revision 1.6  2000/07/06 18:07:42  gdsjaar
C Fix assumption that variables are saved between subroutine calls
C
C Revision 1.5  1999/04/20 22:33:59  gdsjaar
C Two arrays in different areas of the code used the same name for an
C array (JACOB). If the code was run in the wrong order, there would be
C a supes error when the array was reserved for the second time.
C
C Revision 1.4  1999/02/16 21:38:00  gdsjaar
C Converted to read exodusII database format.  Somewhat tested, not
C ready for production yet.
C
C Revision 1.3  1998/03/22 05:34:37  gdsjaar
C General cleanp of unused variables. Reordered DATA statements in
C command.f so would compile with f2c.
C
C Revision 1.2  1993/07/21 22:36:54  gdsjaar
C Removed unused variable--error
C
c Revision 1.1.1.1  1991/02/21  15:44:20  gdsjaar
c NUMBERS: Greg Sjaardema, initial Unix release
c
c Revision 1.1  1991/02/21  15:44:19  gdsjaar
c Initial revision
c
      SUBROUTINE MASSPR (A, TIME, ITMSEL, DENS, MAT, DISP,
     *   NQUAD, LABEL)
C
      DIMENSION A(*), TIME(*), DENS(*), MAT(6,*),
     *   DISP(NUMNP,*)
      LOGICAL ITMSEL(*), ISABRT
      CHARACTER*16  LABEL(32)
      include 'nu_ptim.blk'
      include 'nu_numg.blk'
      include 'nu_mass.blk'
      include 'nu_logs.blk'
C
      DIMENSION XI2(2,4), XI3(3,8)
      LOGICAL FIRST, HAVDEN
      DATA FIRST / .TRUE. /
      DATA XI2/ -1.,-1.,  1.,-1.,  1.,1.,  -1.,1./
      DATA XI3/ 1.,-1.,-1.,  -1.,-1.,-1.,  -1.,-1.,1.,  1.,-1.,1.,
     *   1.,1.,-1.,   -1.,1.,-1.,   -1.,1.,1.,   1.,1.,1./
C
      save

      IF (FIRST) THEN
         FIRST = .FALSE.
         CALL MDRSRV ('MASS'  , IS, NELBLK)
         CALL MDRSRV ('VOLUME', IV, NELBLK)
         CALL MDRSRV ('CENTER', IC, 3)
         CALL MDRSRV ('INERTA', IZ, 6)
         NNODES = 2**NDIM
         NQMAX  = 2**NDIM
         CALL MDRSRV ('XXX'   , IXXX,  (NDIM+1)*NNODES*NQMAX)
         CALL MDRSRV ('XG'    , IXG,   NDIM*NQMAX)
         CALL MDRSRV ('XINI'  , IXINI, NDIM)
C ... 'JACOB' conflicts with jacob in command.f, renamed to jacob1
         CALL MDRSRV ('JACOB1', IAJ,   NDIM*NDIM)
         CALL MDRSRV ('VOL'   , IVM,   4*NELBLK)
         CALL MDRSRV ('IEL'   , IEM,   4*NELBLK)
         CALL MDSTAT (NERRS, NUSED)
         IF (NERRS .GT. 0) THEN
            CALL MEMERR
            STOP
         END IF
      END IF
C
      HAVDEN = .FALSE.
      DO 20 I=1,NELBLK
         IF (DENS(I) .NE. 0.0) HAVDEN = .TRUE.
   20 CONTINUE

      IF (.NOT. HAVDEN) CALL GETDEN (MAT, DENS, NELBLK, LABEL)
C
      IF (EXODUS .AND. ISDIS) THEN
         CALL GETDSP (A(IR), DISP, NDIM, NUMNP, TIME, ITMSEL, 
     *      'R', ISTAT)
         IF (ISTAT .NE. 0) GO TO 40

   30    CONTINUE
         IF (ISABRT()) RETURN
         CALL GETDSP (A(IR), DISP, NDIM, NUMNP, TIME, ITMSEL, 
     *      'A', ISTAT)
         IF (ISTAT .NE. 0) GO TO 40
         IF (NDIM .EQ. 2) THEN
            CALL CGCAL2 (DISP,A(IX),MAT,A(IS),VOL,A(ID),
     *         A(IV),A(IC),A(IZ),A(IXXX),A(IXG),XI2,
     *         A(IXINI),A(IAJ),NNODES,NDIM,NQUAD,
     *         A(IVM),A(IEM),NELBLK,AXI,NUMNP)
         ELSE IF (NDIM .EQ. 3) THEN
            CALL CGCAL3 (DISP,A(IX),MAT,A(IS),VOL,A(ID),
     *         A(IV),A(IC),A(IZ),A(IXXX),A(IXG),XI3,
     *         A(IXINI),A(IAJ),NNODES,NDIM,NQUAD,
     *         A(IVM),A(IEM),NELBLK,NUMNP)
         END IF
C
         CALL OUTPUT (A(IS), A(ID), A(IV), A(IC), A(IZ), MAT,
     *      NDIM,NELBLK, VOL, A(IVM), A(IEM), 
     *      NQUAD, LABEL, AXI, TREAD)
C
         GO TO 30
   40    CONTINUE
      ELSE
         IF (NDIM .EQ. 2) THEN
            CALL CGCAL2 (A(IR),A(IX),MAT,A(IS),VOL,A(ID),
     *         A(IV),A(IC),A(IZ),A(IXXX),A(IXG),XI2,
     *         A(IXINI),A(IAJ),NNODES,NDIM,NQUAD,
     *         A(IVM),A(IEM),NELBLK,AXI,NUMNP)
         ELSE 
            CALL CGCAL3 (A(IR),A(IX),MAT,A(IS),VOL,A(ID),
     *         A(IV),A(IC),A(IZ),A(IXXX),A(IXG),XI3,
     *         A(IXINI),A(IAJ),NNODES,NDIM,NQUAD,
     *         A(IVM),A(IEM),NELBLK,NUMNP)
         END IF
         CALL OUTPUT (A(IS), A(ID), A(IV), A(IC), A(IZ), MAT,
     *      NDIM,NELBLK, VOL, A(IVM), A(IEM), NQUAD, LABEL, 
     *      AXI, TREAD)
C
      END IF
      RETURN
      END
