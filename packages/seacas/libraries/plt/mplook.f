C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION MPLOOK(VX,VY,VZ,PX,PY,PZ,TWIST)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPLOOK')
      PARAMETER (DPR=57.2958)

      MPLOOK = .FALSE.
      DENTHE = SQRT((PX-VX)**2+ (PZ-VZ)**2)
      IF (DENTHE.EQ.0.) THEN
         SINTHE = 0.
         COSTHE = 1.

      ELSE
         SINTHE = (PX-VX)/DENTHE
         COSTHE = (VZ-PZ)/DENTHE
      END IF

      DENPHI = SQRT((PX-VX)**2+ (PY-VY)**2+ (PZ-VZ)**2)
      IF (DENPHI.EQ.0.) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the same eye position as the reference positio
     *n',2)
         RETURN

      END IF

      MPLOOK = .TRUE.
      SINPHI = (VY-PY)/DENPHI
      COSPHI = DENTHE/DENPHI
      CALL LDTRAN(-VX,-VY,-VZ,TMAT1)
      CALL LDROTA('y',COSTHE,SINTHE,TMAT2)
      CALL MXMULT(4,TMAT1,TMAT2,TMAT3)
      CALL LDROTA('x',COSPHI,SINPHI,TMAT1)
      CALL MXMULT(4,TMAT3,TMAT1,TMAT2)
      ANG = -TWIST/DPR
      CALL LDROTA('z',COS(ANG),SIN(ANG),TMAT1)
      CALL MXMULT(4,TMAT2,TMAT1,VIEW)
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      PEYE(1) = VX
      PEYE(2) = VY
      PEYE(3) = VZ
      PLOOK(1) = PX
      PLOOK(2) = PY
      PLOOK(3) = PZ
      ETWIST = TWIST
      RETURN

      END
