C-----------------------------------------------------------------------
C            ------------<   3 . h   >------------
C-----------------------------------------------------------------------
C  Arrays used for the solution of annular poisson problem
C  written so operator P could be arbitrary second order operator, of
C  bandwidth <= l_max
C-----------------------------------------------------------------------
        double precision     HP(0:m_max     ,  0:1       , 0:n_max)
        double precision     CP(0:m_max     ,  1:2       , 0:n_max)
        double precision     VP(1:2         ,  0:1       , 0:n_max)
        double precision     VPa(1:2         ,  0:1       , 0:n_max)
        double precision    CHP(1:2         ,  0:1       , 0:n_max)
        double precision     HP2(0:m_max     ,  0:1       , 0:n_max)
        double precision     CP2(0:m_max     ,  1:2       , 0:n_max)
        double precision     VPb(1:2         ,  0:1       , 0:n_max)
        double precision     VP2(1:2         ,  0:1       , 0:n_max)
        double precision    CHP2(1:2         ,  0:1       , 0:n_max)
C-----------------------------------------------------------------------
	common/aux/HP,CP,VP,CHP, HP2, CP2, VPa, VPb ,CHP2, VP2
C-----------------------------------------------------------------------
c       double precision     O3(0:m_max+buff,  1:3       , 0:n_max)
        double precision     O5(0:m_max+buff,  1:5       , 0:n_max)
        double precision     O25(0:m_max+buff,  1:5       , 0:n_max)
c       double precision     O7(0:m_max+buff,  1:7       , 0:n_max)
c       double precision     O9(0:m_max+buff,  1:9       , 0:n_max)
c       double precision     O13(0:m_max+buff, 1:13      , 0:n_max)
c       double precision     P3(0:m_max+buff,  1:3       , 0:n_max)
        double precision     P5(0:m_max+buff,  1:5       , 0:n_max)
        double precision     P25(0:m_max+buff,  1:5       , 0:n_max)
c       double precision     P7(0:m_max+buff,  1:7       , 0:n_max)
c       double precision     P9(0:m_max+buff,  1:9       , 0:n_max)
c       double precision     P13(0:m_max+buff, 1:13     , 0:n_max)
        double precision     Q5(0:m_max+buff,  1:5       , 0:1)
        double precision     Q25(0:m_max+buff,  1:5       , 0:1)
c       double precision     Q9(0:m_max+buff,  1:9       , 0:1)
C-----------------------------------------------------------------------
c       common/ops3/O3,P3
        common/ops5/O5,P5,Q5, Q25, O25, P25
c       common/ops7/O7,P7,Q9
c       common/ops9/O9,P9
c       common/ops13/O13,P13
C-----------------------------------------------------------------------
        double precision     CD(0:m_max     ,  1:2       , 0:n_max)
        double precision     CN(0:m_max     ,  1:2       , 0:n_max)
        double precision     NS(0:m_max     ,  1:2       , 0:n_max)
C-----------------------------------------------------------------------
	common/constraints/CD,CN,NS
C-----------------------------------------------------------------------
