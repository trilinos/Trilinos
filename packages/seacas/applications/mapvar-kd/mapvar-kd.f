C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK,MAPVAR
      PROGRAM MAPVAR

C     ******************************************************************

C                                --MAPVAR--
C                  A PROGRAM TO MAP FINITE ELEMENT RESULTS
C                 FROM ONE EXODUS-II RESTART FILE TO ANOTHER
C                 EXODUS-II RESTART FILE TO SUPPORT REMESHING

C                             GERALD W. WELLMAN
C                        SANDIA NATIONAL LABORATORIES
C                          ALBUQUERQUE, NEW MEXICO

C     MAPVAR IS BASED ON MERLIN II,A FINITE ELEMENT INTERPOLATION
C     PROGRAM BY DAVID K. GARTLING.

C     THE MERLIN PROGRAM IS DESIGNED TO TRANSFER DATA BETWEEN TWO- AND
C     THREE-DIMENSIONAL FINITE ELEMENT MESHES. GIVEN A FINITE ELEMENT
C     MESH, (MESH-A), AND A SOLUTION FIELD ON THAT MESH, MERLIN
C     INTERPOLATES THE GIVEN SOLUTION ONTO A SECOND FINITE ELEMENT MESH,
C     (MESH-B). THE INTERPOLATION PROCEDURE IS BASED ON USE OF THE
C     ELEMENT SHAPE FUNCTIONS IN MESH-A.

C     MAPVAR IS DESIGNED TO SERVE THE SAME GENERAL PURPOSE AS MERLIN
C     HOWEVER, MAPVAR IS DESIGNED TO PROVIDE THE TRANSLATION IN TERMS
C     OF EXODUS-II-V2 RESTART FILES. MAPVAR ALSO TRANSLATES ELEMENT
C     VARIABLES AS WELL AS NODAL VARIABLES AND READS/WRITES GLOBAL
C     VARIABLES. MAPVAR IS CURRENTLY SET UP FOR THREE ELEMENT TYPES,
C     2-D QUADS, 3-D HEXES, AND 3-D QUAD SHELLS. ALL THE ELEMENTS
C     ORIGINALLY SUPPORTED BY MERLIN CAN BE INCLUDED IN MAPVAR GIVEN
C     THE DESIRE AND RESOURCES. THE SEARCH ENGINE OF MAPVAR HAS BEEN
C     CHANGED TO A BINARY SEARCH FROM THE BIN OR BUCKET SEARCH OF MERLIN.

C     THE INTENT OF MAPVAR IS TO CREATE A RESTART FILE THAT WILL ALLOW
C     A FINITE ELEMENT SOLUTION TO PROCEED WITH A DIFFERENT MESH THAN
C     THE MESH WITH WHICH THE SOLUTION WAS STARTED. THUS, THERE IS AN
C     INHERENT ASSUMPTION THAT MESH-A IS AN EXODUS RESTART FILE.

C     NODAL VARIABLE TRANSFER IS STRAIGHT FORWARD. THE MESH-B NODE IS
C     FOUND INSIDE THE APPROPRIATE MESH-A ELEMENT. THE ELEMENT SHAPE
C     FUNCTIONS ARE USED TO INTERPOLATE RESULTS ONTO THE MESH-B NODE.
C     ELEMENT VARIABLE TRANSFER TAKES TWO STEPS. FIRST THE ELEMENT
C     VARIABLE IS "SCATTERED" TO THE NODES. THEN THE MESH-B ELEMENT
C     CENTROID IS FOUND INSIDE THE MESH-A ELEMENT. TRANSFER TO THE
C     MESH-B ELEMENT CENTROID TAKES PLACE USING THE ELEMENT SHAPE
C     FUNCTIONS AND THE "SCATTERED" ELEMENT DATA (NOW NODAL DATA).
C     ELEMENT DATA IS "SCATTERED" TO THE NODES IN TWO SCHEMES. THE
C     FIRST IS SIMPLE AVERAGING. THIS IS ROBUST AND FAST BUT CAN LEAD
C     TO SIGNIFICANT DISSIPATION IN THE CASE OF GRADIENTS APPROACHING
C     A FREE SURFACE. THE SECOND METHOD EMPLOYS A CONSTRAINED LEAST
C     SQUARES FITTING PROCEDURE TO "SCATTER" ELEMENT RESULTS TO THE
C     NODES. THE CONSTRAINTS ARE BASED ON THE REQUIREMENTS OF THE
C     CONSTITUTIVE MODEL. FINALLY, A DIRECT TRANSFER HAS BEEN IMPLEMENTED
C     BUT IS NOT RECOMMENDED.

C     Special considerations:
C       ELMASS - translated to nodal density, interpolated, translated
C                back to ELMASS
C       ROTATIONS - magnitude of matrix must be 1, translate, compute
C                   magnitude, divide each component by magnitude
C       EQPS - constrained to be .gt. 0.
C       TEARING - constrained to be .gt. 0.
C       DECAY - constrained to be .lt. 1.

      include 'exodusII.inc'
      CHARACTER*(MXSTLN) TYP
      CHARACTER*8   MEMDBG

      include 'aexds1.blk'
      include 'amesh.blk'
      include 'bmesh.blk'
      include 'contrl.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'header.blk'
      include 'rundat.blk'
      include 'steps.blk'
      include 'schdat.blk'
      include 'tapes.blk'
      include 'varnpt.blk'
      include 'varept.blk'

      DIMENSION A(1),IA(1)
      EQUIVALENCE (A(1),IA(1))

C     ******************************************************************

C     MAIN PROGRAM FOR MAPVAR
C     PROGRAM EXECUTION IS DIRECTED FROM THIS ROUTINE

C     ******************************************************************

C     NOTE : ALL ELEMENT DATA,NODAL POINT DATA, SOLUTION FIELDS, WORK
C            SPACE, ETC. ARE STORED IN THE ARRAY "A". THE MEMORY
C            MANAGER REQUESTS THE NEEDED STORAGE SPACE IN "A" DURING
C            EXECUTION. THE POINTERS NA1,NA2,....NAM PARTITION THE
C            ALLOCATED STORAGE INTO SEPARATE ARRAYS. SEE SAND86-0911,
C            "SUPES", FOR DETAILS OF THE MEMORY MANAGER.

C     ******************************************************************

C disable netcdf warning messages

      CALL EXOPTS(EXVRBS,IERR)

C open all disk files

      call debug('MVOPNFIL')
      CALL MVOPNFIL

C get info for QA records

      CALL VERSION(QAINFO)

C initialize memory manager

      CALL MDINIT (A)

C ... If EXT99 Environment variable set, turn on supes memory debugging
C     The numeric value of the variable is used as the unit to write
C     debug information to.
      CALL EXNAME (-99, MEMDBG, L)
      IF (L .GE. 1) THEN
        READ(MEMDBG(:L), '(I8)', ERR=20) IUNIT
        CALL MDDEBG(IUNIT)
      END IF
 20   CONTINUE

C ******************************************************************
C read input parameters needed to control the solution from
C      the screen (time step(s), bins)
C ******************************************************************

      CALL EXINQ(NTP2EX,EXTIMS,NTIMES,RDUM,CDUM,IERR)
      CALL EXGINI (NTP2EX,HED,NDIMA,NODESA,NUMELA,NBLKSA,
     &             NUMNPS,NUMESS,IERR)
      CALL EXGINI (NTP3EX,HED,NDIMB,NODESB,NUMELB,NBLKSB,
     &             NUMNPS,NUMESS,IERR)

C A(NT1)      =   TIMES(1:NTIMES) - times on Mesh-A database
C IA(NAEB)    =   IDA(1:NBLKSA) - Donor mesh element block I.D.'s
C IA(NBEB)    =   IDB(1:NBLKSA) - Recipient mesh element block I.D.'s
C IA(NMAP)    =   MP(1:3,1:MBLK) - Donor to recipient mesh map

      MBLK = NBLKSA * NBLKSB
      CALL MDRSRV ('TIMES', NT1,   NTIMES)
      CALL MDRSRV ('IDA',   NAEB,  NBLKSA)
      CALL MDRSRV ('IDB',   NBEB,  NBLKSB)
      CALL MDRSRV ('MP',    NMAP,  MBLK*3)

C reserve space for storing the search box size for each map operation
      CALL MDRSRV ('MPSEA', NMAPS, MBLK)

      CALL MDSTAT (MNERRS, MNUSED)
      IF (MNERRS .NE. 0) THEN
         CALL MDEROR(NOUT)
         CALL ERROR('MAPVAR',
     &              'MEMORY MANAGER ERROR',
     &              'JUST BEFORE CALL TO RDINPT',0,' ',0,' ',' ',1)
      END IF

      call debug('RDINPT')
      CALL RDINPT (A(NT1),IA(NAEB),IA(NBEB),IA(NMAP),A(NMAPS),IMP,MBLK)

C ******************************************************************
C INITIAL READ OF MESH-A
C ******************************************************************

C read sizing data for mesh A

      WRITE (NOUT, 270) fntp2(:lenstr(fntp2)), HED
      WRITE (NTPOUT, 270) fntp2(:lenstr(fntp2)), HED
      WRITE (NOUT, 280) NDIMA,NODESA,NUMELA,NBLKSA
      WRITE (NTPOUT, 280) NDIMA,NODESA,NUMELA,NBLKSA

C allocate storage for initial read of mesh A

C      A(NAX)     =    XA(1:NODESA) - Mesh-A X-coord
C      A(NAY)     =    YA(1:NODESA) - Mesh-A Y-coord
C      A(NAZ)     =    ZA(1:NODESA) - Mesh-A Z-coord
C      A(NADX)    =    DISXA(1:NODESA) - Mesh-A X-displ
C      A(NADY)    =    DISYA(1:NODESA) - Mesh-A Y-displ
C      A(NADZ)    =    DISZA(1:NODESA) - Mesh-A Z-displ

      CALL MDRSRV ('XA',     NAX,   NODESA)
      CALL MDRSRV ('YA',     NAY,   NODESA)
      IF (NDIMA.EQ.3) THEN
        CALL MDRSRV ('ZA',     NAZ,  NODESA)
      ELSE
        CALL MDRSRV ('ZA',     NAZ,  1)
      END IF

      if (idef .ne. 0) then
        CALL MDRSRV ('DISXA',  NADX,  NODESA)
        CALL MDRSRV ('DISYA',  NADY,  NODESA)
        IF (NDIMA.EQ.3) THEN
          CALL MDRSRV ('DISZA',  NADZ, NODESA)
        ELSE
          CALL MDRSRV ('DISZA',  NADZ, 1)
        END IF
      else
        CALL MDRSRV ('DISXA',  NADX,  1)
        CALL MDRSRV ('DISYA',  NADY,  1)
        CALL MDRSRV ('DISZA',  NADZ,  1)
      end if
      CALL MDSTAT (MNERRS, MNUSED)
      IF (MNERRS .NE. 0) THEN
         CALL MDEROR(NOUT)
         CALL ERROR('MAPVAR',
     &              'MEMORY MANAGER ERROR',
     &              'JUST BEFORE CALL TO RDA1',0,' ',0,' ',' ',1)
      END IF

C Copy "GENESIS" from mesh-B to mesh-C

      CALL EXCOPY(NTP3EX,NTP4EX,IERR)
      IF (IERR .NE. 0)
     &  CALL ERROR('MAPVAR',
     &             'ERROR WITH EXCOPY - GENESIS FILE COPY',
     &             ' ',0,' ',0,' ',' ',1)

c read mesh A (coords,displ,variable names,QA, INFO records),
c write mesh C, (variable names, QA, INFO records)

      call debug('RDA1')
      CALL RDA1 (A(NAX),A(NAY),A(NAZ),A(NADX),A(NADY),A(NADZ))

C ... Mapvar is buggy if mesh contains:
C     * nodal variable(s)
C     * multiple element blocks
C     * multiple timesteps

C     The interpolated mesh will have invalid values for most of the
C     nodes except for those connected to the last element block.
C     Since this is not what the user wants. We check this situation
C     here and then refuse to run...
      if (nvarnp .gt. 0) then
        if (ntimes .gt. 1 .and. nblksa .gt. 1) then
        write (*,500)
 500      format(10x,/,
     *      'FATAL ERROR: Mapvar is buggy and cannot correctly'
     *      ' handle interpolating a mesh with nodal variables,',
     *      /,14x,'multiple timesteps and multiple element blocks.',
     *      /,14x,'See https://github.com/gsjaardema/seacas/packages'
     *      '/seacas/MAPVAR.md')
        stop('INTERNAL ERROR')
      end if

C ... Warn if multiple blocks...
        if (ntimes .gt. 1 .and. nblksa .gt. 1) then
          write (*,600)
 600      format(10x,/,
     *      'WARNING: Mapvar is buggy and might not correctly'
     *      ' handle interpolating a mesh with nodal variables',
     *      /,9x,'and multiple element blocks. Shared nodes',
     *      ' might be incorrect.',
     *      /,9x,'See https://github.com/gsjaardema/seacas/packages'
     *      '/seacas/MAPVAR.md')
        end if
      end if
      call mddel('DISXA')
      call mddel('DISYA')
      call mddel('DISZA')
      IF (MNERRS .NE. 0) THEN
         CALL MDEROR(NOUT)
         CALL ERROR('MAPVAR',
     &              'MEMORY MANAGER ERROR',
     &              'JUST AFTER CALL TO RDA1',0,' ',0,' ',' ',1)
      END IF

      IF (IACCU .EQ. 1)THEN
        IF (IXVEL .NE. 0 .AND. IYVEL .NE. 0 .AND.
     &      (IELMS .NE. 0 .OR. IDENS .NE. 0))THEN

C velocities and mass are available, compute momenta and k.e.
C 1st set up storage for vel's

C  A(NVXA)    = VELXA(1:NODESA) - X-velocity in mesh-A
C  A(NVYA)    = VELYA(1:NODESA) - Y-velocity in mesh-A
C  A(NVZA)    = VELZA(1:NODESA) - Z-velocity in mesh-A
C  A(NNMSA)   = RMSNA(1:NODESA) - nodal mass in mesh-A
C  A(NVXB)    = VELXB(1:NODESB) - X-velocity in mesh-B
C  A(NVYB)    = VELYB(1:NODESB) - Y-velocity in mesh-B
C  A(NVZB)    = VELZB(1:NODESB) - Z-velocity in mesh-B
C  A(NNMSB)   = RMSNB(1:NODESB) - nodal mass in mesh-B

          CALL MDRSRV ('VELXA',   NVXA,   NODESA)
          CALL MDRSRV ('VELYA',   NVYA,   NODESA)
          IF (NDIMA .EQ. 3)THEN
            CALL MDRSRV ('VELZA',   NVZA,   NODESA)
          ELSE
            CALL MDRSRV ('VELZA',   NVZA,   1)
          END IF
          CALL MDRSRV ('RNMSA',   NNMSA,  NODESA)
          CALL MDRSRV ('VELXB',   NVXB,   NODESB)
          CALL MDRSRV ('VELYB',   NVYB,   NODESB)
          IF (NDIMB .EQ. 3)THEN
            CALL MDRSRV ('VELZB',   NVZB,   NODESB)
          ELSE
            CALL MDRSRV ('VELZB',   NVZB,   1)
          END IF
          CALL MDRSRV ('RNMSB',   NNMSB,  NODESB)
          CALL MDSTAT (MNERRS, MNUSED)
          IF (MNERRS .NE. 0) THEN
            CALL ERROR ('MAPVAR',
     &                  'MEMORY MANAGER ERROR',
     &                  'CHECK ACCURACY - VELOCITY',
     &                  0,' ',0,' ',' ',1)
          END IF

C initialization quantities (6 each for now) that need to be
C summed over the element blocks if doing an accuracy check

          IF (ISTEP .EQ. -1)THEN

C need arrays (one slot for each time), else just scalars will do

C A(NTMXA)  =  TMXA(1:NTIMES)  - X-momentum all mesh-A blocks each time
C A(NTMYA)  =  TMYA(1:NTIMES)  - Y-momentum all mesh-A blocks each time
C A(NTMZA)  =  TMZA(1:NTIMES)  - Z-momentum all mesh-A blocks each time
C A(NTKEA)  =  TKEA(1:NTIMES)  - K.E. all mesh-A blocks each time
C A(NTPSQA) =  TPSQA(1:NTIMES) - Pressure squared mesh-A each time
C A(NTJ2A)  =  TJ2A(1:NTIMES)  - J2 mesh-A each time

C A(NTMXB)  =  TMXB(1:NTIMES)  - X-momentum all mesh-B blocks each time
C A(NTMYB)  =  TMYB(1:NTIMES)  - Y-momentum all mesh-B blocks each time
C A(NTMZB)  =  TMZB(1:NTIMES)  - Z-momentum all mesh-B blocks each time
C A(NTKEB)  =  TKEB(1:NTIMES)  - K.E. all mesh-B blocks each time
C A(NTPSQB) =  TPSQB(1:NTIMES) - Pressure squared mesh-B each time
C A(NTJ2B)  =  TJ2B(1:NTIMES)  - J2 mesh-B each time

            CALL MDRSRV('TMXA',  NTMXA,  NTIMES)
            CALL MDRSRV('TMYA',  NTMYA,  NTIMES)
            CALL MDRSRV('TMZA',  NTMZA,  NTIMES)
            CALL MDRSRV('TKEA',  NTKEA,  NTIMES)
            CALL MDRSRV('TPSQA', NTPSQA, NTIMES)
            CALL MDRSRV('TJ2A', NTJ2A, NTIMES)
            CALL MDRSRV('TMXB',  NTMXB,  NTIMES)
            CALL MDRSRV('TMYB',  NTMYB,  NTIMES)
            CALL MDRSRV('TMZB',  NTMZB,  NTIMES)
            CALL MDRSRV('TKEB',  NTKEB,  NTIMES)
            CALL MDRSRV('TPSQB', NTPSQB, NTIMES)
            CALL MDRSRV('TJ2B',  NTJ2B,  NTIMES)
            CALL MDSTAT (MNERRS, MNUSED)
            IF (MNERRS .NE. 0) THEN
              CALL ERROR ('MAPVAR',
     &                    'MEMORY MANAGER ERROR',
     &                    'CHECK ACCURACY - INITIALIZATION',
     &                    0,' ',0,' ',' ',1)
            END IF
          ELSE
            CALL MDRSRV('TMXA',  NTMXA,  1)
            CALL MDRSRV('TMYA',  NTMYA,  1)
            CALL MDRSRV('TMZA',  NTMZA,  1)
            CALL MDRSRV('TKEA',  NTKEA,  1)
            CALL MDRSRV('TPSQA', NTPSQA, 1)
            CALL MDRSRV('TJ2A',  NTJ2A,  1)
            CALL MDRSRV('TMXB',  NTMXB,  1)
            CALL MDRSRV('TMYB',  NTMYB,  1)
            CALL MDRSRV('TMZB',  NTMZB,  1)
            CALL MDRSRV('TKEB',  NTKEB,  1)
            CALL MDRSRV('TPSQB', NTPSQB, 1)
            CALL MDRSRV('TJ2B',  NTJ2B,  1)
            CALL MDSTAT (MNERRS, MNUSED)
            IF (MNERRS .NE. 0) THEN
              CALL ERROR ('MAPVAR',
     &                    'MEMORY MANAGER ERROR',
     &                    'CHECK ACCURACY - INITIALIZATION',
     &                    0,' ',0,' ',' ',1)
            END IF
          END IF
        END IF
      END IF

C ******************************************************************
C INITIAL READ OF MESH-B
C ******************************************************************

C read sizing data for mesh B

      WRITE (NOUT, 290) fntp3(:lenstr(fntp3)), HED
      WRITE (NTPOUT, 290) fntp3(:lenstr(fntp3)), HED
      WRITE (NOUT, 300) NDIMB,NODESB,NUMELB,NBLKSB
      WRITE (NTPOUT, 300) NDIMB,NODESB,NUMELB,NBLKSB

c quick initial check of compatibility mesh-A to mesh-B

      IF (NDIMB .NE. NDIMA) THEN
        CALL ERROR('MAPVAR',
     &             'MESH-B INCOMPATIBLE WITH MESH-A',
     &             'DIMENSION OF MESH-A',NDIMA,
     &             'DIMENSION OF MESH-B',NDIMB,
     &             ' ',' ',1)
      END IF

C allocate storage for mesh B read

C      A(NBX)      =    XB(1:NODESB) - Mesh-B X-coord
C      A(NBY)      =    YB(1:NODESB) - Mesh-B Y-coord
C      A(NBZ)      =    ZB(1:NODESB) - Mesh-B Z-coord

      CALL MDRSRV ('XB',     NBX, NODESB)
      CALL MDRSRV ('YB',     NBY, NODESB)
      IF (NDIMB.EQ.3) THEN
        CALL MDRSRV ('ZB',    NBZ, NODESB)
      ELSE
        CALL MDRSRV ('ZB',    NBZ, 1)
      END IF
      CALL MDSTAT (MNERRS, MNUSED)
      IF (MNERRS .NE. 0) THEN
         CALL MDEROR(NOUT)
         CALL ERROR('MAPVAR',
     &              'MEMORY MANAGER ERROR',
     &              'JUST BEFORE CALL TO RDB1',
     &              0,' ',0,' ',' ',1)
      END IF

c read coordinates for mesh B

      call debug('RDB1')
      CALL RDB1 (A(NBX),A(NBY),A(NBZ))

C *********************************************************
C START INTERPOLATION
C *********************************************************

C set up memory for arrays for nodal results and truth table
C these arrays stay around forever - they don't get deleted
C after each element block is processed like the arrays
C set up within the element block loop

C   A(NASOLN)    =    SOLNA(1:NODESA,1:NVARNP) - Mesh-A nodal data
C   A(NBSOLN)    =    SOLNB(1:NODESB,1:NVARNP) - Mesh-B interpolated
C                                                nodal data
C   IA(ITTA)     =    ITRTA(1:NVAREL,1:NBLKSA) - Mesh-A truth table
C   IA(ITTB)     =    ITRTB(1:NVAREL,1:NBLKSB) - Mesh-B truth table

C    A(NSN)      = SN(1:NODESB)  - storage for nodal vars in ininod
C    A(NSE)      = SE(1:NODESB)  - storage for element vars in ininod

      CALL MDRSRV ('SOLNA',  NASOLN,   NODESA*NVARNP)
      CALL MDRSRV ('SOLNB',  NBSOLN,   NODESB*NVARNP)
      CALL MDRSRV ('ITRTA',  ITTA,     NVAREL*NBLKSA)
      CALL MDRSRV ('ITRTB',  ITTB,     NVAREL*NBLKSB)
      CALL MDRSRV ('SN',     NSN,      NODESB)
      CALL MDRSRV ('SE',     NSE,      NODESB)

      CALL MDSTAT (MNERRS, MNUSED)
      IF (MNERRS .NE. 0) THEN
         CALL MDEROR(NOUT)
         CALL ERROR('MAPVAR',
     &              'MEMORY MANAGER ERROR',
     &              'JUST BEFORE INTERPOLATION LOOP',
     &              0,' ',0,' ',' ',1)
      END IF

      call inirea(nodesb*nvarnp, 0.0, a(nbsoln))

      call debug('TRUTBL')
      CALL TRUTBL(IA(NMAP),IMP,IA(NAEB),IA(NBEB),IA(ITTA),IA(ITTB))

C *********************************************************************

C START OF ELEMENT BLOCK-BY-ELEMENT BLOCK INTERPOLATION LOOP

C *********************************************************************

C store default values of search box tolerances per element type
      TOLSHC = TOLSHL
      TOLQAC = TOLQAD
      TOLHEC = TOLHEX
      TOLTEC = TOLTET

C     A(NAGV)    =    GVAR(1:NVARGP) Global variables

      CALL MDRSRV ('GVAR',   NAGV, NVARGP)

      call debug('WRTC')
      CALL WRTC(A(NBX),A(NBY),A(NBZ),A(NAGV),A(NBSOLN))

      DO 50 IM = 1, IMP
        call getval(ia(nmap), im, idblka, idblkb, ischem)
        TOLSEA = A(NMAPS+IM-1)

        do 15 i=1, nblksa
          if (idblka .eq. ia(naeb-1+i)) then
            iblka = i
            go to 16
          endif
 15     continue
 16     continue

        do 25 i=1, nblksb
          if (idblkb .eq. ia(nbeb-1+i)) then
            IBLKB = i
            go to 26
          endif
 25     continue
 26     continue

C set up controls for many to 1 map
C if first time recipient mesh element block called, insub = 1
C else insub = 2
C if last time recipient mesh element block called, icompl = 1
C else icompl = 0

        INSUB  = 1
        ICOMPL = 1
        IF (IM .GT. 1)THEN
          IDBBM1 = IA(NMAP-5)
          IF (IDBLKB .EQ. IDBBM1)THEN
            INSUB = 2
          END IF
        END IF
        IF (IM .LT. IMP)THEN
          IDBBP1 = IA(NMAP+1)
          IF (IDBLKB .EQ. IDBBP1)THEN
            ICOMPL = 0
          END IF
        END IF

C     **********************************************************
C     ELEMENT BLOCK BY ELEMENT BLOCK INTERPOLATION
C     REQUIRED FOR ELEMENT DATA BUT ALSO USED FOR NODAL DATA
C     **********************************************************

        WRITE(NOUT,330)IM,IMP,IDBLKB
        WRITE(NOUT,320)IDBLKA
        WRITE(NTPOUT,330)IM,IMP,IDBLKB
        WRITE(NTPOUT,320)IDBLKA

        CALL EXGELB(NTP2EX,IDBLKA,TYP,NUMEBA,NELNDA,NATRIB,
     &              IERR)
        CALL EXGELB(NTP3EX,IDBLKB,TYP,NUMEBB,NELNDB,NATRIB,
     &              IERR)
        IF (NUMEBB .EQ. 0)THEN
          GO TO 50
        END IF

C set up arrays for element block-by-element block preliminaries
C these arrays will be deleted at the end of the loop

C IA(NACON)   =    ICONA(1:NELNDA,1:NUMEBA) - Mesh-A connectivity
C IA(NBCON)   =    ICONB(1:NELNDB,1:NUMEBB) - Mesh-B connectivity
C IA(NANDLST) =    NDLSTA(1:NODESA) - Mesh-A nodes in element block
C IA(NBNDLST) =    NDLSTB(1:NODESB) - Mesh-A nodes in element block
C  A(NASTAT)  =    STATUS(1:NUMEBA) - Mesh-A element status

        CALL MDRSRV ('ICONA',  NACON,   NELNDA*NUMEBA)
        CALL MDRSRV ('ICONB',  NBCON,   NELNDB*NUMEBB)
        CALL MDRSRV ('NDLSTA', NANDLST, NODESA)
        CALL MDRSRV ('NDLSTB', NBNDLST, NODESB)
        CALL MDRSRV ('STATUS', NASTAT,  NUMEBA)

        CALL MDSTAT (MNERRS, MNUSED)
        IF (MNERRS .NE. 0) THEN
          CALL MDEROR(NOUT)
          CALL ERROR('MAPVAR',
     &               'MEMORY MANAGER ERROR',
     &               'BLOCKS LOOP PRELIMINARIES',
     &               0,' ',0,' ',' ',1)
        END IF

c 2nd read of mesh A

        call debug('RDA2')
        CALL RDA2 (IDBLKA,IA(NACON),IA(NANDLST),A(NASTAT),
     &             MAXLN)

C Set the search box tolerance for the current mapping

        IF ( ITYPE .EQ. 13) THEN
C shell
          IF ( TOLSEA .GT. 0.0) THEN
            TOLSHL = TOLSEA
          ELSE
            TOLSHL = TOLSHC
          END IF
          WRITE( NOUT, 321) TOLSHL
          WRITE( NTPOUT, 321) TOLSHL
        ELSE IF ( ( ITYPE .EQ. 3) .OR.
     &            ( ITYPE .EQ. 4) .OR.
     &            ( ITYPE .EQ. 5)) THEN
C quad-4
          IF ( TOLSEA .GT. 0.0) THEN
            TOLQAD = TOLSEA
          ELSE
            TOLQAD = TOLQAC
          END IF
          WRITE( NOUT, 321) TOLQAD
          WRITE( NTPOUT, 321) TOLQAD
        ELSE IF ( ( ITYPE .EQ.  6) .OR.
     &            ( ITYPE .EQ. 10)) THEN
C hex-8 or tet-8
          IF ( TOLSEA .GT. 0.0) THEN
            TOLHEX = TOLSEA
            TOLTET = TOLSEA
          ELSE
            TOLHEX = TOLHEC
            TOLTET = TOLTEC
          END IF
          IF ( ITYPE .EQ. 6) THEN
            WRITE( NOUT, 321) TOLTET
            WRITE( NTPOUT, 321) TOLTET
          ELSE
            WRITE( NOUT, 321) TOLHEX
            WRITE( NTPOUT, 321) TOLHEX
          END IF
        ELSE
          CALL ERROR ('MAPVAR','INCORRECT ELEMENT TYPE',
     &                'ELEMENT TYPE =',ITYPE,
     &                'NOT YET IMPLEMENTED',0,' ',' ',1)
        END IF

c 2nd read of mesh-b

        call debug('RDB2')
        CALL RDB2(IDBLKB,IDBLKA,IA(NBCON),IA(NBNDLST))

C set up arrays for element block-by-element block processing
C these arrays will be deleted at the end of the loop

C IA(NS1)     =  ISRCHR(1:1(NISR),1:NUMNDB) Integer search results
C  A(NS2)     =  RSRCHR(1:6(NRSR),1:NUMNDB) Real search results
C IA(NS3)     =    LIST(1:NUMNDB)         Potential contacts
C  A(NS16)    =    XYZSRF(1:NODESA,1:3)   Coords defining element
C  A(NS17)    =    XYZPTS(1:NUMNDB,1:3)   Coords of points searched

C  A(NASOLE)  =    SOLEA(1:NUMEBA,1:NVAREL) - Mesh-A element data
C  A(NBSOLE)  =    SOLEB(1:NUMEBB,1:NVAREL) - Mesh-B interpolated
C                                              element data
C  A(NASOLEN) =    SOLENA(1:NODESA,1:NVAREL) - Mesh-A element data
C                                               "scattered" to nodes
C IA(NANELTN) =    NELTN(1:NODESA) - Mesh-A number of elements per
C                                     node Used for computing averages
C IA(NAINVLN) =    INVLNA(1:NODESA) - Mesh-Anumber of elements per
C                                      node in inverse connectivity
C IA(NAINVC)  =    INVCN(1:MAXLN,1:NODESA) Inverse connectivity
C IA(NACTR)   =    CNTRA(1:NUMEBA,1:3) - Mesh-A element centroids
        IDIM = MAX(NUMNDB, NUMEBB)
        CALL MDRSRV ('ISRCHR', NS1,  1*IDIM)
        CALL MDRSRV ('RSRCHR', NS2,  6*IDIM)
        CALL MDRSRV ('LIST',   NS3,  IDIM)
        CALL MDRSRV ('XYZSRF', NS16, NODESA*3)
        CALL MDRSRV ('XYZPTS', NS17, IDIM*3)

        CALL MDRSRV ('SOLEA',  NASOLE,  NUMEBA*NVAREL)
        CALL MDRSRV ('SOLEB',  NBSOLE,  NUMEBB*NVAREL)

        IF (ISCHEM .EQ. 0)THEN

          CALL MDRSRV ('SOLENA', NASOLEN, NODESA*NVAREL)
          CALL MDRSRV ('NELTN',  NANELTN, NODESA)
        ELSE IF (ISCHEM .EQ. 1)THEN

          CALL MDRSRV ('SOLENA', NASOLEN, NODESA*NVAREL)
          CALL MDRSRV ('INVLNA', NAINVLN, NODESA)
          CALL MDRSRV('INVCN', NAINVC, MAXLN*NODESA)
          CALL MDRSRV ('CNTRA',  NACTR,   NUMEBA*3)
        ELSE IF (ISCHEM .EQ. 2)THEN
          CONTINUE
        ELSE IF (ISCHEM .EQ. 3)THEN

          CALL MDRSRV ('CNTRA',  NACTR,   NUMEBA*3)
          CALL MDRSRV ('INVLNA', NAINVLN, NODESA)
          CALL MDRSRV ('INVCN',  NAINVC,  MAXLN*NODESA)
          CALL MDRSRV ('ICHKEL', NICHKE,  NUMEBA)
          CALL MDRSRV ('SOLGRA', NSOLGR,  NDIMA*NUMEBA*NVAREL)
        ELSE

          CALL ERROR('MAPVAR',' ','ISCHEM',
     &               ischem,'INCORRECT ARGUMENT',0,' ',' ',1)
        END IF

        CALL MDSTAT (MNERRS, MNUSED)
        IF (MNERRS .NE. 0) THEN
          CALL MDEROR(NOUT)
          CALL ERROR('MAPVAR',
     &               'MEMORY MANAGER ERROR',
     &               'JUST AFTER START OF BLOCKS LOOP',
     &               0,' ',0,' ',' ',1)
        END IF

        IF (ITYPE .EQ. 13)THEN

C     **********************************************************
C     Path through code for shells
C     **********************************************************

          call debug('BLDSRF')
          CALL BLDSRF(A(NAX),A(NAY),A(NAZ),A(NS16))

          IF (NVARNP .GT. 0)THEN

            call debug('BLDPTN')
            CALL BLDPTN(A(NBX),A(NBY),A(NBZ),IA(NBNDLST),A(NS17))

            call debug('SRCHS-nodes')
            CALL SRCHS (NODESA,NUMEBA,IA(NACON),A(NS16),
     1       NUMNDB,A(NS17),TOLSHL,1,6,
     2       NISS,NRSS,IA(NS1),A(NS2),
     3       IA(NS3),IERR)

            call debug('SINTPN')
            CALL SINTPN(IA(NACON),A(NASOLN),IA(NS1),1,A(NS2),6,
     1                  A(NBSOLN),IA(NBNDLST),A(NBX),A(NBY),A(NBZ),
     2                  IDBLKB,A(NT1),INSUB,A(NSN))

          END IF
          IF (NVAREL .GT. 0)THEN

            call debug('BLDPTE')
            CALL BLDPTE(A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NS17))

            call debug('SRCHS-element centroids')
            CALL SRCHS (NODESA,NUMEBA,IA(NACON),A(NS16),
     1       NUMEBB,A(NS17),TOLSHL,1,6,
     2       NISS,NRSS,IA(NS1),A(NS2),
     3       IA(NS3),IERR)

c element centroid values to nodes by averaging

            IF (ISCHEM .EQ. 0)THEN

C Set up time steps

              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF

              DO 610 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF

                call debug('SETON0')
                CALL SETON0(IA(NACON),IA(NANELTN),A(NASOLE),
     &           A(NASOLEN),IDBLKA,A(NAX),A(NAY),A(NAZ),ISTP,
     &           IA(ITTB),iblkb)

                call debug('SINTPE')
                CALL SINTPE(IA(NACON),A(NASOLEN),IA(NS1),1,A(NS2),6,
     &                      A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     &                      IA(NBCON),IA(ITTB),IBLKB,A(NT1),A(NS17),
     &                      ISTP,IST,INSUB,ICOMPL,A(NSE))
  610         CONTINUE

            ELSE IF (ISCHEM .EQ. 1) THEN

              call debug('INVCON')
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))

C Set up time steps

              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF

              DO 620 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF

                call debug('SETON1')
                CALL SETON1(A(NACTR),A(NASOLE),A(NASOLEN),IDBLKA,
     &                    A(NAX),A(NAY),A(NAZ),IA(NACON),IA(NANDLST),
     &                    IA(NAINVLN),IA(NAINVC),MAXLN,ISTP,
     &                    IA(ITTB),iblkb)

                call debug('SINTPE')
                CALL SINTPE(IA(NACON),A(NASOLEN),IA(NS1),1,A(NS2),6,
     &                      A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     &                      IA(NBCON),IA(ITTB),IBLKB,A(NT1),A(NS17),
     &                      ISTP,IST,INSUB,ICOMPL,A(NSE))
  620         CONTINUE

            ELSE IF (ISCHEM .EQ. 2) THEN

c direct transfer, does not require scatter to nodes

              call debug('STRAN')
              CALL STRAN(IA(NS1),1,A(NASOLE),A(NBSOLE),
     &                    IDBLKA,IDBLKB,
     &                    IA(ITTB),IBLKB,A(NT1),A(NS17),
     &                    INSUB,ICOMPL,
     &                    A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NSE))

            ELSE IF (ISCHEM .EQ. 3)THEN

              call debug('INVCON')
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))

C Set up time steps

              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF

              DO 630 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF

                call debug('ELGRAD')
                CALL ELGRAD(A(NACTR),A(NAX),A(NAY),A(NAZ),
     &                      A(NASOLE),A(NSOLGR),IA(NICHKE),
     &                      IDBLKA,IA(NACON),IA(NAINVLN),IA(NAINVC),
     &                      MAXLN,ISTP,IA(ITTB),IBLKB)

                call debug('INTRP3')
                CALL INTRP3(A(NACTR),A(NS17),IA(NS1),
     &                      A(NBSOLE),A(NASOLE),A(NSOLGR),
     &                      IDBLKB,IA(ITTB),IBLKB,A(NT1),
     &                      ISTP,IST,INSUB,ICOMPL,
     &                      A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NSE))
 630          CONTINUE

            ELSE
              CALL ERROR('MAPVAR',' ','ISCHEM =',
     &                 ischem,'INCORRECT ARGUMENT',0,' ',' ',1)
            END IF

          END IF

C ITYPE = 3 - 4 node quad
C ITYPE = 4 - 8 node quad
C ITYPE = 5 - 9 node quad

        ELSE IF (ITYPE .EQ. 3 .OR. ITYPE .EQ. 4 .OR.
     &           ITYPE .EQ. 5) THEN

C     *****************************************************
C     PATH THROUGH CODE FOR CONTINUUM ELEMENTS
C              (QUAD-4)
C     *****************************************************

C find and store location of mesh-b nodes within mesh-a

          call debug('BLDSRF')
          CALL BLDSRF(A(NAX),A(NAY),A(NAZ),A(NS16))

          IF (NVARNP .GT. 0)THEN

            call debug('BLDPTN')
            CALL BLDPTN(A(NBX),A(NBY),A(NBZ),IA(NBNDLST),A(NS17))

            call debug('SRCHQ-nodes')
            CALL SRCHQ (NODESA,NUMEBA,IA(NACON),A(NS16),
     1       NUMNDB,A(NS17),TOLQAD,1,3,
     2       NISS,NRSS,IA(NS1),A(NS2),
     3       IA(NS3),IERR)
           DO 530 I = 1, NUMNDB
             IF (IA(NS1-1+I) .EQ. 0)THEN
               WRITE(NOUT,430)IA(NBNDLST-1+I),IDBLKB
               WRITE(NTPOUT,430)IA(NBNDLST-1+I),IDBLKB
             END IF
 530       CONTINUE

c   interpolate nodal variables

           call debug('INTRPN')
           CALL INTRPN(IA(NACON),A(NASOLN),IA(NS1),A(NS2),
     &       A(NBSOLN),IA(NBNDLST),A(NBX),A(NBY),A(NBZ),
     &       IDBLKB,A(NT1),INSUB,A(NSN))

c start element variable interpolation

c locate Mesh-B element centroid in Mesh-A

         END IF
          IF (NVAREL .GT. 0)THEN

            call debug('BLDPTE')
            CALL BLDPTE(A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NS17))

            call debug('SRCHQ-element centroids')
            CALL SRCHQ (NODESA,NUMEBA,IA(NACON),A(NS16),
     1       NUMEBB,A(NS17),TOLQAD,1,3,
     2       NISS,NRSS,IA(NS1),A(NS2),
     3       IA(NS3),IERR)

c element centroid variables averaged to nodes

            IF (ISCHEM .EQ. 0)THEN

C Set up time steps

              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF

              DO 640 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF

                call debug('ELTON0')
                CALL ELTON0(IA(NACON),IA(NANELTN),A(NASOLE),
     &           A(NASOLEN),IDBLKA,A(NAX),A(NAY),A(NAZ),ISTP,
     &           IA(ITTB),IBLKB)

c interpolate element vars

                call debug('INTRPE')
                CALL INTRPE(IA(NACON),A(NASOLEN),IA(NS1),A(NS2),
     1                    A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     2                    IA(NBCON),IA(ITTB),IBLKB, A(NT1),
     3                    A(NS17),ISTP,IST,INSUB,ICOMPL,A(NSE))
  640         CONTINUE

c element centroid variables linear least squares to nodes

            ELSE IF (ISCHEM .EQ. 1)THEN

              call debug('INVCON')
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))

C Set up time steps

              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF

              DO 650 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF

                call debug('ELTON1')
                CALL ELTON1(A(NACTR),A(NASOLE),A(NASOLEN),IDBLKA,
     &                    A(NAX),A(NAY),A(NAZ),IA(NACON),IA(NANDLST),
     &                    IA(NAINVLN),IA(NAINVC),MAXLN,ISTP,
     &                    IA(ITTB),IBLKB)

c interpolate element vars

                call debug('INTRPE')
                CALL INTRPE(IA(NACON),A(NASOLEN),IA(NS1),A(NS2),
     1                    A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     2                    IA(NBCON),IA(ITTB),IBLKB, A(NT1),
     3                    A(NS17),ISTP,IST,INSUB,ICOMPL,A(NSE))
  650         CONTINUE

            ELSE IF (ISCHEM .EQ. 2)THEN

c direct transfer from Mesh-A to Mesh-B

              call debug('TRANAB')
              CALL TRANAB(IA(NS1),A(NASOLE),A(NBSOLE),
     &                  IDBLKA,IDBLKB,
     &                  IA(ITTB),IBLKB,A(NT1),A(NS17),
     &                  INSUB,ICOMPL,
     &                  A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NSE))

            ELSE IF (ISCHEM .EQ. 3)THEN

              call debug('INVCON')
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))

C Set up time steps

              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF

              DO 660 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF

                call debug('ELGRAD')
                CALL ELGRAD(A(NACTR),A(NAX),A(NAY),A(NAZ),
     &                    A(NASOLE),A(NSOLGR),IA(NICHKE),
     &                    IDBLKA,IA(NACON),IA(NAINVLN),IA(NAINVC),
     &                    MAXLN,ISTP,IA(ITTB),IBLKB)

                call debug('INTRP3')
                CALL INTRP3(A(NACTR),A(NS17),IA(NS1),
     &                    A(NBSOLE),A(NASOLE),A(NSOLGR),
     &                    IDBLKB,IA(ITTB),IBLKB,A(NT1),
     &                    ISTP,IST,INSUB,ICOMPL,
     &                    A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NSE))
 660          CONTINUE
            ELSE
              CALL ERROR('MAPVAR',' ','ISCHEM =',
     &                 ischem,'INCORRECT ARGUMENT',0,' ',' ',1)
            END IF

          END IF

        ELSE IF (ITYPE .EQ. 10 .OR. ITYPE .EQ. 6) THEN

C     *****************************************************
C     PATH THROUGH CODE FOR 3-D CONTINUUM ELEMENTS
C              (HEX-8) OR (TET-8)
C     *****************************************************

C     FIND AND STORE LOCATION OF MESH-B NODES WITHIN MESH-A

          call debug('BLDSRF')
          CALL BLDSRF(A(NAX),A(NAY),A(NAZ),A(NS16))

          IF (NVARNP .GT. 0)THEN

            call debug('BLDPTN')
             CALL BLDPTN(A(NBX),A(NBY),A(NBZ),IA(NBNDLST),A(NS17))

            IF (ITYPE .EQ. 10)THEN
              call debug('SRCHH-nodes')
              CALL SRCHH (NODESA,NUMEBA,IA(NACON),A(NS16),
     1          NUMNDB,A(NS17),TOLHEX,1,3,
     2          NISS,NRSS,IA(NS1),A(NS2),IA(NS3),
     5          IERR)

            ELSEIF (ITYPE .EQ. 6)THEN
              call debug('SRCHT-nodes')
              CALL SRCHT (NODESA,NUMEBA,IA(NACON),A(NS16),
     1          NUMNDB,A(NS17),TOLTET,1,3,
     2          NISS,NRSS,IA(NS1),A(NS2),
     3          IA(NS3),IERR)
            END IF

c interpolate nodal variables

            call debug('INTRPN')
            CALL INTRPN(IA(NACON),A(NASOLN),IA(NS1),A(NS2),
     &                A(NBSOLN),IA(NBNDLST),A(NBX),A(NBY),A(NBZ),
     &                IDBLKB,A(NT1),INSUB,A(NSN))

c start element variable interpolation

c locate Mesh-B element centroid in Mesh-A

          END IF
          IF (NVAREL .GT. 0)THEN

            call debug('BLDPTE')
            CALL BLDPTE(A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NS17))

            IF (ITYPE .EQ. 10)THEN
              call debug('SRCHH-element centroids')
              CALL SRCHH (NODESA,NUMEBA,IA(NACON),A(NS16),
     1          NUMEBB,A(NS17),TOLHEX,1,3,
     2          NISS,NRSS,IA(NS1),A(NS2),
     3          IA(NS3),IERR)

            ELSEIF (ITYPE .EQ. 6)THEN
              call debug('SRCHT-element centroids')
              CALL SRCHT (NODESA,NUMEBA,IA(NACON),A(NS16),
     1          NUMEBB,A(NS17),TOLTET,1,3,
     2          NISS,NRSS,IA(NS1),A(NS2),
     3          IA(NS3),IERR)
            END IF

            IF (ISCHEM .EQ. 0)THEN

C Set up time steps

              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF

              DO 670 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF

                call debug('ELTON0')
                CALL ELTON0(IA(NACON),IA(NANELTN),A(NASOLE),
     &               A(NASOLEN),IDBLKA,A(NAX),A(NAY),A(NAZ),ISTP,
     &               IA(ITTB),IBLKB)

c interpolate element vars

                call debug('INTRPE')
                CALL INTRPE(IA(NACON),A(NASOLEN),IA(NS1),A(NS2),
     1                    A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     2                    IA(NBCON),IA(ITTB),IBLKB, A(NT1),
     3                    A(NS17),ISTP,IST,INSUB,ICOMPL,A(NSE))
  670         CONTINUE

            ELSE IF (ISCHEM .EQ. 1)THEN

              call debug('INVCON')
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))

C Set up time steps

              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF

              DO 680 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF

                call debug('ELTON1')
                CALL ELTON1(A(NACTR),A(NASOLE),A(NASOLEN),IDBLKA,
     &                    A(NAX),A(NAY),A(NAZ),IA(NACON),IA(NANDLST),
     &                    IA(NAINVLN),IA(NAINVC),MAXLN,ISTP,
     &                    IA(ITTB),IBLKB)

c interpolate element vars

                call debug('INTRPE')
                CALL INTRPE(IA(NACON),A(NASOLEN),IA(NS1),A(NS2),
     1                    A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     2                    IA(NBCON),IA(ITTB),IBLKB, A(NT1),
     3                    A(NS17),ISTP,IST,INSUB,ICOMPL,A(NSE))
  680         CONTINUE

            ELSE IF (ISCHEM .EQ. 2)THEN

c direct transfer from Mesh-A to Mesh-B

              call debug('TRANAB')
              CALL TRANAB(IA(NS1),A(NASOLE),A(NBSOLE),
     &                IDBLKA,IDBLKB,
     &                IA(ITTB),IBLKB,A(NT1),A(NS17),
     &                INSUB,ICOMPL,
     &                A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NSE))

            ELSE IF (ISCHEM .EQ. 3)THEN

              call debug('INVCON')
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))

C Set up time steps

              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF

              DO 690 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF

                call debug('ELGRAD')
                CALL ELGRAD(A(NACTR),A(NAX),A(NAY),A(NAZ),
     &                    A(NASOLE),A(NSOLGR),IA(NICHKE),
     &                    IDBLKA,IA(NACON),IA(NAINVLN),IA(NAINVC),
     &                    MAXLN,ISTP,IA(ITTB),IBLKB)

                call debug('INTRP3')
                CALL INTRP3(A(NACTR),A(NS17),IA(NS1),
     &                    A(NBSOLE),A(NASOLE),A(NSOLGR),
     &                    IDBLKB,IA(ITTB),IBLKB,A(NT1),
     &                    ISTP,IST,INSUB,ICOMPL,
     &                    A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NSE))
 690          CONTINUE

            ELSE
              CALL ERROR('MAPVAR',' ','ISCHEM =',
     &                 ischem,'INCORRECT ARGUMENT',0,' ',' ',1)
            END IF

          END IF
        ELSE
          CALL ERROR ('MAPVAR','INCORRECT ELEMENT TYPE',
     &              'ELEMENT TYPE =',ITYPE,
     &              'NOT YET IMPLEMENTED',0,' ',' ',1)

        END IF

        IF (IACCU .EQ. 1)THEN

C velocities and mass are available, compute momenta and k.e.
C 1st set up storage for mass

          IF (IELMS .NE. 0 .AND. IDENS .EQ. 0)THEN

C  A(NEMSA)   = EMSSA(1:NODESA) - element mass in mesh-A
C  A(NEMSB)   = EMSSB(1:NODESB) - element mass in mesh-B

            CALL MDRSRV ('EMSSA',   NEMSA,   NUMEBA)
            CALL MDRSRV ('DENSA',   NDENA,   1)
            CALL MDRSRV ('EMSSB',   NEMSB,   NUMEBB)
            CALL MDRSRV ('DENSB',   NDENB,   1)

            CALL MDSTAT (MNERRS, MNUSED)
            IF (MNERRS .NE. 0) THEN
              CALL ERROR ('MAPVAR',
     &                   'MEMORY MANAGER ERROR',
     &                   'CHECK ACCURACY - ELEMENT MASS',
     &                   0,' ',0,' ',' ',1)
            END IF
          ELSE IF(IDENS .NE. 0)THEN

C  A(NEMSA)   = EMSSA(1:NUMEBA) - element mass in mesh-A
C  A(NEMSB)   = EMSSB(1:NUMEBB) - element mass in mesh-B
C  A(NDENA)   = DENSA(1:NUMEBA) - element density in mesh-A
C  A(NDENB)   = DENSB(1:NUMEBB) - element density in mesh-B

            CALL MDRSRV ('EMSSA',   NEMSA,   NUMEBA)
            CALL MDRSRV ('DENSA',   NDENA,   NUMEBA)
            CALL MDRSRV ('EMSSB',   NEMSB,   NUMEBB)
            CALL MDRSRV ('DENSB',   NDENB,   NUMEBB)

            CALL MDSTAT (MNERRS, MNUSED)
            IF (MNERRS .NE. 0) THEN
              CALL ERROR ('MAPVAR',
     &                   'MEMORY MANAGER ERROR',
     &                   'CHECK ACCURACY - ELEMENT MASS',
     &                    0,' ',0,' ',' ',1)
            END IF
          END IF

          IF (NDIMA .EQ. 3)THEN

C  A(NSXXA)  = SIGXXA(1:NUMEBA) - XX component of stress tensor
C  A(NSYYA)  = SIGYYA(1:NUMEBA) - YY component of stress tensor
C  A(NSZZA)  = SIGZZA(1:NUMEBA) - ZZ component of stress tensor
C  A(NSXYA)  = SIGXYA(1:NUMEBA) - XY component of stress tensor
C  A(NSYZA)  = SIGYZA(1:NUMEBA) - YZ component of stress tensor
C  A(NSZXA)  = SIGZXA(1:NUMEBA) - ZX component of stress tensor
C  A(NSXXB)  = SIGXXB(1:NUMEBB) - XX component of stress tensor
C  A(NSYYB)  = SIGYYB(1:NUMEBB) - YY component of stress tensor
C  A(NSZZB)  = SIGZZB(1:NUMEBB) - ZZ component of stress tensor
C  A(NSXYB)  = SIGXYB(1:NUMEBB) - XY component of stress tensor
C  A(NSYZB)  = SIGYZB(1:NUMEBB) - YZ component of stress tensor
C  A(NSZXB)  = SIGZXB(1:NUMEBB) - ZX component of stress tensor

            CALL MDRSRV ('SIGXXA' , NSXXA,  NUMEBA)
            CALL MDRSRV ('SIGYYA' , NSYYA,  NUMEBA)
            CALL MDRSRV ('SIGZZA' , NSZZA,  NUMEBA)
            CALL MDRSRV ('SIGXYA' , NSXYA,  NUMEBA)
            CALL MDRSRV ('SIGYZA' , NSYZA,  NUMEBA)
            CALL MDRSRV ('SIGZXA' , NSZXA,  NUMEBA)
            CALL MDRSRV ('SIGXXB' , NSXXB,  NUMEBB)
            CALL MDRSRV ('SIGYYB' , NSYYB,  NUMEBB)
            CALL MDRSRV ('SIGZZB' , NSZZB,  NUMEBB)
            CALL MDRSRV ('SIGXYB' , NSXYB,  NUMEBB)
            CALL MDRSRV ('SIGYZB' , NSYZB,  NUMEBB)
            CALL MDRSRV ('SIGZXB' , NSZXB,  NUMEBB)

            CALL MDSTAT (MNERRS, MNUSED)
            IF (MNERRS .NE. 0) THEN
              CALL ERROR ('MAPVAR',
     &                   'MEMORY MANAGER ERROR',
     &                   'CHECK ACCURACY - ELEMENT MASS',
     &                    0,' ',0,' ',' ',1)
            END IF
          ELSE IF (NDIMA .EQ. 2)THEN

C  A(NSXXA)  = SIGXXA(1:NUMEBA) - XX component of stress tensor
C  A(NSYYA)  = SIGYYA(1:NUMEBA) - YY component of stress tensor
C  A(NSZZA)  = SIGZZA(1:NUMEBA) - ZZ component of stress tensor
C  A(NSXYA)  = SIGXYA(1:NUMEBA) - XY component of stress tensor
C  A(NSXXB)  = SIGXXB(1:NUMEBB) - XX component of stress tensor
C  A(NSYYB)  = SIGYYB(1:NUMEBB) - YY component of stress tensor
C  A(NSZZB)  = SIGZZB(1:NUMEBB) - ZZ component of stress tensor
C  A(NSXYB)  = SIGXYB(1:NUMEBB) - XY component of stress tensor

            CALL MDRSRV ('SIGXXA' , NSXXA,  NUMEBA)
            CALL MDRSRV ('SIGYYA' , NSYYA,  NUMEBA)
            CALL MDRSRV ('SIGZZA' , NSZZA,  NUMEBA)
            CALL MDRSRV ('SIGXYA' , NSXYA,  NUMEBA)
            CALL MDRSRV ('SIGYZA' , NSYZA,  1)
            CALL MDRSRV ('SIGZXA' , NSZXA,  1)
            CALL MDRSRV ('SIGXXB' , NSXXB,  NUMEBB)
            CALL MDRSRV ('SIGYYB' , NSYYB,  NUMEBB)
            CALL MDRSRV ('SIGZZB' , NSZZB,  NUMEBB)
            CALL MDRSRV ('SIGXYB' , NSXYB,  NUMEBB)
            CALL MDRSRV ('SIGYZB' , NSYZB,  1)
            CALL MDRSRV ('SIGZXB' , NSZXB,  1)

            CALL MDSTAT (MNERRS, MNUSED)
            IF (MNERRS .NE. 0) THEN
              CALL ERROR ('MAPVAR',
     &                   'MEMORY MANAGER ERROR',
     &                   'CHECK ACCURACY - ELEMENT MASS',
     &                    0,' ',0,' ',' ',1)
            END IF
          END IF

C Set up time steps

          IF (ISTEP .EQ. -1)THEN
            NTM = NTIMES
          ELSE
            NTM = 1
          END IF

          DO 710 IST = 1, NTM
            IF (ISTEP .EQ. -1)THEN
              ISTP = IST
            ELSE
              ISTP = ISTEP
            END IF

            CALL MKEI(IST,ISTP,A(NT1),IDBLKA,IA(NACON),IA(NANDLST),
     &               A(NAX),A(NAY),A(NAZ),A(NVXA),A(NVYA),A(NVZA),
     &               A(NEMSA),A(NDENA),A(NNMSA),
     &               A(NTMXA),A(NTMYA),A(NTMZA),A(NTKEA),A(NTPSQA),
     &                                                   A(NTJ2A),
     &               A(NSXXA),A(NSYYA),A(NSZZA),A(NSXYA),A(NSYZA),
     &                                                   A(NSZXA),
     &               IDBLKB,IA(NBCON),IA(NBNDLST),
     &               A(NBX),A(NBY),A(NBZ),A(NVXB),A(NVYB),A(NVZB),
     &               A(NEMSB),A(NDENB),A(NNMSB),
     &               A(NTMXB),A(NTMYB),A(NTMZB),A(NTKEB),A(NTPSQB),
     &                                          A(NTJ2B),ICOMPL,
     &               A(NSXXB),A(NSYYB),A(NSZZB),A(NSXYB),A(NSYZB),
     &                                                   A(NSZXB))
  710     CONTINUE

          CALL MDDEL ('EMSSA')
          CALL MDDEL ('DENSA')
          CALL MDDEL ('EMSSB')
          CALL MDDEL ('DENSB')
          CALL MDDEL ('SIGXXA')
          CALL MDDEL ('SIGYYA')
          CALL MDDEL ('SIGZZA')
          CALL MDDEL ('SIGXYA')
          CALL MDDEL ('SIGYZA')
          CALL MDDEL ('SIGZXA')
          CALL MDDEL ('SIGXXB')
          CALL MDDEL ('SIGYYB')
          CALL MDDEL ('SIGZZB')
          CALL MDDEL ('SIGXYB')
          CALL MDDEL ('SIGYZB')
          CALL MDDEL ('SIGZXB')

        END IF

C     *****************************************************************
C     CLEAN UP STUFF FOR NEXT ELEMENT BLOCK
C     *****************************************************************

        CALL MDDEL ('ISRCHR')
        CALL MDDEL ('RSRCHR')
        CALL MDDEL ('LIST')
        CALL MDDEL ('XYZSRF')
        CALL MDDEL ('XYZPTS')
        CALL MDDEL ('ICONA')
        CALL MDDEL ('ICONB')
        CALL MDDEL ('NDLSTA')
        CALL MDDEL ('NDLSTB')
        CALL MDDEL ('SOLEA')
        CALL MDDEL ('SOLEB')
        IF (ISCHEM .EQ. 0)THEN
          CALL MDDEL ('SOLENA')
          CALL MDDEL ('NELTN')
        ELSE IF (ISCHEM .EQ. 1)THEN
          CALL MDDEL ('SOLENA')
          CALL MDDEL ('INVLNA')
          CALL MDDEL ('INVCN')
          CALL MDDEL ('CNTRA')
        ELSE IF (ISCHEM .EQ. 2)THEN
          CONTINUE
        ELSE IF (ISCHEM .EQ. 3)THEN
          CALL MDDEL ('CNTRA')
          CALL MDDEL ('INVLNA')
          CALL MDDEL ('INVCN')
          CALL MDDEL ('ICHKEL')
          CALL MDDEL ('SOLGRA')
        END IF
        CALL MDDEL ('STATUS')

        CALL MDSTAT (MNERRS, MNUSED)
        IF (MNERRS .NE. 0) THEN
           CALL MDEROR(NOUT)
        END IF

 50   CONTINUE

C     *****************************************************************
C     STOP COMMAND
C     *****************************************************************

C     CLOSE FILES AND STOP

      CALL BANNR2(84,'NORMAL',NTPOUT)
      CALL BANNR2(84,'EXIT',NTPOUT)
      call debug('CLSFIL')
      CALL CLSFIL

      call addlog (qainfo(1))
      call wrapup(qainfo(1))
      STOP

  270 FORMAT (3X,/'DATA FROM MESH "A" (MESH & SOLUTION) FILE - ',A,
     *        //,10X,'HEADING - ',A,/)
  280 FORMAT (10x,I1,'-DIMENSIONAL MODEL',/
     &       ,10X,'NUMBER OF NODES IN MESH A (nodesa) -      ',I9,/
     &       ,10X,'NUMBER OF ELEMENTS IN MESH A (numela) -   ',I9,/
     &       ,10X,'NUMBER OF ELEMENT BLOCKS IN A (nblksa) -  ',I9)
  290 FORMAT (3X,/,'DATA FROM MESH "B" (MESH) FILE - ',A,
     *        //,10X,'HEADING - ',A,/)
  300 FORMAT (10x,I1,'-DIMENSIONAL MODEL',/
     &       ,10X,'NUMBER OF NODES IN MESH B (nodesb) -      ',I9,/
     &       ,10X,'NUMBER OF ELEMENTS IN MESH B (numelb) -   ',I9,/
     &       ,10X,'NUMBER OF ELEMENT BLOCKS IN B (nblksb) -  ',I9)
  320 FORMAT (10X,'CORRESPONDS TO MESH-A ELEMENT BLOCK ID',I9)
  321 FORMAT (10X,'USING SEARCH BOX TOLERANCE            ',F14.6,/)
  330 FORMAT (/,10X,'WORKING ON MESH-B ELEMENT BLOCK       ',I9,'/',I9,/
     &       ,10X,'ELEMENT BLOCK ID                      ',I9)
  430 FORMAT (/,5X,'MESH-B NODE NUMBER ',I9,/
     &       ,5x,' ELEMENT BLOCK     ',I9,/
     &       ,5X,' WAS NOT FOUND IN MESH-A BY SRCHQ')
      END
      BLOCK DATA INITLZ

C     ******************************************************************

C     BLOCK DATA SUBROUTINE TO INITIALIZE VARIABLES STORED IN
C     NAMED COMMON BLOCKS

C     ******************************************************************

C...NOTE: Cannot include exodusII.inc in a block data routine.
      PARAMETER (MXSTLN=32)
c      CHARACTER*10 ELTYPE

      include 'header.blk'
      include 'ntpdat.blk'
      include 'contrl.blk'
      include 'amesh.blk'
      include 'bmesh.blk'
      include 'aexds1.blk'
      include 'schdat.blk'
      include 'tapes.blk'
      include 'toldat.blk'
      include 'varnpt.blk'
      include 'varept.blk'
      include 'debg.blk'
      include 'inival.blk'

      DATA HED/' '/
      DATA NOUT,NTPOUT,NTP2,NTP3,NTP4/
     1     6,7,12,13,14/
      DATA (IFILES(I),I=1,5)/5*0/
      DATA ISCHEM/1/
      DATA IDEF/1/
      DATA IDEBUG/0/
      DATA IACCU/0/
      DATA IXDIS,IYDIS,IZDIS,IXVEL,IYVEL,IZVEL/6*0/
      DATA ISXX,ISYY,ISZZ,ISXY,ISYZ,ISZX,IELMS,IDENS/8*0/
      DATA NUMELA,NODESA,NBLKSA,NDIMA,NELNDA/5*0/
      DATA NUMELB,NODESB,NBLKSB,NDIMB,NELNDB/5*0/
      DATA NQAREC,NVARGP,NVARNP,NVAREL/4*0/
c      DATA (ELTYPE(I),I=1,13)/'TRI3','TRI6','QUAD4','QUAD8','QUAD9',
c     1     'TETRA4','TETRA10','PRISM6','PRISM15','HEX8','HEX20',
c     2     'HEX27','SHELL'/
C      DATA (NNELM(I),I=1,13)/3,6,4,8,9,4,10,6,15,8,20,27,4/

      DATA VALINI /0.0/
      DATA TOLSHL,TOLQAD,TOLHEX,TOLTET,NISS,NRSS
     *  /0.01,0.01,0.01,0.01, 5,10/

C TOLSHL=extension of box around MESH-A shell element
C TOLQAD=extension of box around MESH-A quad element
C TOLHEX=extension of box around MESH-A hex element
C NISS=number of integer search scratch  (=5)
C NRSS=number of    real search scratch (=10)

      DATA TOL,EPS,STRLMT,ITERMX/0.01,0.01,20.,20/

C TOL=difference in isoparametric coords after newton iteration (don't change)
C EPS=tolerance used in checking if point is within element or coincident
C     with a node
C STRLMT=tolerance for isoparametric coords to lie within an element

      END

      subroutine getval(IMAP, IM, idblka, idblkb, ischem)
      integer imap(3,*)
      idblka = imap(1,im)
      idblkb = imap(2,im)
      ischem = imap(3,im)
      return
      end
