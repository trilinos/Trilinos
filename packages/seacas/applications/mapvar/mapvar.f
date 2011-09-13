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


C=======================================================================
*DECK,MAPVAR
      PROGRAM MAPVAR
C
C     ******************************************************************
C
C                                --MAPVAR--
C                  A PROGRAM TO MAP FINITE ELEMENT RESULTS
C                 FROM ONE EXODUS-II RESTART FILE TO ANOTHER
C                 EXODUS-II RESTART FILE TO SUPPORT REMESHING
C
C                             GERALD W. WELLMAN
C                        SANDIA NATIONAL LABORATORIES
C                          ALBUQUERQUE, NEW MEXICO
C
C     MAPVAR IS BASED ON MERLIN II,A FINITE ELEMENT INTERPOLATION 
C     PROGRAM BY DAVID K. GARTLING.
C
C     THE MERLIN PROGRAM IS DESIGNED TO TRANSFER DATA BETWEEN TWO- AND
C     THREE-DIMENSIONAL FINITE ELEMENT MESHES. GIVEN A FINITE ELEMENT
C     MESH, (MESH-A), AND A SOLUTION FIELD ON THAT MESH, MERLIN
C     INTERPOLATES THE GIVEN SOLUTION ONTO A SECOND FINITE ELEMENT MESH,
C     (MESH-B). THE INTERPOLATION PROCEDURE IS BASED ON USE OF THE
C     ELEMENT SHAPE FUNCTIONS IN MESH-A.
C
C     MAPVAR IS DESIGNED TO SERVE THE SAME GENERAL PURPOSE AS MERLIN
C     HOWEVER, MAPVAR IS DESIGNED TO PROVIDE THE TRANSLATION IN TERMS
C     OF EXODUS-II-V2 RESTART FILES. MAPVAR ALSO TRANSLATES ELEMENT 
C     VARIABLES AS WELL AS NODAL VARIABLES AND READS/WRITES GLOBAL 
C     VARIABLES. MAPVAR IS CURRENTLY SET UP FOR THREE ELEMENT TYPES,
C     2-D QUADS, 3-D HEXES, AND 3-D QUAD SHELLS. ALL THE ELEMENTS 
C     ORIGINALLY SUPPORTED BY MERLIN CAN BE INCLUDED IN MAPVAR GIVEN 
C     THE DESIRE AND RESOURCES. THE SEARCH ENGINE OF MAPVAR HAS BEEN
C     CHANGED TO A BINARY SEARCH FROM THE BIN OR BUCKET SEARCH OF MERLIN.
C
C     THE INTENT OF MAPVAR IS TO CREATE A RESTART FILE THAT WILL ALLOW
C     A FINITE ELEMENT SOLUTION TO PROCEED WITH A DIFFERENT MESH THAN 
C     THE MESH WITH WHICH THE SOLUTION WAS STARTED. THUS, THERE IS AN
C     INHERENT ASSUMPTION THAT MESH-A IS AN EXODUS RESTART FILE.
C
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
C
C     Special considerations:
C       ELMASS - translated to nodal density, interpolated, translated
C                back to ELMASS
C       ROTATIONS - magnitude of matrix must be 1, translate, compute
C                   magnitude, divide each component by magnitude
C       EQPS - constrained to be .gt. 0.
C       TEARING - constrained to be .gt. 0.
C       DECAY - constrained to be .lt. 1.
C
C       Code Tree:
C (SUPES, EXODUSII, and multiple calls to ERROR not included)
C 
C SET-UP
C                 MAPVAR - OPNFIL
C                          ERROR -  CLSFIL
C                          RDINPT - BANNER
C                                   ERROR
C                          RDA1 -   LENSTR
C                          RDB1
C                          RDA2 -   ERROR
C                          RDB2 -   ERROR
C SHELLS
C                          BLDSRF
C                          BLDPTN
C                          SRCHS -  MKRNK -  INDEXX
C                                            RANK
C                                -  GETBND - SRCHGE
C                                            SRCHGT
C                                -  MKLSTV
C                                -  SHLSRC
C                          SINTPN - SHAPEF
C                          BLDPTE - CNTR
C                          SRCHS - 
C      SCHEME 0
C                          SETON0 - VOL
C                          SINTPE - SHAPEF
C                                   VOL
C      SCHEME 1
C                          SETON1 - CNTR
C                                   VOL
C                                   INVCON
C                                   EXTS -   FRGE
C                                            BS
C                                   AVG
C                          SINTPE -
C      SCHEME 2
C                          STRAN
C      SCHEME 3
C                          INVCON
C                          ELGRAD - CNTR
C                                 - VOL
C                                 - FLGRAD
C                          INTRP3
C QUADS
C                          BLDSRF
C                          BLDPTN
C                          SRCHQ -  MKRNK -  INDEXX
C                                            RANK
C                                -  GETBND - SRCHGE
C                                            SRCHGT
C                                -  MKLSTV
C                                -  QADSRC - NODE
C                                            JACOBN
C                          INTRPN - SHAPEF
C                          BLDPTE -
C                          SRCHQ -
C      SCHEME 0 
C                          ELTON0 - VOL
C                          INTRPE - SHAPEF
C                                   VOL
C      SCHEME 1
C                          ELTON1 - CNTR
C                                   VOL
C                                   INVCON
C                                   EXTQ -   FRGE
C                                            BS
C                                   AVG
C                                   EXTH -   FRGE
C                                            BS
C                          INTRPE -
C      SCHEME 2
C                          TRANAB
C      SCHEME 3
C                          INVCON
C                          ELGRAD - CNTR
C                                 - VOL
C                                 - FLGRAD
C                          INTRP3
C HEXES
C                          BLDSRF
C                          BLDPTN
C                          SRCHH -   MKRNK -  INDEXX
C                                            RANK
C                                -  GETBND - SRCHGE
C                                            SRCHGT
C                                -  MKLSTV
C                                -  HEXSRC - NODE
C                                            JACOBN
C                          INTRPN -
C                          BLDPTE -
C                          SRCHH -
C      SCHEME 0
C                          ELTON0 -
C                          INTRPE -
C      SCHEME 1
C                          ELTON1 -
C                          INTRPE -
C      SCHEME 2
C                          TRANAB
C      SCHEME 3
C                          INVCON
C                          ELGRAD - CNTR
C                                 - VOL
C                                 - FLGRAD
C                          INTRP3
C WRAP-UP
C                          WRTC
C                          BANNER
C                          CLSFIL
C
C
C SUPES CALLS:
C              EXDATE, EXTIME, FREFLD, MDDEL, MDEROR, MDGET,
C              MDINIT, MDRSRV, MDSTAT, STRIPB
C
C EXODUSII CALLS:
C              EXCLOS, EXINQ,  (EXOPN)-fcn,    
C              EXGATM, EXGCON, EXGCOR, EXGEAT, EXGEBI, EXGELB,
C              EXGELC, EXGEV,  EXGGV,  EXGINI, EXGNSI, EXGNP,
C              EXGNV,  EXGNS,  EXGNSD, EXGP,   EXGPN,  EXGQA,
C              EXGSS,  EXGSSD, EXGSSI, EXGSP,  EXGTIM, EXGVAN,
C              EXGVP,
C              EXPCON, EXPCOR, EXPEAT, EXPELB, EXPELC, EXPEV,
C              EXPGV,  EXPINI, EXPNS,  EXPNSD, EXPNV,  EXPNP,
C              EXPP,   EXPPN,  EXPQA,  EXPSP,  EXPSS,  EXPSSD,
C              EXPTIM, EXPVAN, 
C
C     ******************************************************************
C
C     THE BASIC REFERENCE DOCUMENT FOR THIS CODE IS SAND 99-0466
C
C     ******************************************************************
C
C     COMPUTER CODE MANAGEMENT SYSTEM INFORMATION --
C
C     CURRENT VERSION DESIGNATOR- $Revision: 1.12 $
C
C     ******************************************************************
C
      include 'exodusII.inc'
C
      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'amesh.blk'
      include 'bmesh.blk'
      include 'contrl.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'header.blk'
      include 'ntpdat.blk'
      include 'rundat.blk'
      include 'steps.blk'
      include 'schdat.blk'
      include 'tapes.blk'
      include 'toldat.blk'
      include 'varnpt.blk'
      include 'varept.blk'
C
      EXTERNAL INITLZ
C
      DIMENSION A(1),IA(1)
      EQUIVALENCE (A(1),IA(1))
      CHARACTER*(MXSTLN) TYP
C
C     ******************************************************************
C
C     MAIN PROGRAM FOR MAPVAR
C     PROGRAM EXECUTION IS DIRECTED FROM THIS ROUTINE
C
C     ******************************************************************
C
C     NOTE : ALL ELEMENT DATA,NODAL POINT DATA, SOLUTION FIELDS, WORK
C            SPACE, ETC. ARE STORED IN THE ARRAY "A". THE MEMORY
C            MANAGER REQUESTS THE NEEDED STORAGE SPACE IN "A" DURING
C            EXECUTION. THE POINTERS NA1,NA2,....NAM PARTITION THE
C            ALLOCATED STORAGE INTO SEPARATE ARRAYS. SEE SAND86-0911,
C            "SUPES", FOR DETAILS OF THE MEMORY MANAGER.
C
C     ******************************************************************
C
C open all disk files
C
c      call excpus(timin)
c      write(nout,1001)
c      write(ntpout,1001)
C
      call init

C disable netcdf warning messages
C
      CALL EXOPTS(0,IERR)
C
      CALL OPNFIL
c      call excpus(timout)
c      timsub= timout - timin
c      write(nout,2001)timsub
c      write(ntpout,2001)timsub
c 2001 format('cpu time in subroutine',1x,e14.6)
C
C get info for QA records
C
      CALL VERSION(QAINFO)
C
C initialize memory manager
C
      CALL MDINIT (A)
C
C
C ******************************************************************
C read input parameters needed to control the solution from
C      the screen (time step(s), bins)
C ******************************************************************
C
      CALL EXINQ(NTP2EX,EXTIMS,NTIMES,RDUM,CDUM,IERR)
      CALL EXGINI (NTP2EX,HED,NDIMA,NODESA,NUMELA,NBLKSA,
     &             NUMNPS,NUMESS,IERR)
      CALL EXGINI (NTP3EX,HED,NDIMB,NODESB,NUMELB,NBLKSB,
     &             NUMNPS,NUMESS,IERR)
      MBLK = NBLKSA+NBLKSB
C
C A(NT1)      =   TIMES(1:NTIMES) - times on Mesh-A database
C IA(NAEB)    =   IDA(1:NBLKSA) - Donor mesh element block I.D.'s
C IA(NBEB)    =   IDB(1:NBLKSA) - Recipient mesh element block I.D.'s
C IA(NMAP)    =   MP(1:3,1:MBLK) - Donor to recipient mesh map
C
      CALL MDRSRV ('TIMES', NT1,   NTIMES)
      CALL MDRSRV ('IDA',   NAEB,  NBLKSA)
      CALL MDRSRV ('IDB',   NBEB,  NBLKSB)
      CALL MDRSRV ('MP',    NMAP,  MBLK*3)
C
      CALL MDSTAT (MNERRS, MNUSED)
      IF (MNERRS .NE. 0) THEN
         CALL MDEROR(NOUT)
         CALL ERROR('MAPVAR',
     &              'MEMORY MANAGER ERROR',
     &              'JUST BEFORE CALL TO RDINPT',0,' ',0,' ',' ',1)
      END IF
C
C
c      write(nout,1003)
c      write(ntpout,1003)
      CALL RDINPT (A(NT1),IA(NAEB),IA(NBEB),IA(NMAP),IMP)
C
C
C ******************************************************************
C INITIAL READ OF MESH-A
C ******************************************************************
C
C read sizing data for mesh A
C
      WRITE (NOUT, 270) HED
      WRITE (NTPOUT, 270) HED
      WRITE (NOUT, 280) NDIMA,NODESA,NUMELA,NBLKSA
      WRITE (NTPOUT, 280) NDIMA,NODESA,NUMELA,NBLKSA
C
C allocate storage for initial read of mesh A
C
C      A(NAX)     =    XA(1:NODESA) - Mesh-A X-coord
C      A(NAY)     =    YA(1:NODESA) - Mesh-A Y-coord
C      A(NAZ)     =    ZA(1:NODESA) - Mesh-A Z-coord
C      A(NADX)    =    DISXA(1:NODESA) - Mesh-A X-displ
C      A(NADY)    =    DISYA(1:NODESA) - Mesh-A Y-displ
C      A(NADZ)    =    DISZA(1:NODESA) - Mesh-A Z-displ
C
      CALL MDRSRV ('XA',     NAX,   NODESA)
      CALL MDRSRV ('DISXA',  NADX,  NODESA)
      CALL MDRSRV ('YA',     NAY,   NODESA)
      CALL MDRSRV ('DISYA',  NADY,  NODESA)
      IF (NDIMA.EQ.3) THEN
        CALL MDRSRV ('ZA',     NAZ,  NODESA)
        CALL MDRSRV ('DISZA',  NADZ, NODESA)
      ELSE
        CALL MDRSRV ('ZA',     NAZ,  1)
        CALL MDRSRV ('DISZA',  NADZ, 1)
      END IF
      CALL MDSTAT (MNERRS, MNUSED)
      IF (MNERRS .NE. 0) THEN
         CALL MDEROR(NOUT)
         CALL ERROR('MAPVAR',
     &              'MEMORY MANAGER ERROR',
     &              'JUST BEFORE CALL TO RDA1',0,' ',0,' ',' ',1)
      END IF
C
C
C Copy "GENESIS" from mesh-B to mesh-C
C
      CALL EXCOPY(NTP3EX,NTP4EX,IERR)
      IF (IERR .NE. 0)
     &  CALL ERROR('MAPVAR',
     &             'ERROR WITH EXCOPY - GENESIS FILE COPY',
     &             ' ',0,' ',0,' ',' ',1)
C
C
c read mesh A (coords,displ,variable names,QA, INFO records), 
c write mesh C, (variable names, QA, INFO records)
c
c
c      write(nout,1005)
c      write(ntpout,1005)
      CALL RDA1 (A(NAX),A(NAY),A(NAZ),A(NADX),A(NADY),A(NADZ))
C
      IF (IACCU .EQ. 1)THEN
        IF (IXVEL .NE. 0 .AND. IYVEL .NE. 0 .AND. 
     &      (IELMS .NE. 0 .OR. IDENS .NE. 0))THEN
C
C velocities and mass are available, compute momenta and k.e.
C 1st set up storage for vel's
C
C  A(NVXA)    = VELXA(1:NODESA) - X-velocity in mesh-A
C  A(NVYA)    = VELYA(1:NODESA) - Y-velocity in mesh-A
C  A(NVZA)    = VELZA(1:NODESA) - Z-velocity in mesh-A
C  A(NNMSA)   = RMSNA(1:NODESA) - nodal mass in mesh-A
C  A(NVXB)    = VELXB(1:NODESB) - X-velocity in mesh-B
C  A(NVYB)    = VELYB(1:NODESB) - Y-velocity in mesh-B
C  A(NVZB)    = VELZB(1:NODESB) - Z-velocity in mesh-B
C  A(NNMSB)   = RMSNB(1:NODESB) - nodal mass in mesh-B
C
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
C
C initialization quantities (6 each for now) that need to be 
C summed over the element blocks if doing an accuracy check
C
          IF (ISTEP .EQ. -1)THEN
C
C need arrays (one slot for each time), else just scalars will do
C
C A(NTMXA)  =  TMXA(1:NTIMES)  - X-momentum all mesh-A blocks each time
C A(NTMYA)  =  TMYA(1:NTIMES)  - Y-momentum all mesh-A blocks each time
C A(NTMZA)  =  TMZA(1:NTIMES)  - Z-momentum all mesh-A blocks each time
C A(NTKEA)  =  TKEA(1:NTIMES)  - K.E. all mesh-A blocks each time
C A(NTPSQA) =  TPSQA(1:NTIMES) - Pressure squared mesh-A each time
C A(NTJ2A)  =  TJ2A(1:NTIMES)  - J2 mesh-A each time
C
C A(NTMXB)  =  TMXB(1:NTIMES)  - X-momentum all mesh-B blocks each time
C A(NTMYB)  =  TMYB(1:NTIMES)  - Y-momentum all mesh-B blocks each time
C A(NTMZB)  =  TMZB(1:NTIMES)  - Z-momentum all mesh-B blocks each time
C A(NTKEB)  =  TKEB(1:NTIMES)  - K.E. all mesh-B blocks each time
C A(NTPSQB) =  TPSQB(1:NTIMES) - Pressure squared mesh-B each time
C A(NTJ2B)  =  TJ2B(1:NTIMES)  - J2 mesh-B each time
C
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
C
C
C ******************************************************************
C INITIAL READ OF MESH-B
C ******************************************************************
C
C
C read sizing data for mesh B
C
      WRITE (NOUT, 290) HED
      WRITE (NTPOUT, 290) HED
      WRITE (NOUT, 300) NDIMB,NODESB,NUMELB,NBLKSB
      WRITE (NTPOUT, 300) NDIMB,NODESB,NUMELB,NBLKSB
c
c quick initial check of compatibility mesh-A to mesh-B 
c
      IF (NDIMB .NE. NDIMA) THEN
        CALL ERROR('MAPVAR',
     &             'MESH-B INCOMPATIBLE WITH MESH-A',
     &             'DIMENSION OF MESH-A',NDIMA,
     &             'DIMENSION OF MESH-B',NDIMB,
     &             ' ',' ',1)
      END IF
C
C allocate storage for mesh B read
C
C      A(NBX)      =    XB(1:NODESB) - Mesh-B X-coord
C      A(NBY)      =    YB(1:NODESB) - Mesh-B Y-coord
C      A(NBZ)      =    ZB(1:NODESB) - Mesh-B Z-coord
C
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
C
C
c read coordinates for mesh B
c
c      write(nout,1009)
c      write(ntpout,1009)
      CALL RDB1 (A(NBX),A(NBY),A(NBZ))
C
C
C *********************************************************
C START INTERPOLATION
C *********************************************************
C
C set up memory for arrays for nodal results and truth table
C these arrays stay around forever - they don't get deleted
C after each element block is processed like the arrays 
C set up within the element block loop
C
C   A(NASOLN)    =    SOLNA(1:NODESA,1:NVARNP) - Mesh-A nodal data
C   A(NBSOLN)    =    SOLNB(1:NODESB,1:NVARNP) - Mesh-B interpolated  
C                                                nodal data
C   IA(ITTA)     =    ITRTA(1:NVAREL,1:NBLKSA) - Mesh-A truth table
C   IA(ITTB)     =    ITRTB(1:NVAREL,1:NBLKSB) - Mesh-B truth table
C
C    A(NSN)      = SN(1:NODESB)  - storage for nodal vars in ininod
C    A(NSE)      = SE(1:NODESB)  - storage for element vars in ininod
C

      CALL MDRSRV ('SOLNA',  NASOLN,   NODESA*NVARNP)
      CALL MDRSRV ('SOLNB',  NBSOLN,   NODESB*NVARNP)
      CALL MDRSRV ('ITRTA',  ITTA,     NVAREL*NBLKSA)
      CALL MDRSRV ('ITRTB',  ITTB,     NVAREL*NBLKSB)
      CALL MDRSRV ('SN',     NSN,      NODESB)
      CALL MDRSRV ('SE',     NSE,      NODESB)
C
      CALL MDSTAT (MNERRS, MNUSED)
      IF (MNERRS .NE. 0) THEN
         CALL MDEROR(NOUT)
         CALL ERROR('MAPVAR',
     &              'MEMORY MANAGER ERROR',
     &              'JUST BEFOR INTERPOLATION LOOP',
     &              0,' ',0,' ',' ',1)
      END IF
C
c      write(nout,1010)
c      write(ntpout,1010)
      CALL TRUTBL(IA(NMAP),IMP,IA(NAEB),IA(NBEB),IA(ITTA),IA(ITTB))
C
C *********************************************************************
C
C START OF ELEMENT BLOCK-BY-ELEMENT BLOCK INTERPOLATION LOOP
C
C *********************************************************************
C
      DO 50 IM = 1, IMP
        IMOFF = IM * 3
        IDBLKA = IA(NMAP-3+IMOFF)
        IDBLKB = IA(NMAP-2+IMOFF)
        ISCHEM = IA(NMAP-1+IMOFF)

        do 15 i=1, nblksa
          if (idblka .eq. ia(naeb-1+i)) then
            iblka = i
            go to 16
          endif
 15     continue
 16     continue
        
        do 25 i=1, nblksb
          if (idblkb .eq. ia(nbeb-1+i)) then
            iblkb = i
            go to 26
          endif
 25     continue
 26     continue
          
C
C set up controls for many to 1 map
C if first time recipient mesh element block called, insub = 1
C else insub = 2
C if last time recipient mesh element block called, icompl = 1
C else icompl = 0
C
        INSUB  = 1
        ICOMPL = 1
        IF (IM .GT. 1)THEN
          IDBBM1 = IA(NMAP-5+IMOFF)
          IF (IDBLKB .EQ. IDBBM1)THEN
            INSUB = 2
          END IF
        END IF
        IF (IM .LT. IMP)THEN
          IDBBP1 = IA(NMAP+1+IMOFF)
          IF (IDBLKB .EQ. IDBBP1)THEN
            ICOMPL = 0
          END IF
        END IF
C
C
C     **********************************************************
C     ELEMENT BLOCK BY ELEMENT BLOCK INTERPOLATION
C     REQUIRED FOR ELEMENT DATA BUT ALSO USED FOR NODAL DATA
C     **********************************************************
C
        WRITE(NOUT,330)IM,IDBLKB
        WRITE(NOUT,320)IDBLKA
        WRITE(NTPOUT,330)IM,IDBLKB
        WRITE(NTPOUT,320)IDBLKA
C
        CALL EXGELB(NTP2EX,IDBLKA,TYP,NUMEBA,NELNDA,NATRIB,
     &              IERR)
        CALL EXGELB(NTP3EX,IDBLKB,TYP,NUMEBB,NELNDB,NATRIB,
     &              IERR)
        IF (NUMEBB .EQ. 0)THEN
          GO TO 50
        END IF
C
C set up arrays for element block-by-element block preliminaries
C these arrays will be deleted at the end of the loop
C
C
C IA(NACON)   =    ICONA(1:NELNDA,1:NUMEBA) - Mesh-A connectivity
C IA(NBCON)   =    ICONB(1:NELNDB,1:NUMEBB) - Mesh-B connectivity
C IA(NANDLST) =    NDLSTA(1:NODESA) - Mesh-A nodes in element block
C IA(NBNDLST) =    NDLSTB(1:NODESB) - Mesh-A nodes in element block
C  A(NASTAT)  =    STATUS(1:NUMEBA) - Mesh-A element status
C
        CALL MDRSRV ('ICONA',  NACON,   NELNDA*NUMEBA)
        CALL MDRSRV ('ICONB',  NBCON,   NELNDB*NUMEBB)
        CALL MDRSRV ('NDLSTA', NANDLST, NODESA)
        CALL MDRSRV ('NDLSTB', NBNDLST, NODESB)
        CALL MDRSRV ('STATUS', NASTAT,  NUMEBA)
        CALL MDSTAT (MNERRS, MNUSED)
C
        IF (MNERRS .NE. 0) THEN
          CALL MDEROR(NOUT)
          CALL ERROR('MAPVAR',
     &               'MEMORY MANAGER ERROR',
     &               'BLOCKS LOOP PRELIMINARIES',
     &               0,' ',0,' ',' ',1)
        END IF
C
c 2nd read of mesh A
c
c        write(nout,1011)
c        write(ntpout,1011)
        CALL RDA2 (IDBLKA,IA(NACON),IA(NANDLST),A(NASTAT),
     &             MAXLN)
c
c 2nd read of mesh-b
c
c        write(nout,1013)
c        write(ntpout,1013)
        CALL RDB2(IDBLKB,IDBLKA,IA(NBCON),IA(NBNDLST))
C
C
C set up arrays for element block-by-element block processing
C these arrays will be deleted at the end of the loop
C
C IA(NS1)     =  ISRCHR(1:1(NISR),1:NUMNDB) Integer search results
C  A(NS2)     =  RSRCHR(1:6(NRSR),1:NUMNDB) Real search results
C IA(NS3)     =    LIST(1:NUMNDB)         Potential contacts 
C IA(NS4)     =     IND(1:NUMNDB,1:3)     Index array point order
C IA(NS5)     =    IRNK(1:NUMNDB,1:3)     Rank array
C IA(NS6)     =   IRNK2(1:NUMNDB,1:3,1:2) Indirect rank array
C IA(NS7)     =    INDX(1:NUMNDB)         Intermediate list potential
C                                         pairs (intersection lists)
C IA(NS8)     =     ILO(1:LBLK,1:3)       Min index search box
C IA(NS9)     =     IUP(1:LBLK,1:3)       Max index search box
C IA(NS10)    =     IDP(1:LBLK)           Points paired with element
C IA(NS11)    =     IDS(1:LBLK)           Elements paired with point
C  A(NS12)    =    XMIN(1:LBLK,1:3)       Min dim search box
C  A(NS13)    =    XMAX(1:LBLK,1:3)       Max dim search box
C IA(NS14)    =    ISCR(1:NISS,1:LBLK)    Integer scratch          
C  A(NS15)    =    RSCR(1:NRSS,1:LBLK)    Real scratch             
C  A(NS16)    =    XYZSRF(1:NODESA,1:3)   Coords defining element  
C  A(NS17)    =    XYZPTS(1:NUMNDB,1:3)   Coords of points searched
C
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
        CALL MDRSRV ('ISRCHR', NS1,  1*idim)
        CALL MDRSRV ('RSRCHR', NS2,  6*idim)
        CALL MDRSRV ('LIST',   NS3,  idim)
        CALL MDRSRV ('IND',    NS4,  idim*3)
        CALL MDRSRV ('IRNK',   NS5,  idim*3)
        CALL MDRSRV ('IRNK2',  NS6,  idim*6)
        CALL MDRSRV ('INDX',   NS7,  idim)
        CALL MDRSRV ('ILO',    NS8,  LBLK*3)
        CALL MDRSRV ('IUP',    NS9,  LBLK*3)
        CALL MDRSRV ('IDP',    NS10, LBLK)
        CALL MDRSRV ('IDS',    NS11, LBLK)
        CALL MDRSRV ('XMIN',   NS12, LBLK*3)
        CALL MDRSRV ('XMAX',   NS13, LBLK*3)
        CALL MDRSRV ('ISCR',   NS14, NISS*LBLK)
        CALL MDRSRV ('RSCR',   NS15, NRSS*LBLK)
        CALL MDRSRV ('XYZSRF', NS16, NODESA*3)
        CALL MDRSRV ('XYZPTS', NS17, idim*3)
C
        CALL MDRSRV ('SOLEA',  NASOLE,  NUMEBA*NVAREL)
        CALL MDRSRV ('SOLEB',  NBSOLE,  NUMEBB*NVAREL)
C
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
C
        CALL MDSTAT (MNERRS, MNUSED)
        IF (MNERRS .NE. 0) THEN
          CALL MDEROR(NOUT)
          CALL ERROR('MAPVAR',
     &               'MEMORY MANAGER ERROR',
     &               'JUST AFTER START OF BLOCKS LOOP',
     &               0,' ',0,' ',' ',1)
        END IF
C
C
        IF (ITYPE .EQ. 13)THEN
C
C     **********************************************************
C     Path through code for shells
C     **********************************************************
C
c          write(nout,1014)
c          write(ntpout,1014)
          CALL BLDSRF(A(NAX),A(NAY),A(NAZ),A(NS16))
C
          IF (NVARNP .GT. 0)THEN
C
c            write(nout,1015)
c            write(ntpout,1015)
            CALL BLDPTN(A(NBX),A(NBY),A(NBZ),IA(NBNDLST),A(NS17))
C
c            write(nout,1017)
c            write(ntpout,1017)
            CALL SRCHS (NODESA,NUMEBA,IA(NACON),A(NS16),
     1       NUMNDB,A(NS17),TOLSHL,1,6,
     2       NISS,NRSS,IA(NS1),A(NS2),LBLK,
     3       IA(NS3),IA(NS4),IA(NS5),IA(NS6),IA(NS7),IA(NS8),IA(NS9),
     4       IA(NS10),IA(NS11),A(NS12),A(NS13),IA(NS14),A(NS15),IERR)
C
c            write(nout,1019)
c            write(ntpout,1019)
            CALL SINTPN(IA(NACON),A(NASOLN),IA(NS1),1,A(NS2),6,
     1                  A(NBSOLN),IA(NBNDLST),A(NBX),A(NBY),A(NBZ),
     2                  IDBLKB,A(NT1),INSUB,A(NSN))
C
          END IF
          IF (NVAREL .GT. 0)THEN
C
c            write(nout,1021)
c            write(ntpout,1021)
            CALL BLDPTE(A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NS17))
C
c            write(nout,1023)
c            write(ntpout,1023)
            CALL SRCHS (NODESA,NUMEBA,IA(NACON),A(NS16),
     1       NUMEBB,A(NS17),TOLSHL,1,6,
     2       NISS,NRSS,IA(NS1),A(NS2),LBLK,
     3       IA(NS3),IA(NS4),IA(NS5),IA(NS6),IA(NS7),IA(NS8),IA(NS9),
     4       IA(NS10),IA(NS11),A(NS12),A(NS13),IA(NS14),A(NS15),IERR)
C
c element centroid values to nodes by averaging
c
            IF (ISCHEM .EQ. 0)THEN
C
C Set up time steps
C
              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF
C
              DO 610 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF
C     
c                write(nout,1025)
c                write(ntpout,1025)
                CALL SETON0(IA(NACON),IA(NANELTN),A(NASOLE),
     &           A(NASOLEN),IDBLKA,A(NAX),A(NAY),A(NAZ),ISTP,
     &           IA(ITTB),IBLKB)
C
c                write(nout,1027)
c                write(ntpout,1027)
                CALL SINTPE(IA(NACON),A(NASOLEN),IA(NS1),1,A(NS2),6,
     &                      A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     &                      IA(NBCON),IA(ITTB),IBLKB,A(NT1),A(NS17),
     &                      ISTP,IST,INSUB,ICOMPL,A(NSE))
  610         CONTINUE
C
            ELSE IF (ISCHEM .EQ. 1) THEN
C
c              write(nout,1028)
c              write(ntpout,1028)
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))
C
C Set up time steps
C
              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF
C
              DO 620 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF
C     
c                write(nout,1029)
c                write(ntpout,1029)
                CALL SETON1(A(NACTR),A(NASOLE),A(NASOLEN),IDBLKA,
     &                    A(NAX),A(NAY),A(NAZ),IA(NACON),IA(NANDLST),
     &                    IA(NAINVLN),IA(NAINVC),MAXLN,ISTP,
     &                    IA(ITTB),IBLKB)
C
c                write(nout,1027)
c                write(ntpout,1027)
                CALL SINTPE(IA(NACON),A(NASOLEN),IA(NS1),1,A(NS2),6,
     &                      A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     &                      IA(NBCON),IA(ITTB),IBLKB,A(NT1),A(NS17),
     &                      ISTP,IST,INSUB,ICOMPL,A(NSE))
  620         CONTINUE
C
            ELSE IF (ISCHEM .EQ. 2) THEN
c
c direct transfer, does not require scatter to nodes
c
c              write(nout,1031)
c              write(ntpout,1031)
              CALL STRAN(IA(NS1),1,A(NASOLE),A(NBSOLE),
     &                    IDBLKA,IDBLKB,
     &                    IA(ITTB),IBLKB,A(NT1),A(NS17),
     &                    INSUB,ICOMPL,
     &                    A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NSE))
C
            ELSE IF (ISCHEM .EQ. 3)THEN
C
c              write(nout,1028)
c              write(ntpout,1028)
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))
C
C Set up time steps
C
              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF
C
              DO 630 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF
C
c                write(nout,1065)
c                write(ntpout,1065)
                CALL ELGRAD(A(NACTR),A(NAX),A(NAY),A(NAZ),
     &                      A(NASOLE),A(NSOLGR),IA(NICHKE),
     &                      IDBLKA,IA(NACON),IA(NAINVLN),IA(NAINVC),
     &                      MAXLN,ISTP,IA(ITTB),IBLKB)
C
c                write(nout,1067)
c                write(ntpout,1067)
                CALL INTRP3(A(NACTR),A(NS17),IA(NS1),
     &                      A(NBSOLE),A(NASOLE),A(NSOLGR),
     &                      IDBLKB,IA(ITTB),IBLKB,A(NT1),
     &                      ISTP,IST,INSUB,ICOMPL,
     &                      A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NSE))
 630          CONTINUE
C
            ELSE
              CALL ERROR('MAPVAR',' ','ISCHEM =',
     &                 ischem,'INCORRECT ARGUMENT',0,' ',' ',1)
            END IF
C
          END IF
C
C ITYPE = 3 - 4 node quad
C ITYPE = 4 - 8 node quad
C ITYPE = 5 - 9 node quad
C
        ELSE IF (ITYPE .EQ. 3 .OR. ITYPE .EQ. 4 .OR. 
     &           ITYPE .EQ. 5) THEN
C
C
C     *****************************************************
C     PATH THROUGH CODE FOR CONTINUUM ELEMENTS
C              (QUAD-4)
C     *****************************************************
C
C find and store location of mesh-b nodes within mesh-a
C
c          write(nout,1014)
c          write(ntpout,1014)
          CALL BLDSRF(A(NAX),A(NAY),A(NAZ),A(NS16))
C
          IF (NVARNP .GT. 0)THEN
C
c            write(nout,1015)
c            write(ntpout,1015)
            CALL BLDPTN(A(NBX),A(NBY),A(NBZ),IA(NBNDLST),A(NS17))
C   
c            write(nout,1035)
c            write(ntpout,1035)
            CALL SRCHQ (NODESA,NUMEBA,IA(NACON),A(NS16),
     1       NUMNDB,A(NS17),TOLQAD,1,3,
     2       NISS,NRSS,IA(NS1),A(NS2),LBLK,
     3       IA(NS3),IA(NS4),IA(NS5),IA(NS6),IA(NS7),IA(NS8),IA(NS9),
     4       IA(NS10),IA(NS11),A(NS12),A(NS13),IA(NS14),A(NS15),IERR)
c           DO 530 I = 1, NUMNDB
c             IF (IA(NS1-1+I) .EQ. 0)THEN
c               WRITE(NOUT,430)IA(NBNDLST-1+I),IDBLKB
c               WRITE(NTPOUT,430)IA(NBNDLST-1+I),IDBLKB
c             END IF
c 530       CONTINUE
C   
c   interpolate nodal variables
c   
c            write(nout,1037)
c            write(ntpout,1037)
            CALL INTRPN(IA(NACON),A(NASOLN),IA(NS1),A(NS2),
     &                A(NBSOLN),IA(NBNDLST),A(NBX),A(NBY),A(NBZ),
     &                IDBLKB,A(NT1),INSUB,A(NSN))
c
c start element variable interpolation
c
c
c locate Mesh-B element centroid in Mesh-A
c
          END IF
          IF (NVAREL .GT. 0)THEN
C
c            write(nout,1021)
c            write(ntpout,1021)
            CALL BLDPTE(A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NS17))
C
c            write(nout,1039)
c            write(ntpout,1039)
            CALL SRCHQ (NODESA,NUMEBA,IA(NACON),A(NS16),
     1       NUMEBB,A(NS17),TOLQAD,1,3,
     2       NISS,NRSS,IA(NS1),A(NS2),LBLK,
     3       IA(NS3),IA(NS4),IA(NS5),IA(NS6),IA(NS7),IA(NS8),IA(NS9),
     4       IA(NS10),IA(NS11),A(NS12),A(NS13),IA(NS14),A(NS15),IERR)
C
c
c element centroid variables averaged to nodes
c
            IF (ISCHEM .EQ. 0)THEN
C
C Set up time steps
C
              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF
C
              DO 640 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF
C
c                write(nout,1041)
c                write(ntpout,1041)
                CALL ELTON0(IA(NACON),IA(NANELTN),A(NASOLE),
     &           A(NASOLEN),IDBLKA,A(NAX),A(NAY),A(NAZ),ISTP,
     &           IA(ITTB),IBLKB)
c
c interpolate element vars
c
c                write(nout,1043)
c                write(ntpout,1043)
                CALL INTRPE(IA(NACON),A(NASOLEN),IA(NS1),A(NS2),
     1                    A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     2                    IA(NBCON),IA(ITTB),IBLKB,A(NT1),
     3                    A(NS17),ISTP,IST,INSUB,ICOMPL,A(NSE))
  640         CONTINUE
C
c element centroid variables linear least squares to nodes
c
            ELSE IF (ISCHEM .EQ. 1)THEN
C
c              write(nout,1028)
c              write(ntpout,1028)
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))
C
C Set up time steps
C
              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF
C
              DO 650 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF
C
c                write(nout,1045)
c                write(ntpout,1045)
                CALL ELTON1(A(NACTR),A(NASOLE),A(NASOLEN),IDBLKA,
     &                    A(NAX),A(NAY),A(NAZ),IA(NACON),IA(NANDLST),
     &                    IA(NAINVLN),IA(NAINVC),MAXLN,ISTP,
     &                    IA(ITTB),IBLKB)
c
c interpolate element vars
c
c                write(nout,1043)
c                write(ntpout,1043)
                CALL INTRPE(IA(NACON),A(NASOLEN),IA(NS1),A(NS2),
     1                    A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     2                    IA(NBCON),IA(ITTB),IBLKB,A(NT1),
     3                    A(NS17),ISTP,IST,INSUB,ICOMPL,A(NSE))
  650         CONTINUE
C
            ELSE IF (ISCHEM .EQ. 2)THEN
C
c direct transfer from Mesh-A to Mesh-B
c
c              write(nout,1047)
c              write(ntpout,1047)
              CALL TRANAB(IA(NS1),A(NASOLE),A(NBSOLE),
     &                  IDBLKA,IDBLKB,
     &                  IA(ITTB),IBLKB,A(NT1),A(NS17),
     &                  INSUB,ICOMPL,
     &                  A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NSE))
c
            ELSE IF (ISCHEM .EQ. 3)THEN
C
c              write(nout,1028)
c              write(ntpout,1028)
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))
C
C Set up time steps
C
              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF
C
              DO 660 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF
C 
c                write(nout,1065)
c                write(ntpout,1065)
                CALL ELGRAD(A(NACTR),A(NAX),A(NAY),A(NAZ),
     &                    A(NASOLE),A(NSOLGR),IA(NICHKE),
     &                    IDBLKA,IA(NACON),IA(NAINVLN),IA(NAINVC),
     &                    MAXLN,ISTP,IA(ITTB),IBLKB)
C
c                write(nout,1067)
c                write(ntpout,1067)
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
C
          END IF
C
        ELSE IF (ITYPE .EQ. 10 .OR. ITYPE .EQ. 6) THEN
C
C
C     *****************************************************
C     PATH THROUGH CODE FOR 3-D CONTINUUM ELEMENTS
C              (HEX-8) OR (TET-8)
C     *****************************************************
C
C     FIND AND STORE LOCATION OF MESH-B NODES WITHIN MESH-A
C
c          write(nout,1014)
c          write(ntpout,1014)
          CALL BLDSRF(A(NAX),A(NAY),A(NAZ),A(NS16))
C
          IF (NVARNP .GT. 0)THEN
C
c             write(nout,1015)
c             write(ntpout,1015)
             CALL BLDPTN(A(NBX),A(NBY),A(NBZ),IA(NBNDLST),A(NS17))
C
            IF (ITYPE .EQ. 10)THEN
c              write(nout,1049)
c              write(ntpout,1049)
              CALL SRCHH (NODESA,NUMEBA,IA(NACON),A(NS16),
     1          NUMNDB,A(NS17),TOLHEX,1,3,
     2          NISS,NRSS,IA(NS1),A(NS2),LBLK,
     3          IA(NS3),IA(NS4),IA(NS5),IA(NS6),IA(NS7),IA(NS8),
     4          IA(NS9),
     5          IA(NS10),IA(NS11),A(NS12),A(NS13),IA(NS14),A(NS15),
     6          IERR)
C
            ELSEIF (ITYPE .EQ. 6)THEN
c              write(nout,1050)
c              write(ntpout,1050)
              CALL SRCHT (NODESA,NUMEBA,IA(NACON),A(NS16),
     1          NUMNDB,A(NS17),TOLHEX,1,3,
     2          NISS,NRSS,IA(NS1),A(NS2),LBLK,
     3          IA(NS3),IA(NS4),IA(NS5),IA(NS6),IA(NS7),IA(NS8),
     4          IA(NS9),
     5          IA(NS10),IA(NS11),A(NS12),A(NS13),IA(NS14),A(NS15),
     6          IERR)
            END IF
C
c interpolate nodal variables
c
c            write(nout,1051)
c            write(ntpout,1051)
            CALL INTRPN(IA(NACON),A(NASOLN),IA(NS1),A(NS2),
     &                A(NBSOLN),IA(NBNDLST),A(NBX),A(NBY),A(NBZ),
     &                IDBLKB,A(NT1),INSUB,A(NSN))
c
c start element variable interpolation
c
c locate Mesh-B element centroid in Mesh-A
c
          END IF
          IF (NVAREL .GT. 0)THEN
C
c            write(nout,1021)
c            write(ntpout,1021)
            CALL BLDPTE(A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NS17))
C
            IF (ITYPE .EQ. 10)THEN
c              write(nout,1053)
c              write(ntpout,1053)
              CALL SRCHH (NODESA,NUMEBA,IA(NACON),A(NS16),
     1          NUMEBB,A(NS17),TOLHEX,1,3,
     2          NISS,NRSS,IA(NS1),A(NS2),LBLK,
     3          IA(NS3),IA(NS4),IA(NS5),IA(NS6),IA(NS7),IA(NS8),
     4          IA(NS9),
     5          IA(NS10),IA(NS11),A(NS12),A(NS13),IA(NS14),A(NS15),
     6          IERR)
C
            ELSEIF (ITYPE .EQ. 6)THEN
c              write(nout,1054)
c              write(ntpout,1054)
              CALL SRCHT (NODESA,NUMEBA,IA(NACON),A(NS16),
     1          NUMEBB,A(NS17),TOLHEX,1,3,
     2          NISS,NRSS,IA(NS1),A(NS2),LBLK,
     3          IA(NS3),IA(NS4),IA(NS5),IA(NS6),IA(NS7),IA(NS8),
     4          IA(NS9),
     5          IA(NS10),IA(NS11),A(NS12),A(NS13),IA(NS14),A(NS15),
     6          IERR)
            END IF
C
            IF (ISCHEM .EQ. 0)THEN
C
C Set up time steps
C
              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF
C
              DO 670 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF
C     
c                write(nout,1055)
c                write(ntpout,1055)
                CALL ELTON0(IA(NACON),IA(NANELTN),A(NASOLE),
     &               A(NASOLEN),IDBLKA,A(NAX),A(NAY),A(NAZ),ISTP,
     &               IA(ITTB),IBLKB)
C
c interpolate element vars
c
c                write(nout,1043)
c                write(ntpout,1043)
                CALL INTRPE(IA(NACON),A(NASOLEN),IA(NS1),A(NS2),
     1                    A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     2                    IA(NBCON),IA(ITTB),IBLKB,A(NT1),
     3                    A(NS17),ISTP,IST,INSUB,ICOMPL,A(NSE))
  670         CONTINUE
C
            ELSE IF (ISCHEM .EQ. 1)THEN
C
c              write(nout,1028)
c              write(ntpout,1028)
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))
C
C Set up time steps
C
              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF
C
              DO 680 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF
C     
c                write(nout,1045)
c                write(ntpout,1045)
                CALL ELTON1(A(NACTR),A(NASOLE),A(NASOLEN),IDBLKA,
     &                    A(NAX),A(NAY),A(NAZ),IA(NACON),IA(NANDLST),
     &                    IA(NAINVLN),IA(NAINVC),MAXLN,ISTP,
     &                    IA(ITTB),IBLKB)
C
c interpolate element vars
c
c                write(nout,1043)
c                write(ntpout,1043)
                CALL INTRPE(IA(NACON),A(NASOLEN),IA(NS1),A(NS2),
     1                    A(NBSOLE),IDBLKB,A(NBX),A(NBY),A(NBZ),
     2                    IA(NBCON),IA(ITTB),IBLKB,A(NT1),
     3                    A(NS17),ISTP,IST,INSUB,ICOMPL,A(NSE))
  680         CONTINUE
c
            ELSE IF (ISCHEM .EQ. 2)THEN
C
c direct transfer from Mesh-A to Mesh-B
c
c              write(nout,1047)
c              write(ntpout,1047)
              CALL TRANAB(IA(NS1),A(NASOLE),A(NBSOLE),
     &                IDBLKA,IDBLKB,
     &                IA(ITTB),IBLKB,A(NT1),A(NS17),
     &                INSUB,ICOMPL,
     &                A(NBX),A(NBY),A(NBZ),IA(NBCON),A(NSE))
c
            ELSE IF (ISCHEM .EQ. 3)THEN
C
c              write(nout,1028)
c              write(ntpout,1028)
              CALL INVCON(IA(NAINVLN),MAXLN,IA(NAINVC),IA(NACON))
C
C Set up time steps
C
              IF (ISTEP .EQ. -1)THEN
                NTM = NTIMES
              ELSE
                NTM = 1
              END IF
C
              DO 690 IST = 1, NTM
                IF (ISTEP .EQ. -1)THEN
                  ISTP = IST
                ELSE
                  ISTP = ISTEP
                END IF
C
c                write(nout,1065)
c                write(ntpout,1065)
                CALL ELGRAD(A(NACTR),A(NAX),A(NAY),A(NAZ),
     &                    A(NASOLE),A(NSOLGR),IA(NICHKE),
     &                    IDBLKA,IA(NACON),IA(NAINVLN),IA(NAINVC),
     &                    MAXLN,ISTP,IA(ITTB),IBLKB)
C
c                write(nout,1067)
c                write(ntpout,1067)
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
C
          END IF
        ELSE
          CALL ERROR ('MAPVAR','INCORRECT ELEMENT TYPE',
     &              'ELEMENT TYPE =',ITYPE,
     &              'NOT YET IMPLEMENTED',0,' ',' ',1)
C
        END IF
C
        IF (IACCU .EQ. 1)THEN
C
C velocities and mass are available, compute momenta and k.e.
C 1st set up storage for mass
C
          IF (IELMS .NE. 0 .AND. IDENS .EQ. 0)THEN
C
C  A(NEMSA)   = EMSSA(1:NODESA) - element mass in mesh-A
C  A(NEMSB)   = EMSSB(1:NODESB) - element mass in mesh-B
C
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
C
C  A(NEMSA)   = EMSSA(1:NUMEBA) - element mass in mesh-A
C  A(NEMSB)   = EMSSB(1:NUMEBB) - element mass in mesh-B
C  A(NDENA)   = DENSA(1:NUMEBA) - element density in mesh-A
C  A(NDENB)   = DENSB(1:NUMEBB) - element density in mesh-B
C
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
C
          IF (NDIMA .EQ. 3)THEN
C
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
C
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
C
C  A(NSXXA)  = SIGXXA(1:NUMEBA) - XX component of stress tensor
C  A(NSYYA)  = SIGYYA(1:NUMEBA) - YY component of stress tensor
C  A(NSZZA)  = SIGZZA(1:NUMEBA) - ZZ component of stress tensor
C  A(NSXYA)  = SIGXYA(1:NUMEBA) - XY component of stress tensor
C  A(NSXXB)  = SIGXXB(1:NUMEBB) - XX component of stress tensor
C  A(NSYYB)  = SIGYYB(1:NUMEBB) - YY component of stress tensor
C  A(NSZZB)  = SIGZZB(1:NUMEBB) - ZZ component of stress tensor
C  A(NSXYB)  = SIGXYB(1:NUMEBB) - XY component of stress tensor
C
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
C
C Set up time steps
C
          IF (ISTEP .EQ. -1)THEN
            NTM = NTIMES
          ELSE
            NTM = 1
          END IF
C
          DO 710 IST = 1, NTM
            IF (ISTEP .EQ. -1)THEN
              ISTP = IST
            ELSE
              ISTP = ISTEP
            END IF
C
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
C
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
C
        END IF
C
C
C     *****************************************************************
C     CLEAN UP STUFF FOR NEXT ELEMENT BLOCK
C     *****************************************************************
C
        CALL MDDEL ('ISRCHR')
        CALL MDDEL ('RSRCHR')
        CALL MDDEL ('LIST')
        CALL MDDEL ('IND')
        CALL MDDEL ('IRNK')
        CALL MDDEL ('IRNK2')
        CALL MDDEL ('INDX')
        CALL MDDEL ('ILO')
        CALL MDDEL ('IUP')
        CALL MDDEL ('IDP')
        CALL MDDEL ('IDS')
        CALL MDDEL ('XMIN')
        CALL MDDEL ('XMAX')
        CALL MDDEL ('ISCR')
        CALL MDDEL ('RSCR')
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
C
c
 50   CONTINUE
c      
C     A(NAGV)    =    GVAR(1:NVARGP) Global variables
C
      CALL MDRSRV ('GVAR',   NAGV, NVARGP)
C
c      write(nout,1061)
c      write(ntpout,1061)
      CALL WRTC(A(NBX),A(NBY),A(NBZ),A(NAGV),A(NBSOLN))
C
C
C
C     *****************************************************************
C     STOP COMMAND
C     *****************************************************************
C
C     CLOSE FILES AND STOP
C
c
      CALL BANNR2(84,'NORMAL',NTPOUT)
      CALL BANNR2(84,'EXIT',NTPOUT)
c      write(nout,1063)
c      write(ntpout,1063)
      CALL CLSFIL
C
      call addlog (qainfo(1))
      call wrapup(qainfo(1))
      STOP
C
  210 FORMAT (//,3X,'DATA FOR RECTANGULAR GRID USED IN COARSE SEARCH -',
     1//,10X,'NUMBER OF ZONES (X COORDINATE) - ',I7,/,10X,'NUMBER OF ZON
     2ES (Y COORDINATE) - ',I7,/,10X,'NUMBER OF ZONES (Z COORDINATE) - '
     3,I7,/,10X,'TOTAL NUMBER OF ZONES (BINS) -   ',I7,/,10X,'NUMBER OF 
     4ELEMENTS ALLOWED IN EACH BIN - ',I7,/)
  270 FORMAT (3X,'DATA FROM MESH "A" FILE-',//,10X,'HEADING - ',A,/)
  280 FORMAT (10x,I1,'-DIMENSIONAL MODEL',/
     &       ,10X,'NUMBER OF NODES IN MESH A (nodesa) -      ',I7,/
     &       ,10X,'NUMBER OF ELEMENTS IN MESH A (numela) -   ',I7,/
     &       ,10X,'NUMBER OF ELEMENT BLOCKS IN A (nblksa) -  ',I7)
  290 FORMAT (3X,'DATA FROM MESH "B" FILE-',//,10X,'HEADING - ',A,/)
  300 FORMAT (10x,I1,'-DIMENSIONAL MODEL',/
     &       ,10X,'NUMBER OF NODES IN MESH B (nodesb) -      ',I7,/
     &       ,10X,'NUMBER OF ELEMENTS IN MESH B (numelb) -   ',I7,/
     &       ,10X,'NUMBER OF ELEMENT BLOCKS IN B (nblksb) -  ',I7)
  310 FORMAT (10x,'MESH-B INCOMPATIBLE WITH MESH-A',/
     &       ,10X,'NUMBER OF BLOCKS IN MESH-A (nblksa)',I7,/
     &       ,10X,'NUMBER OF BLOCKS IN MESH-B (nblksb)',I7,/
     &       ,10X,'****THIS IS ONLY A WARNING****')
  320 FORMAT (10X,'CORRESPONDS TO MESH-A ELEMENT BLOCK ID',I7,/)
  330 FORMAT (10X,'WORKING ON MESH-B ELEMENT BLOCK       ',I7,/
     &       ,10X,'ELEMENT BLOCK ID                      ',I7)
  410 FORMAT (/,5X,'MESH-B NODE NUMBER ',I7,/
     &       ,5x,' ELEMENT BLOCK     ',I7,/
     &       ,5X,' WAS NOT FOUND IN MESH-A BY SRCHS')
  420 FORMAT (/,5X,'MESH-B CENTROID ELEMENT NUMBER ',I7,/
     &       ,5X,' OF ELEMENT BLOCK                    ',I7,/
     &       ,5X,' WAS NOT FOUND IN MESH-A BY SRCHS')
  430 FORMAT (/,5X,'MESH-B NODE NUMBER ',I7,/
     &       ,5x,' ELEMENT BLOCK     ',I7,/
     &       ,5X,' WAS NOT FOUND IN MESH-A BY SRCHQ')
  440 FORMAT (/,5X,'MESH-B CENTROID ELEMENT NUMBER ',I7,/
     &       ,5X,' OF ELEMENT BLOCK                    ',I7,/
     &       ,5X,' WAS NOT FOUND IN MESH-A BY SRCHQ')
  450 FORMAT (/,5X,'MESH-B NODE NUMBER ',I7,/
     &       ,5x,' ELEMENT BLOCK     ',I7,/
     &       ,5X,' WAS NOT FOUND IN MESH-A BY SRCHH')
  460 FORMAT (/,5X,'MESH-B CENTROID ELEMENT NUMBER ',I7,/
     &       ,5X,' OF ELEMENT BLOCK                    ',I7,/
     &       ,5X,' WAS NOT FOUND IN MESH-A BY SRCHH')
 1001 format(5x,'into  OPNFIL')
 1003 format(5x,'into  RDINPT')
 1005 format(5x,'into  RDA1')
 1009 format(5x,'into  RDB1')
 1010 format(5x,'into  TRUTBL')
 1011 format(5x,'into  RDA2')
 1013 format(5x,'into  RDB2')
 1014 format(5x,'into  BLDSRF')
 1015 format(5x,'into  BLDPTN')
 1017 format(5x,'into  SRCHS - nodes')
 1019 format(5x,'into  SINTPN')
 1021 format(5x,'into  BLDPTE')
 1023 format(5x,'into  SRCHS - element centroids')
 1025 format(5x,'into  SETON0')
 1027 format(5x,'into  SINTPE')
 1028 format(5x,'into  INVCON')
 1029 format(5x,'into  SETON1')
 1031 format(5x,'into  STRAN')
 1035 format(5x,'into  SRCHQ - nodes')
 1037 format(5x,'into  INTRPN')
 1039 format(5x,'into  SRCHQ - element centroids')
 1041 format(5x,'into  ELTON0')      
 1043 format(5x,'into  INTRPE')
 1045 format(5x,'into  ELTON1')
 1047 format(5x,'into  TRANAB')
 1049 format(5x,'into  SRCHH - nodes')
 1050 format(5x,'into  SRCHT - nodes')
 1051 format(5x,'into  INTRPN')
 1053 format(5x,'into  SRCHH - element centroids')
 1054 format(5x,'into  SRCHT - element centroids')
 1055 format(5x,'into  ELTON0')      
 1061 format(5x,'into  WRTC')
 1063 format(5x,'into  CLSFIL',//,5x,'NORMAL TERMINATION')
 1065 format(5x,'into  ELGRAD')
 1067 format(5x,'into  INTRP3')
      END
