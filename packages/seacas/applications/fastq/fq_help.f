C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE HELP_FQ (I)
C************************************************************************

C  SUBROUTINE HELP = WRITES HELP MESSAGES ONTO THE SCREEN

C************************************************************************

C  SUBROUTINE CALLED BY ANY ROUTINE NEEDED HELP MESSAGES

C************************************************************************

C  VARIABLES USED:
C     I = THE POINTER TO DESIGNATE WHICH MESSAGEIS NEEDED

C***********************************************************************

      IF (I .EQ. 1) THEN
         CALL MESSAGE('        ')
         CALL MESSAGE('THE FOLLOWING MAIN OPTIONS ARE AVAILABLE:')
         CALL MESSAGE('     R*EAD       = READS FASTQ DATA FILE')
         CALL MESSAGE('     W*RITE      = WRITES FASTQ DATA FILE')
         CALL MESSAGE('     RW*RITE     = WRITES SPECIFIED REGIONS '//
     &      'TO FASTQ DATA FILE')
         CALL MESSAGE('     BW*RITE     = WRITES SPECIFIED BARSETS '//
     &      'TO FASTQ DATA FILE')
         CALL MESSAGE('     T*ABLET     = GEOMETRY INPUT FROM A '//
     &      'DIGITIZING TABLET')
         CALL MESSAGE('     S*TRAIGHTEN = STRAIGHTENS LINES IN X OR '//
     &      'Y DIRECTION')
         CALL MESSAGE('     K*EY-IN     = INPUTS GEOMETRY FROM '//
     &      'KEYBOARD')
         CALL MESSAGE('     G*RAPHICS   = DRAWS CURRENT FASTQ DATA')
         CALL MESSAGE('     L*IST       = LISTS FASTQ DATA')
         CALL MESSAGE('     M*ESH       = GENERATES THE MESH')
         CALL MESSAGE('     D*ELETE     = DELETES PORTIONS OF '//
     &      'CURRENT FASTQ DATA')
         CALL MESSAGE('     F*LUSH      = CLEARS ALL FASTQ DATA')
         CALL MESSAGE('     EX*IT       = EXITS FASTQ')
         CALL MESSAGE('     SP*AWN      = SPAWNS A SUBPROCESS')

      ELSE IF (I .EQ. 2) THEN
         CALL MESSAGE(' ')
         CALL MESSAGE( '|-----------------+-----------------+'//
     &      '-----------------+-----------------+')
         CALL MESSAGE( '| 0    CLOSEST    | 1     POINT     |'//
     &      ' 2 STRAIGHT LINE | 3   CURSOR      |')
         CALL MESSAGE( '|      (PREFIX)   |0- CLOSEST POINT |'//
     &      '0- CLOSEST LINE  |                 |')
         CALL MESSAGE( '|                 |C- POINT AT POINT|'//
     &      'C- LINE ON LINE  |                 |')
         CALL MESSAGE( '|                 |D- DELETE POINT  |'//
     &      'D- DELETE LINE   |                 |')
         CALL MESSAGE( '+-----------------+-----------------+'//
     &      '-----------------+-----------------+')
         CALL MESSAGE( '| 4  BISECT LINE  | 5    CCW ARC    |'//
     &      ' 6    CW ARC     | 7   REGION      |')
         CALL MESSAGE( '|                 |0- CLOSEST CCW   |'//
     &      '0- CLOSEST CCW   |                 |')
         CALL MESSAGE( '|                 |C- CCW ARC ON ARC|'//
     &      'C- CW ARC ON ARC |                 |')
         CALL MESSAGE( '|                 |                 |'//
     &      '                 |                 |')
         CALL MESSAGE( '+-----------------+-----------------+'//
     &      '-----------------+-----------------+')
         CALL MESSAGE( '| 8  MOVE POINT   | 9 REPAINT SCREEN|'//
     &      ' A TOGGLE GRID   | B ZOOM (BOX)    |')
         CALL MESSAGE( '|                 |                 |'//
     &      '                 |0- PREVIOUS ZOOM |')
         CALL MESSAGE( '|                 |                 |'//
     &      'D- DELETE GRID   |D- RESET ZOOM    |')
         CALL MESSAGE( '|                 |                 |'//
     &      '                 |                 |')
         CALL MESSAGE( '+-----------------+-----------------+'//
     &      '-----------------+-----------------+')
         CALL MESSAGE( '| C               | D               |'//
     &      ' E               | F               |')
         CALL MESSAGE( '|  SLIDE-LINE     |  DELETE         |'//
     &      '  EXIT           |  (NOT USED)     |')
         CALL MESSAGE( '|   (PREFIX)      |   (PREFIX)      |'//
     &      '                 |                 |')
         CALL MESSAGE( '|                 |                 |'//
     &      '                 |                 |')
         CALL MESSAGE( '+-----------------+-----------------+'//
     &      '-----------------+-----------------+')

      ELSE IF (I .EQ. 3) THEN
         CALL MESSAGE(' ')
         CALL MESSAGE('THE FOLLOWING SCHEME AND STEP PROCESSING '//
     &      'CONTROL OPTIONS ARE AVAILABLE:')
         CALL MESSAGE('    +A OR A = INCREASES SMOOTHING PARAMETER '//
     &      'TO FORCE')
         CALL MESSAGE('              MORE EQUAL ELEMENT AREAS '//
     &      '(DEFAULT = 0.7)')
         CALL MESSAGE('    -A      = DECREASES ABOVE SMOOTHING '//
     &      'PARAMETER')
         CALL MESSAGE('     D      = DELETES THE ELEMENT WITH THE '//
     &      'SMALLEST')
         CALL MESSAGE('              INTERIOR ANGLE (MUST BE < '//
     &      'CURRENT MINIMUM ANGLE')
         CALL MESSAGE('              TO BE DELETED - SEE "V"')
         CALL MESSAGE('     E      = EXIT STEP PROCESSING SAVING '//
     &      'REGION')
         CALL MESSAGE('    +F OR F = FREES (INCREASES) SMOOTHING '//
     &      'RELAXATION PARAMETER')
         CALL MESSAGE('              (DEFAULT IS 1.)')
         CALL MESSAGE('    -F      = DECREASES SMOOTHING RELAXATION '//
     &      'PARAMETER')
         CALL MESSAGE('    +I OR I = INCREASES MAXIMUM SMOOTHING '//
     &      'ITERATIONS BY 50%')
         CALL MESSAGE('              (DEFAULT IS 5 * NO. OF ELEMENTS)')
         CALL MESSAGE('    -I      = DECREASES MAXIMUM SMOOTHING '//
     &      'ITERATIONS BY 33%')
         CALL MESSAGE('    +J OR J = INCREASES SMOOTHING NODE '//
     &      'MOVEMENT TOLERANCE')
         CALL MESSAGE('              BY 2**.333 (DEFAULT IS 3% OF '//
     &      'ELEMENT SIDE LENGTH)')
         CALL MESSAGE('    -J      = DECREASES SMOOTHING NODE '//
     &      'MOVEMENT TOLERANCE')
         CALL MESSAGE('              BY 2**.333')
         CALL MESSAGE('     L      = INSERT AN INNER NECKLACES OF '//
     &      'ELEMENTS AROUND A HOLE')
         CALL MESSAGE('     N      = NECKLACES THE REGION WITH A '//
     &      'NEW ROW OF ELEMENTS')
         CALL MESSAGE('     O      = RETURNS PROCESSING TO ORIGINAL '//
     &      'MESH')
         CALL MESSAGE('     P      = PLOTS THE CURRENT MESH')
         CALL MESSAGE('     Q      = ENDS PROCESSING WITHOUT SAVING '//
     &      'MESH (QUIT)')
         CALL MESSAGE('     R      = TRIES A RESTRUCTURE OF THE MESH')
         CALL MESSAGE('     S      = SMOOTHS THE MESH IF SOME '//
     &      'ACTIVITY HAS CHANGED IT')
         CALL MESSAGE('    +V OR V = INCREASES MAXIMUM ANGLE FOR '//
     &      'DELETING ELEMENTS')
         CALL MESSAGE('              (DEFAULT IS SET TO 45 DEGREES)')
         CALL MESSAGE('    -V      = DECREASES MAXIMUM ANGLE FOR '//
     &      'DELETING ELEMENTS')
         CALL MESSAGE('     W      = ATTEMPT RESTRUCTURE OF WORST '//
     &      'ELEMENT ONLY')
         CALL MESSAGE('    +Y OR Y = INCREASES WEIGHT FOR ISOPARAMETRIC'
     &      //' SMOOTHING')
         CALL MESSAGE('              (DEFAULT IS SET TO 0.7)')
         CALL MESSAGE('    -Y      = DECREASES WEIGHT FOR ISOPARAMETRIC'
     &      //' SMOOTHING')
         CALL MESSAGE('     1      = SETS SMOOTHING TO EQUIPOTENTIAL '//
     &      'IF MESH IS')
         CALL MESSAGE('              STRUCTURALLY THE SAME AS THE '//
     &      'ORIGINAL')
         CALL MESSAGE('              (OTHERWISE #2 IS USED)')
         CALL MESSAGE('     2      = SETS SMOOTHING TO AREA PULL AND '//
     &      'LAPLACIAN')
         CALL MESSAGE('     3      = SETS SMOOTHING TO CENTROID '//
     &      'INVERSE AREA PUSH')
         CALL MESSAGE('              AND LAPLACIAN')
         CALL MESSAGE('     4      = SETS SMOOTHING TO CENTROID AREA '//
     &      'PULL')
         CALL MESSAGE('     5      = SETS SMOOTHING TO LAPLACIAN')
         CALL MESSAGE('     6      = SETS SMOOTHING TO LENGTH-'//
     &      'WEIGHTED LAPLACIAN')
         CALL MESSAGE('     7      = SETS SMOOTHING TO WEIGHTED '//
     &      'ISOPARAMETRIC AND LAPLACIAN')
         CALL MESSAGE('     (      = MARKS THE BEGINNING OF A LOOP')
         CALL MESSAGE('     )      = MARKS THE ENDING OF A LOOP')
         CALL MESSAGE('              (LOOP IS DONE WHEN NO ACTIVITY '//
     &      'HAS OCCURRED')
         CALL MESSAGE('              SUCH AS A SMOOTH, DELETION, '//
     &      'RESTRUCTURE, ETC.)')

      ELSE IF (I .EQ. 4) THEN
         CALL MESSAGE(' ')
         CALL MESSAGE('THE FOLLOWING LIST OPTIONS ARE AVAILABLE:')
         CALL MESSAGE('    P*OINT      = LISTS POINT DATA')
         CALL MESSAGE('    L*INE       = LISTS LINE DATA')
         CALL MESSAGE('    SI*DE       = LISTS SIDE DATA')
         CALL MESSAGE('    BA*R SETS   = LISTS BAR SET DATA')
         CALL MESSAGE('    R*EGION     = LISTS REGION DATA')
         CALL MESSAGE('    HO*LE       = LISTS REGION''S HOLE DATA')
         CALL MESSAGE('    S*CHEME     = LISTS SCHEME DATA')
         CALL MESSAGE('    BOD*Y       = LISTS REGIONS IN THE BODY')
         CALL MESSAGE('    B*OUNDARY   = LISTS BOUNDARY FLAGS')
         CALL MESSAGE('    REN*UM      = LISTS RENUM CARDS')
         CALL MESSAGE('    T*HREE      = LISTS 3 NODE GENERATION '//
     &      'TOGGLE')
         CALL MESSAGE('    EI*GHT      = LISTS 8 NODE GENERATION '//
     &      'TOGGLE')
         CALL MESSAGE('    N*INE       = LISTS 9 NODE GENERATION '//
     &      'TOGGLE')
         CALL MESSAGE('    SP*AWN      = SPAWNS A SUBPROCESS')
         CALL MESSAGE('    EX*IT       = EXITS FASTQ')
         CALL MESSAGE('                  (CARRIAGE RETURN TO EXIT '//
     &      'LISTING)')

      ELSE IF (I .EQ. 5) THEN
         CALL MESSAGE(' ')
         CALL MESSAGE('THE FOLLOWING GRAPHICS OPTIONS ARE AVAILABLE:')
         CALL MESSAGE('    P*LOT       = PLOTS THE CURRENT DATA')
         CALL MESSAGE('    R*PLOT      = PLOTS THE DATA BY REGION')
         CALL MESSAGE('    B*PLOT      = PLOTS THE DATA BY BARSET')
         CALL MESSAGE('    SP*LOT      = PLOTS THE DATA BY SIDE')
         CALL MESSAGE('    LP*LOT      = PLOTS THE DATA BY LINE')
         CALL MESSAGE('    A*XIS       = TOGGLES DRAWING OF X-Y AXIS')
         CALL MESSAGE('    PO*INT      = TOGGLES DISPLAY OF POINT NO.S')
         CALL MESSAGE('    L*INE       = TOGGLES DISPLAY OF LINE NO.S')
         CALL MESSAGE('    RE*GION     = TOGGLES DISPLAY OF REGION '//
     &      'NO.S')
         CALL MESSAGE('    F*ULL       = TOGGLES FULL PROPERTY DISPLAY')
         CALL MESSAGE('    PB*OUNDARY  = TOGGLES DISPLAY OF POINBC '//
     &      'FLAGS')
         CALL MESSAGE('    N*BOUNDARY  = TOGGLES DISPLAY OF NODEBC '//
     &      'FLAGS')
         CALL MESSAGE('    EB*OUNDARY  = TOGGLES DISPLAY OF ELEMBC '//
     &      'FLAGS')
         CALL MESSAGE('    I*NTERVAL   = TOGGLES DISPLAY OF '//
     &      'INTERVALS ON LINES')
         CALL MESSAGE('    FA*CTOR     = TOGGLES DISPLAY OF FACTORS '//
     &      'ON LINES')
         CALL MESSAGE('    M*ATERIAL   = TOGGLES DISPLAY OF REGION '//
     &      'BLOCK ID NO.S')
         CALL MESSAGE('    SC*HEME     = TOGGLES DISPLAY OF REGION '//
     &      'SCHEMES')
         CALL MESSAGE('    S*TATUS     = DISPLAYS STATUS OF ALL '//
     &      'TOGGLES')
         CALL MESSAGE('    Z*OOM       = SETS UP NEW PLOT LIMITS')
         CALL MESSAGE('    H*ARDCOPY   = HARDCOPY PLOT FILE OUTPUT')
         CALL MESSAGE('    II*NTERVAL  = INPUTS LINE INTERVALS')
         CALL MESSAGE('    IF*ACTOR    = INPUTS LINE FACTORS')
         CALL MESSAGE('    SPA*WN      = SPAWNS A SUBPROCESS')
         CALL MESSAGE('    EX*IT       = EXITS FASTQ')
         CALL MESSAGE('                  (CARRIAGE RETURN TO EXIT '//
     &      'GRAPHICS)')

      ELSE IF (I .EQ. 6) THEN
         CALL MESSAGE(' ')
         CALL MESSAGE('THE FOLLOWING DELETE OPTIONS ARE AVAILABLE:')
         CALL MESSAGE('    P*OINT    = DELETES POINT DATA')
         CALL MESSAGE('    L*INE     = DELETES LINE DATA')
         CALL MESSAGE('    S*IDE     = DELETES SIDE DATA')
         CALL MESSAGE('    R*EGION   = DELETES REGION DATA')
         CALL MESSAGE('    BA*RSET   = DELETES BARSET DATA')
         CALL MESSAGE('    SC*HEME   = DELETES SCHEME DATA')
         CALL MESSAGE('    B*OUNDARY = DELETES BOUNDARY DATA')
         CALL MESSAGE('    REN*UM    = DELETES RENUMBERING CARDS')
         CALL MESSAGE('    SP*AWN    = SPAWNS A SUBPROCESS')
         CALL MESSAGE('    EX*IT     = EXITS FASTQ')
         CALL MESSAGE('                (CARRIAGE RETURN TO EXIT '//
     &      'DELETE)')

      ELSE IF (I .EQ. 7) THEN
         CALL MESSAGE(' ')
         CALL MESSAGE('THE FOLLOWING DEVICES ARE AVAILABLE:')
         CALL MESSAGE('    AED = AED 512')
         CALL MESSAGE('    ALP = ALPHA NUMERIC TERMINAL')
         CALL MESSAGE('    LS5 = LEAR SIEGLER 220/230 (ENVISION)')
         CALL MESSAGE('    MET = METAFILE')
         CALL MESSAGE('    TEK = TEKTRONICS 4010')
         CALL MESSAGE('    TK4 = TEKTRONICS 4014')
         CALL MESSAGE('    T05 = TEKTRONICS 4105')
         CALL MESSAGE('    T07 = TEKTRONICS 4107, 4109, 4207, 4208')
         CALL MESSAGE('    T13 = TEKTRONICS 4113')
         CALL MESSAGE('    T15 = TEKTRONICS 4115')
         CALL MESSAGE('    V25 = VT 125')
         CALL MESSAGE('    V40 = VT 240')
         CALL MESSAGE('    R25 = RASTER TECH ONE-25')
         CALL MESSAGE('    RET = RETROGRAPHICS')

      ELSE IF (I .EQ. 8) THEN
         CALL MESSAGE(' ')
         CALL MESSAGE('THE FOLLOWING CORRECTION OPTIONS ARE AVAILABLE')
         CALL MESSAGE('    X     = CONSTANT X VALUES ALONG LINE(S)')
         CALL MESSAGE('    Y     = CONSTANT Y VALUES ALONG LINE(S)')
         CALL MESSAGE('    Z*ERO = ZERO X VALUES (CENTERLINES)')

      ELSE IF (I .EQ. 9) THEN
         CALL MESSAGE(' ')
         CALL MESSAGE('THE FOLLOWING KEYIN OPTIONS ARE AVAILABLE:')
         CALL MESSAGE('    T*ITLE    = ENTERS TITLE')
         CALL MESSAGE('    P*OINT    = ENTERS POINT DATA')
         CALL MESSAGE('    L*INE     = ENTERS LINE DATA')
         CALL MESSAGE('    SI*DE     = ENTERS SIDE DATA')
         CALL MESSAGE('    R*EGION   = ENTERS REGION DATA')
         CALL MESSAGE('    HO*LE     = ENTERS REGIONS''S HOLE DATA')
         CALL MESSAGE('    BA*R SET  = ENTERS BARSET DATA')
         CALL MESSAGE('    S*CHEME   = ENTERS SCHEME DATA')
         CALL MESSAGE('    BOD*Y     = ENTERS THE BODY LIST')
         CALL MESSAGE('    I*NTERVAL = ENTERS LINE INTERVALS')
         CALL MESSAGE('    F*ACTOR   = ENTERS LINE FACTORS')
         CALL MESSAGE('    B*OUNDARY = ENTERS LINE/POINT BOUNDARY DATA')
         CALL MESSAGE('    W*EIGHT   = ENTERS A BOUNDARY FLAG '//
     &      'WEIGHTING')
         CALL MESSAGE('    M*ATERIAL = ENTERS REGION MATERIAL NUMBERS')
         CALL MESSAGE('    REN*UM    = ENTERS RENUMBERING CARDS')
         CALL MESSAGE('    O*PTIMIZE = TOGGLES NUMBERING OPTIMIZATION')
         CALL MESSAGE('    TH*REE    = TOGGLES 3 NODE BAR GENERATION')
         CALL MESSAGE('    EI*GHT    = TOGGLES 8 NODE QUAD GENERATION')
         CALL MESSAGE('    N*INE     = TOGGLES 9 NODE QUAD GENERATION')
         CALL MESSAGE('    SP*AWN    = SPAWNS A SUBPROCESS')
         CALL MESSAGE('    EX*IT     = EXITS FASTQ')
         CALL MESSAGE('                (CARRIAGE RETURN TO EXIT KEYIN)')

      ELSE IF (I .EQ. 10) THEN
         CALL MESSAGE(' ')
         CALL MESSAGE('THE FOLLOWING NUMBERING OPTIONS ARE AVAILABLE:')
         CALL MESSAGE('    P*-L-P      = ENTERS POINT-LINE-POINT '//
     &      'SEQUENCE')
         CALL MESSAGE('    X*-Y        = ENTERS X-Y LOCATION TO '//
     &      'START FROM')
         CALL MESSAGE('    N*ODE       = ENTERS NODE NUID''S '//
     &      'LOCATION TO START FROM')

      ELSE IF (I .EQ. 11) THEN
         CALL MESSAGE(' ')
         CALL MESSAGE('THE FOLLOWING MESH GRAPHICS OPTIONS ARE '//
     &      'AVAILABLE:')
         CALL MESSAGE('    A*XIS       = TOGGLES DRAWING OF X-Y AXIS')
         CALL MESSAGE('    N*NUMBERING = TOGGLES DISPLAY OF NODE '//
     &      'NUMBERS')
         CALL MESSAGE('    EN*UMBERING = TOGGLES DISPLAY OF ELEMENT '//
     &      'NUMBERS')
         CALL MESSAGE('    MN*UMBERING = TOGGLES DISPLAY OF BLOCK ID '//
     &      'NUMBERS')
         CALL MESSAGE('    O*ORDER     = TOGGLES DISPLAY OF '//
     &      'OPTIMIZED ELEMENT ORDER')
         CALL MESSAGE('    NB*OUNDARY  = TOGGLES DISPLAY OF NODAL '//
     &      'BOUNDARIES')
         CALL MESSAGE('    EB*OUNDARY  = TOGGLES DISPLAY OF ELEMENT '//
     &      'BOUNDARIES')
         CALL MESSAGE('    W*EIGHT     = TOGGLES DISPLAY OF '//
     &      'WEIGHTING FACTORS')
         CALL MESSAGE('    S*TATUS     = DISPLAYS STATUS OF ALL '//
     &      'TOGGLES')
         CALL MESSAGE('    P*LOT       = PLOTS THE MESH USING '//
     &      'CURRENT ZOOM')
         CALL MESSAGE('    EP*LOT      = PLOTS THE MESH BASED ON '//
     &      'ELEMENT NUMBERS')
         CALL MESSAGE('    R*PLOT      = PLOTS THE MESH BASED ON '//
     &      'REGION NUMBERS')
         CALL MESSAGE('    B*PLOT      = PLOTS THE MESH BASED ON '//
     &      'BARSET NUMBERS')
         CALL MESSAGE('    M*PLOT      = PLOTS THE MESH BASED ON '//
     &      'MATERIAL NUMBERS')
         CALL MESSAGE('    Z*OOM       = SETS UP NEW PLOT LIMITS')
         CALL MESSAGE('    H*ARDCOPY   = HARDCOPY PLOT FILE OUTPUT')
         CALL MESSAGE('    SP*AWN      = SPAWNS A SUBPROCESS')
         CALL MESSAGE('    EX*IT       = EXITS FASTQ')
         CALL MESSAGE('                  (CARRIAGE RETURN TO EXIT '//
     &      'MESH GRAPHICS)')

      ELSE IF (I .EQ. 12) THEN
         CALL MESSAGE(' ')
         CALL MESSAGE('THE FOLLOWING MESH OPTIONS ARE AVAILABLE:')
         CALL MESSAGE('    P*ROCESS    = PROCESSES THE MESH')
         CALL MESSAGE('    S*TEP       = STEPS THROUGH PROCESSING '//
     &      'INTERACTIVELY')
         CALL MESSAGE('    G*RAPHICS   = GRAPHICALLY DISPLAYS MESH')
         CALL MESSAGE('    I*NTERVALS  = ENTERS LINE INTERVALS')
         CALL MESSAGE('    F*ACTOR     = ENTERS LINE FACTORS')
         CALL MESSAGE('    SI*ZE       = ENTERS SIZES FOR REGIONS')
         CALL MESSAGE('    O*PTIMIZE   = TOGGLE FOR BANDWIDTH '//
     &      'OPTIMIZATION')
         CALL MESSAGE('    T*HREE      = TOGGLE FOR 3 NODE BAR '//
     &      'GENERATION')
         CALL MESSAGE('    EI*GHT      = TOGGLE FOR 8 NODE QUAD '//
     &      'GENERATION')
         CALL MESSAGE('    NI*NE       = TOGGLE FOR 9 NODE QUAD '//
     &      'GENERATION')
         CALL MESSAGE('    R*EAD       = READS MESH DATA FROM A FILE')
         CALL MESSAGE('    R*MESH      = REMESHES BASED ON AN ERROR'//
     &      ' ESTIMATE')
         CALL MESSAGE('    AD*JUST     = ADJUSTS A GENERATED MESH')
         CALL MESSAGE('    W*RITE      = WRITES A GENESIS MESH FILE')
         CALL MESSAGE('    A*BAQUS     = WRITES AN ABAQUS MESH FILE')
         CALL MESSAGE('    N*ASTRAN    = WRITES A NASTRAN MESH FILE')
         CALL MESSAGE('    SP*AWN      = SPAWNS A SUBPROCESS')
         CALL MESSAGE('    EX*IT       = EXITS FASTQ')
         CALL MESSAGE('    D*ISTORTION = CALCULATES A DISTORTION '//
     &      'INDEX')
         CALL MESSAGE('                  (CARRIAGE RETURN TO EXIT '//
     &      'MESH)')

      ELSE IF (I .EQ. 13) THEN
         CALL MESSAGE('THE FOLLOWING INITIAL MESH GENERATION SCHEMES '//
     &      'ARE AVAILABLE:')
         CALL MESSAGE('    B      = TRANSITION REGION GENERATION')
         CALL MESSAGE('    C      = SEMICIRCULAR REGION GENERATION')
         CALL MESSAGE('    M      = AUTOMATED RECTANGULAR MESH '//
     &      'GENERATION')
         CALL MESSAGE('    T      = TRIANGULAR REGION GENERATION')
         CALL MESSAGE('    U      = PENATGON REGION GENERATION')
         CALL MESSAGE('    X      = PAVING MESH GENERATION')
         CALL MESSAGE('    Z      = HOLE REGION GENERATION')
         CALL MESSAGE('             NO SCHEME (CARRIAGE RETURN) '//
     &      'DEFAULTS TO A')
         CALL MESSAGE('             FORCED RECTANGULAR SCHEME')

      ELSE IF (I .EQ. 14) THEN
         CALL MESSAGE('THE FOLLOWING TABLET OPTIONS ARE '//
     &      'AVAILABLE: ')
         CALL MESSAGE('    A*XIS       = TOGGLES DRAWING OF X-Y AXIS')
         CALL MESSAGE('    B*UTTONS    = SHOWS THE MOUSE BUTTON'//
     &      ' DEFINITIONS')
         CALL MESSAGE('    C*LEAR GRID = ERASES ALL GRID LINES')
         CALL MESSAGE('    XC*LEAR     = ERASES ALL X GRID LINES')
         CALL MESSAGE('    YC*LEAR     = ERASES ALL Y GRID LINES')
         CALL MESSAGE('    D*IGITIZE   = DIGITZES GEOMETRY WITH MOUSE'//
     &      '/TABLET')
         CALL MESSAGE('    DE*FAULT    = USES THE ZOOM LIMITS TO'//
     &      ' INITIALIZE TABLET')
         CALL MESSAGE('    I*NITIALIZE = INITIALIZES A DRAWING TO THE'//
     &      ' TABLET')
         CALL MESSAGE('    S*NAP       = TOGGLES THE SNAP MODE FOR'//
     &      ' GRID INTERSECTIONS')
         CALL MESSAGE('    U*NIFORM    = ADDS UNIFORMLY SPACED SQUARE'//
     &      ' GRID LINES')
         CALL MESSAGE('    UX          = ADDS UNIFORMLY X SPACED GRID'//
     &      ' LINES')
         CALL MESSAGE('    UY          = ADDS UNIFORMLY Y SPACED GRID'//
     &      ' LINES')
         CALL MESSAGE('    X*GRID      = ADDS ARBITRARY X GRID LINES')
         CALL MESSAGE('    Y*GRID      = ADDS ARBITRARY Y GRID LINES')
         CALL MESSAGE('    P*OINT GRID = ADDS X AND Y GRIDS THROUGH'//
     &      ' ALL EXISTING POINTS')
         CALL MESSAGE('    Z*OOM       = SETS PLOTTING (AND '//
     &      ' TABLET) LIMITS')
         CALL MESSAGE('    SP*AWN      = SPAWNS A SUBPROCESS')
         CALL MESSAGE('    EX*IT       = EXITS FASTQ')
         CALL MESSAGE('                  (CARRIAGE RETURN TO EXIT '//
     &      'TABLET)')
      END IF

      RETURN

      END
