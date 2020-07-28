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
C     I = THE POINTER TO DESIGNATE WHICH MESSAGE IS NEEDED

C***********************************************************************

      IF (I .EQ. 1) THEN
         CALL MESAGE ('        ')
         CALL MESAGE ('THE FOLLOWING MAIN OPTIONS ARE AVAILABLE:')
         CALL MESAGE ('     R*EAD       = READS FASTQ DATA FILE')
         CALL MESAGE ('     W*RITE      = WRITES FASTQ DATA FILE')
         CALL MESAGE ('     RW*RITE     = WRITES SPECIFIED REGIONS '//
     &      'TO FASTQ DATA FILE')
         CALL MESAGE ('     BW*RITE     = WRITES SPECIFIED BARSETS '//
     &      'TO FASTQ DATA FILE')
         CALL MESAGE ('     T*ABLET     = GEOMETRY INPUT FROM A '//
     &      'DIGITIZING TABLET')
         CALL MESAGE ('     S*TRAIGHTEN = STRAIGHTENS LINES IN X OR '//
     &      'Y DIRECTION')
         CALL MESAGE ('     K*EY-IN     = INPUTS GEOMETRY FROM '//
     &      'KEYBOARD')
         CALL MESAGE ('     G*RAPHICS   = DRAWS CURRENT FASTQ DATA')
         CALL MESAGE ('     L*IST       = LISTS FASTQ DATA')
         CALL MESAGE ('     M*ESH       = GENERATES THE MESH')
         CALL MESAGE ('     D*ELETE     = DELETES PORTIONS OF '//
     &      'CURRENT FASTQ DATA')
         CALL MESAGE ('     F*LUSH      = CLEARS ALL FASTQ DATA')
         CALL MESAGE ('     EX*IT       = EXITS FASTQ')
         CALL MESAGE ('     SP*AWN      = SPAWNS A SUBPROCESS')

      ELSE IF (I .EQ. 2) THEN
         CALL MESAGE (' ')
         CALL MESAGE ( '|-----------------+-----------------+'//
     &      '-----------------+-----------------+')
         CALL MESAGE ( '| 0    CLOSEST    | 1     POINT     |'//
     &      ' 2 STRAIGHT LINE | 3   CURSOR      |')
         CALL MESAGE ( '|      (PREFIX)   |0- CLOSEST POINT |'//
     &      '0- CLOSEST LINE  |                 |')
         CALL MESAGE ( '|                 |C- POINT AT POINT|'//
     &      'C- LINE ON LINE  |                 |')
         CALL MESAGE ( '|                 |D- DELETE POINT  |'//
     &      'D- DELETE LINE   |                 |')
         CALL MESAGE ( '+-----------------+-----------------+'//
     &      '-----------------+-----------------+')
         CALL MESAGE ( '| 4  BISECT LINE  | 5    CCW ARC    |'//
     &      ' 6    CW ARC     | 7   REGION      |')
         CALL MESAGE ( '|                 |0- CLOSEST CCW   |'//
     &      '0- CLOSEST CCW   |                 |')
         CALL MESAGE ( '|                 |C- CCW ARC ON ARC|'//
     &      'C- CW ARC ON ARC |                 |')
         CALL MESAGE ( '|                 |                 |'//
     &      '                 |                 |')
         CALL MESAGE ( '+-----------------+-----------------+'//
     &      '-----------------+-----------------+')
         CALL MESAGE ( '| 8  MOVE POINT   | 9 REPAINT SCREEN|'//
     &      ' A TOGGLE GRID   | B ZOOM (BOX)    |')
         CALL MESAGE ( '|                 |                 |'//
     &      '                 |0- PREVIOUS ZOOM |')
         CALL MESAGE ( '|                 |                 |'//
     &      'D- DELETE GRID   |D- RESET ZOOM    |')
         CALL MESAGE ( '|                 |                 |'//
     &      '                 |                 |')
         CALL MESAGE ( '+-----------------+-----------------+'//
     &      '-----------------+-----------------+')
         CALL MESAGE ( '| C               | D               |'//
     &      ' E               | F               |')
         CALL MESAGE ( '|  SLIDE-LINE     |  DELETE         |'//
     &      '  EXIT           |  (NOT USED)     |')
         CALL MESAGE ( '|   (PREFIX)      |   (PREFIX)      |'//
     &      '                 |                 |')
         CALL MESAGE ( '|                 |                 |'//
     &      '                 |                 |')
         CALL MESAGE ( '+-----------------+-----------------+'//
     &      '-----------------+-----------------+')

      ELSE IF (I .EQ. 3) THEN
         CALL MESAGE (' ')
         CALL MESAGE ('THE FOLLOWING SCHEME AND STEP PROCESSING '//
     &      'CONTROL OPTIONS ARE AVAILABLE:')
         CALL MESAGE ('    +A OR A = INCREASES SMOOTHING PARAMETER '//
     &      'TO FORCE')
         CALL MESAGE ('              MORE EQUAL ELEMENT AREAS '//
     &      '(DEFAULT = 0.7)')
         CALL MESAGE ('    -A      = DECREASES ABOVE SMOOTHING '//
     &      'PARAMETER')
         CALL MESAGE ('     D      = DELETES THE ELEMENT WITH THE '//
     &      'SMALLEST')
         CALL MESAGE ('              INTERIOR ANGLE (MUST BE < '//
     &      'CURRENT MINIMUM ANGLE')
         CALL MESAGE ('              TO BE DELETED - SEE "V"')
         CALL MESAGE ('     E      = EXIT STEP PROCESSING SAVING '//
     &      'REGION')
         CALL MESAGE ('    +F OR F = FREES (INCREASES) SMOOTHING '//
     &      'RELAXATION PARAMETER')
         CALL MESAGE ('              (DEFAULT IS 1.)')
         CALL MESAGE ('    -F      = DECREASES SMOOTHING RELAXATION '//
     &      'PARAMETER')
         CALL MESAGE ('    +I OR I = INCREASES MAXIMUM SMOOTHING '//
     &      'ITERATIONS BY 50%')
         CALL MESAGE ('              (DEFAULT IS 5 * NO. OF ELEMENTS)')
         CALL MESAGE ('    -I      = DECREASES MAXIMUM SMOOTHING '//
     &      'ITERATIONS BY 33%')
         CALL MESAGE ('    +J OR J = INCREASES SMOOTHING NODE '//
     &      'MOVEMENT TOLERANCE')
         CALL MESAGE ('              BY 2**.333 (DEFAULT IS 3% OF '//
     &      'ELEMENT SIDE LENGTH)')
         CALL MESAGE ('    -J      = DECREASES SMOOTHING NODE '//
     &      'MOVEMENT TOLERANCE')
         CALL MESAGE ('              BY 2**.333')
         CALL MESAGE ('     L      = INSERT AN INNER NECKLACES OF '//
     &      'ELEMENTS AROUND A HOLE')
         CALL MESAGE ('     N      = NECKLACES THE REGION WITH A '//
     &      'NEW ROW OF ELEMENTS')
         CALL MESAGE ('     O      = RETURNS PROCESSING TO ORIGINAL '//
     &      'MESH')
         CALL MESAGE ('     P      = PLOTS THE CURRENT MESH')
         CALL MESAGE ('     Q      = ENDS PROCESSING WITHOUT SAVING '//
     &      'MESH (QUIT)')
         CALL MESAGE ('     R      = TRIES A RESTRUCTURE OF THE MESH')
         CALL MESAGE ('     S      = SMOOTHS THE MESH IF SOME '//
     &      'ACTIVITY HAS CHANGED IT')
         CALL MESAGE ('    +V OR V = INCREASES MAXIMUM ANGLE FOR '//
     &      'DELETING ELEMENTS')
         CALL MESAGE ('              (DEFAULT IS SET TO 45 DEGREES)')
         CALL MESAGE ('    -V      = DECREASES MAXIMUM ANGLE FOR '//
     &      'DELETING ELEMENTS')
         CALL MESAGE ('     W      = ATTEMPT RESTRUCTURE OF WORST '//
     &      'ELEMENT ONLY')
         CALL MESAGE ('    +Y OR Y = INCREASES WEIGHT FOR ISOPARAMETRIC'
     &      //' SMOOTHING')
         CALL MESAGE ('              (DEFAULT IS SET TO 0.7)')
         CALL MESAGE ('    -Y      = DECREASES WEIGHT FOR ISOPARAMETRIC'
     &      //' SMOOTHING')
         CALL MESAGE ('     1      = SETS SMOOTHING TO EQUIPOTENTIAL '//
     &      'IF MESH IS')
         CALL MESAGE ('              STRUCTURALLY THE SAME AS THE '//
     &      'ORIGINAL')
         CALL MESAGE ('              (OTHERWISE #2 IS USED)')
         CALL MESAGE ('     2      = SETS SMOOTHING TO AREA PULL AND '//
     &      'LAPLACIAN')
         CALL MESAGE ('     3      = SETS SMOOTHING TO CENTROID '//
     &      'INVERSE AREA PUSH')
         CALL MESAGE ('              AND LAPLACIAN')
         CALL MESAGE ('     4      = SETS SMOOTHING TO CENTROID AREA '//
     &      'PULL')
         CALL MESAGE ('     5      = SETS SMOOTHING TO LAPLACIAN')
         CALL MESAGE ('     6      = SETS SMOOTHING TO LENGTH-'//
     &      'WEIGHTED LAPLACIAN')
         CALL MESAGE ('     7      = SETS SMOOTHING TO WEIGHTED '//
     &      'ISOPARAMETRIC AND LAPLACIAN')
         CALL MESAGE ('     (      = MARKS THE BEGINNING OF A LOOP')
         CALL MESAGE ('     )      = MARKS THE ENDING OF A LOOP')
         CALL MESAGE ('              (LOOP IS DONE WHEN NO ACTIVITY '//
     &      'HAS OCCURRED')
         CALL MESAGE ('              SUCH AS A SMOOTH, DELETION, '//
     &      'RESTRUCTURE, ETC.)')

      ELSE IF (I .EQ. 4) THEN
         CALL MESAGE (' ')
         CALL MESAGE ('THE FOLLOWING LIST OPTIONS ARE AVAILABLE:')
         CALL MESAGE ('    P*OINT      = LISTS POINT DATA')
         CALL MESAGE ('    L*INE       = LISTS LINE DATA')
         CALL MESAGE ('    SI*DE       = LISTS SIDE DATA')
         CALL MESAGE ('    BA*R SETS   = LISTS BAR SET DATA')
         CALL MESAGE ('    R*EGION     = LISTS REGION DATA')
         CALL MESAGE ('    HO*LE       = LISTS REGION''S HOLE DATA')
         CALL MESAGE ('    S*CHEME     = LISTS SCHEME DATA')
         CALL MESAGE ('    BOD*Y       = LISTS REGIONS IN THE BODY')
         CALL MESAGE ('    B*OUNDARY   = LISTS BOUNDARY FLAGS')
         CALL MESAGE ('    REN*UM      = LISTS RENUM CARDS')
         CALL MESAGE ('    T*HREE      = LISTS 3 NODE GENERATION '//
     &      'TOGGLE')
         CALL MESAGE ('    EI*GHT      = LISTS 8 NODE GENERATION '//
     &      'TOGGLE')
         CALL MESAGE ('    N*INE       = LISTS 9 NODE GENERATION '//
     &      'TOGGLE')
         CALL MESAGE ('    SP*AWN      = SPAWNS A SUBPROCESS')
         CALL MESAGE ('    EX*IT       = EXITS FASTQ')
         CALL MESAGE ('                  (CARRIAGE RETURN TO EXIT '//
     &      'LISTING)')

      ELSE IF (I .EQ. 5) THEN
         CALL MESAGE (' ')
         CALL MESAGE ('THE FOLLOWING GRAPHICS OPTIONS ARE AVAILABLE:')
         CALL MESAGE ('    P*LOT       = PLOTS THE CURRENT DATA')
         CALL MESAGE ('    R*PLOT      = PLOTS THE DATA BY REGION')
         CALL MESAGE ('    B*PLOT      = PLOTS THE DATA BY BARSET')
         CALL MESAGE ('    SP*LOT      = PLOTS THE DATA BY SIDE')
         CALL MESAGE ('    LP*LOT      = PLOTS THE DATA BY LINE')
         CALL MESAGE ('    A*XIS       = TOGGLES DRAWING OF X-Y AXIS')
         CALL MESAGE ('    PO*INT      = TOGGLES DISPLAY OF POINT NO.S')
         CALL MESAGE ('    L*INE       = TOGGLES DISPLAY OF LINE NO.S')
         CALL MESAGE ('    RE*GION     = TOGGLES DISPLAY OF REGION '//
     &      'NO.S')
         CALL MESAGE ('    F*ULL       = TOGGLES FULL PROPERTY DISPLAY')
         CALL MESAGE ('    PB*OUNDARY  = TOGGLES DISPLAY OF POINBC '//
     &      'FLAGS')
         CALL MESAGE ('    N*BOUNDARY  = TOGGLES DISPLAY OF NODEBC '//
     &      'FLAGS')
         CALL MESAGE ('    EB*OUNDARY  = TOGGLES DISPLAY OF ELEMBC '//
     &      'FLAGS')
         CALL MESAGE ('    I*NTERVAL   = TOGGLES DISPLAY OF '//
     &      'INTERVALS ON LINES')
         CALL MESAGE ('    FA*CTOR     = TOGGLES DISPLAY OF FACTORS '//
     &      'ON LINES')
         CALL MESAGE ('    M*ATERIAL   = TOGGLES DISPLAY OF REGION '//
     &      'BLOCK ID NO.S')
         CALL MESAGE ('    SC*HEME     = TOGGLES DISPLAY OF REGION '//
     &      'SCHEMES')
         CALL MESAGE ('    S*TATUS     = DISPLAYS STATUS OF ALL '//
     &      'TOGGLES')
         CALL MESAGE ('    Z*OOM       = SETS UP NEW PLOT LIMITS')
         CALL MESAGE ('    H*ARDCOPY   = HARDCOPY PLOT FILE OUTPUT')
         CALL MESAGE ('    II*NTERVAL  = INPUTS LINE INTERVALS')
         CALL MESAGE ('    IF*ACTOR    = INPUTS LINE FACTORS')
         CALL MESAGE ('    SPA*WN      = SPAWNS A SUBPROCESS')
         CALL MESAGE ('    EX*IT       = EXITS FASTQ')
         CALL MESAGE ('                  (CARRIAGE RETURN TO EXIT '//
     &      'GRAPHICS)')

      ELSE IF (I .EQ. 6) THEN
         CALL MESAGE (' ')
         CALL MESAGE ('THE FOLLOWING DELETE OPTIONS ARE AVAILABLE:')
         CALL MESAGE ('    P*OINT    = DELETES POINT DATA')
         CALL MESAGE ('    L*INE     = DELETES LINE DATA')
         CALL MESAGE ('    S*IDE     = DELETES SIDE DATA')
         CALL MESAGE ('    R*EGION   = DELETES REGION DATA')
         CALL MESAGE ('    BA*RSET   = DELETES BARSET DATA')
         CALL MESAGE ('    SC*HEME   = DELETES SCHEME DATA')
         CALL MESAGE ('    B*OUNDARY = DELETES BOUNDARY DATA')
         CALL MESAGE ('    REN*UM    = DELETES RENUMBERING CARDS')
         CALL MESAGE ('    SP*AWN    = SPAWNS A SUBPROCESS')
         CALL MESAGE ('    EX*IT     = EXITS FASTQ')
         CALL MESAGE ('                (CARRIAGE RETURN TO EXIT '//
     &      'DELETE)')

      ELSE IF (I .EQ. 7) THEN
         CALL MESAGE (' ')
         CALL MESAGE ('THE FOLLOWING DEVICES ARE AVAILABLE:')
         CALL MESAGE ('    AED = AED 512')
         CALL MESAGE ('    ALP = ALPHA NUMERIC TERMINAL')
         CALL MESAGE ('    LS5 = LEAR SIEGLER 220/230 (ENVISION)')
         CALL MESAGE ('    MET = METAFILE')
         CALL MESAGE ('    TEK = TEKTRONICS 4010')
         CALL MESAGE ('    TK4 = TEKTRONICS 4014')
         CALL MESAGE ('    T05 = TEKTRONICS 4105')
         CALL MESAGE ('    T07 = TEKTRONICS 4107, 4109, 4207, 4208')
         CALL MESAGE ('    T13 = TEKTRONICS 4113')
         CALL MESAGE ('    T15 = TEKTRONICS 4115')
         CALL MESAGE ('    V25 = VT 125')
         CALL MESAGE ('    V40 = VT 240')
         CALL MESAGE ('    R25 = RASTER TECH ONE-25')
         CALL MESAGE ('    RET = RETROGRAPHICS')

      ELSE IF (I .EQ. 8) THEN
         CALL MESAGE (' ')
         CALL MESAGE ('THE FOLLOWING CORRECTION OPTIONS ARE AVAILABLE')
         CALL MESAGE ('    X     = CONSTANT X VALUES ALONG LINE(S)')
         CALL MESAGE ('    Y     = CONSTANT Y VALUES ALONG LINE(S)')
         CALL MESAGE ('    Z*ERO = ZERO X VALUES (CENTERLINES)')

      ELSE IF (I .EQ. 9) THEN
         CALL MESAGE (' ')
         CALL MESAGE ('THE FOLLOWING KEYIN OPTIONS ARE AVAILABLE:')
         CALL MESAGE ('    T*ITLE    = ENTERS TITLE')
         CALL MESAGE ('    P*OINT    = ENTERS POINT DATA')
         CALL MESAGE ('    L*INE     = ENTERS LINE DATA')
         CALL MESAGE ('    SI*DE     = ENTERS SIDE DATA')
         CALL MESAGE ('    R*EGION   = ENTERS REGION DATA')
         CALL MESAGE ('    HO*LE     = ENTERS REGIONS''S HOLE DATA')
         CALL MESAGE ('    BA*R SET  = ENTERS BARSET DATA')
         CALL MESAGE ('    S*CHEME   = ENTERS SCHEME DATA')
         CALL MESAGE ('    BOD*Y     = ENTERS THE BODY LIST')
         CALL MESAGE ('    I*NTERVAL = ENTERS LINE INTERVALS')
         CALL MESAGE ('    F*ACTOR   = ENTERS LINE FACTORS')
         CALL MESAGE ('    B*OUNDARY = ENTERS LINE/POINT BOUNDARY DATA')
         CALL MESAGE ('    W*EIGHT   = ENTERS A BOUNDARY FLAG '//
     &      'WEIGHTING')
         CALL MESAGE ('    M*ATERIAL = ENTERS REGION MATERIAL NUMBERS')
         CALL MESAGE ('    REN*UM    = ENTERS RENUMBERING CARDS')
         CALL MESAGE ('    O*PTIMIZE = TOGGLES NUMBERING OPTIMIZATION')
         CALL MESAGE ('    TH*REE    = TOGGLES 3 NODE BAR GENERATION')
         CALL MESAGE ('    EI*GHT    = TOGGLES 8 NODE QUAD GENERATION')
         CALL MESAGE ('    N*INE     = TOGGLES 9 NODE QUAD GENERATION')
         CALL MESAGE ('    SP*AWN    = SPAWNS A SUBPROCESS')
         CALL MESAGE ('    EX*IT     = EXITS FASTQ')
         CALL MESAGE ('                (CARRIAGE RETURN TO EXIT KEYIN)')

      ELSE IF (I .EQ. 10) THEN
         CALL MESAGE (' ')
         CALL MESAGE ('THE FOLLOWING NUMBERING OPTIONS ARE AVAILABLE:')
         CALL MESAGE ('    P*-L-P      = ENTERS POINT-LINE-POINT '//
     &      'SEQUENCE')
         CALL MESAGE ('    X*-Y        = ENTERS X-Y LOCATION TO '//
     &      'START FROM')
         CALL MESAGE ('    N*ODE       = ENTERS NODE NUID''S '//
     &      'LOCATION TO START FROM')

      ELSE IF (I .EQ. 11) THEN
         CALL MESAGE (' ')
         CALL MESAGE ('THE FOLLOWING MESH GRAPHICS OPTIONS ARE '//
     &      'AVAILABLE:')
         CALL MESAGE ('    A*XIS       = TOGGLES DRAWING OF X-Y AXIS')
         CALL MESAGE ('    N*NUMBERING = TOGGLES DISPLAY OF NODE '//
     &      'NUMBERS')
         CALL MESAGE ('    EN*UMBERING = TOGGLES DISPLAY OF ELEMENT '//
     &      'NUMBERS')
         CALL MESAGE ('    MN*UMBERING = TOGGLES DISPLAY OF BLOCK ID '//
     &      'NUMBERS')
         CALL MESAGE ('    O*ORDER     = TOGGLES DISPLAY OF '//
     &      'OPTIMIZED ELEMENT ORDER')
         CALL MESAGE ('    NB*OUNDARY  = TOGGLES DISPLAY OF NODAL '//
     &      'BOUNDARIES')
         CALL MESAGE ('    EB*OUNDARY  = TOGGLES DISPLAY OF ELEMENT '//
     &      'BOUNDARIES')
         CALL MESAGE ('    W*EIGHT     = TOGGLES DISPLAY OF '//
     &      'WEIGHTING FACTORS')
         CALL MESAGE ('    S*TATUS     = DISPLAYS STATUS OF ALL '//
     &      'TOGGLES')
         CALL MESAGE ('    P*LOT       = PLOTS THE MESH USING '//
     &      'CURRENT ZOOM')
         CALL MESAGE ('    EP*LOT      = PLOTS THE MESH BASED ON '//
     &      'ELEMENT NUMBERS')
         CALL MESAGE ('    R*PLOT      = PLOTS THE MESH BASED ON '//
     &      'REGION NUMBERS')
         CALL MESAGE ('    B*PLOT      = PLOTS THE MESH BASED ON '//
     &      'BARSET NUMBERS')
         CALL MESAGE ('    M*PLOT      = PLOTS THE MESH BASED ON '//
     &      'MATERIAL NUMBERS')
         CALL MESAGE ('    Z*OOM       = SETS UP NEW PLOT LIMITS')
         CALL MESAGE ('    H*ARDCOPY   = HARDCOPY PLOT FILE OUTPUT')
         CALL MESAGE ('    SP*AWN      = SPAWNS A SUBPROCESS')
         CALL MESAGE ('    EX*IT       = EXITS FASTQ')
         CALL MESAGE ('                  (CARRIAGE RETURN TO EXIT '//
     &      'MESH GRAPHICS)')

      ELSE IF (I .EQ. 12) THEN
         CALL MESAGE (' ')
         CALL MESAGE ('THE FOLLOWING MESH OPTIONS ARE AVAILABLE:')
         CALL MESAGE ('    P*ROCESS    = PROCESSES THE MESH')
         CALL MESAGE ('    S*TEP       = STEPS THROUGH PROCESSING '//
     &      'INTERACTIVELY')
         CALL MESAGE ('    G*RAPHICS   = GRAPHICALLY DISPLAYS MESH')
         CALL MESAGE ('    I*NTERVALS  = ENTERS LINE INTERVALS')
         CALL MESAGE ('    F*ACTOR     = ENTERS LINE FACTORS')
         CALL MESAGE ('    SI*ZE       = ENTERS SIZES FOR REGIONS')
         CALL MESAGE ('    O*PTIMIZE   = TOGGLE FOR BANDWIDTH '//
     &      'OPTIMIZATION')
         CALL MESAGE ('    T*HREE      = TOGGLE FOR 3 NODE BAR '//
     &      'GENERATION')
         CALL MESAGE ('    EI*GHT      = TOGGLE FOR 8 NODE QUAD '//
     &      'GENERATION')
         CALL MESAGE ('    NI*NE       = TOGGLE FOR 9 NODE QUAD '//
     &      'GENERATION')
         CALL MESAGE ('    R*EAD       = READS MESH DATA FROM A FILE')
         CALL MESAGE ('    R*MESH      = REMESHES BASED ON AN ERROR'//
     &      ' ESTIMATE')
         CALL MESAGE ('    AD*JUST     = ADJUSTS A GENERATED MESH')
         CALL MESAGE ('    W*RITE      = WRITES A GENESIS MESH FILE')
         CALL MESAGE ('    A*BAQUS     = WRITES AN ABAQUS MESH FILE')
         CALL MESAGE ('    N*ASTRAN    = WRITES A NASTRAN MESH FILE')
         CALL MESAGE ('    SP*AWN      = SPAWNS A SUBPROCESS')
         CALL MESAGE ('    EX*IT       = EXITS FASTQ')
         CALL MESAGE ('    D*ISTORTION = CALCULATES A DISTORTION '//
     &      'INDEX')
         CALL MESAGE ('                  (CARRIAGE RETURN TO EXIT '//
     &      'MESH)')

      ELSE IF (I .EQ. 13) THEN
         CALL MESAGE ('THE FOLLOWING INITIAL MESH GENERATION SCHEMES '//
     &      'ARE AVAILABLE:')
         CALL MESAGE ('    B      = TRANSITION REGION GENERATION')
         CALL MESAGE ('    C      = SEMICIRCULAR REGION GENERATION')
         CALL MESAGE ('    M      = AUTOMATED RECTANGULAR MESH '//
     &      'GENERATION')
         CALL MESAGE ('    T      = TRIANGULAR REGION GENERATION')
         CALL MESAGE ('    U      = PENATGON REGION GENERATION')
         CALL MESAGE ('    X      = PAVING MESH GENERATION')
         CALL MESAGE ('    Z      = HOLE REGION GENERATION')
         CALL MESAGE ('             NO SCHEME (CARRIAGE RETURN) '//
     &      'DEFAULTS TO A')
         CALL MESAGE ('             FORCED RECTANGULAR SCHEME')

      ELSE IF (I .EQ. 14) THEN
         CALL MESAGE ('THE FOLLOWING TABLET OPTIONS ARE '//
     &      'AVAILABLE: ')
         CALL MESAGE ('    A*XIS       = TOGGLES DRAWING OF X-Y AXIS')
         CALL MESAGE ('    B*UTTONS    = SHOWS THE MOUSE BUTTON'//
     &      ' DEFINITIONS')
         CALL MESAGE ('    C*LEAR GRID = ERASES ALL GRID LINES')
         CALL MESAGE ('    XC*LEAR     = ERASES ALL X GRID LINES')
         CALL MESAGE ('    YC*LEAR     = ERASES ALL Y GRID LINES')
         CALL MESAGE ('    D*IGITIZE   = DIGITZES GEOMETRY WITH MOUSE'//
     &      '/TABLET')
         CALL MESAGE ('    DE*FAULT    = USES THE ZOOM LIMITS TO'//
     &      ' INITIALIZE TABLET')
         CALL MESAGE ('    I*NITIALIZE = INITIALIZES A DRAWING TO THE'//
     &      ' TABLET')
         CALL MESAGE ('    S*NAP       = TOGGLES THE SNAP MODE FOR'//
     &      ' GRID INTERSECTIONS')
         CALL MESAGE ('    U*NIFORM    = ADDS UNIFORMLY SPACED SQUARE'//
     &      ' GRID LINES')
         CALL MESAGE ('    UX          = ADDS UNIFORMLY X SPACED GRID'//
     &      ' LINES')
         CALL MESAGE ('    UY          = ADDS UNIFORMLY Y SPACED GRID'//
     &      ' LINES')
         CALL MESAGE ('    X*GRID      = ADDS ARBITRARY X GRID LINES')
         CALL MESAGE ('    Y*GRID      = ADDS ARBITRARY Y GRID LINES')
         CALL MESAGE ('    P*OINT GRID = ADDS X AND Y GRIDS THROUGH'//
     &      ' ALL EXISTING POINTS')
         CALL MESAGE ('    Z*OOM       = SETS PLOTTING (AND '//
     &      ' TABLET) LIMITS')
         CALL MESAGE ('    SP*AWN      = SPAWNS A SUBPROCESS')
         CALL MESAGE ('    EX*IT       = EXITS FASTQ')
         CALL MESAGE ('                  (CARRIAGE RETURN TO EXIT '//
     &      'TABLET)')
      END IF

      RETURN

      END
