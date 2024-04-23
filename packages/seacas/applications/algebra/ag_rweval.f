C    Copyright(C) 1999-2020, 2023, 2024 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RWEVAL (NDBIN, NDBOUT, A, ia, C, NPTIMS, NUMSTO,
     &  LTMENT, MAXSTK, NWRIT, IOERR, MERR)
C=======================================================================

C   --*** RWEVAL *** (ALGEBRA) Read, evaluate, write database
C   --   Written by Amy Gilkey - revised 07/13/89
C   --   Modified for EXODUSIIV2 format 8/28/95
C   --
C   --RWEVAL writes the output database.  It first writes the header
C   --information, then it processes and writes the time steps.
C   --
C   --Parameters:
C   --   A         - IN  - the dynamic memory base array int's, reals, logicals
C   --   C         - IN  - the dynamic memory base array for characters
C   --   NPTIMS    - IN  - the number of selected times
C   --   NUMSTO    - IN  - total number of input and output variable storage
C   --                     locations needed, including 1 for time/global
C   --                     variables
C   --   LTMENT    - IN  - the length of the time/global variable
C   --                     storage entry
C   --   MAXSTK    - IN  - the maximum stack size for any equation
C   --   NWRIT     - OUT - the number of times steps written
C   --   IOERR     - OUT - Input/Output error flag
C   --   MERR      - OUT - memory error flag
C   --
C   --Common Variables:
C   --   Uses QAINFO of /PROGQA/
C   --   Uses NUMEQN, NUMENT, NAMENT, TYPENT, INXENT, VALENT, VSZENT of /ENT../
C   --   Sets NUMNPO, NUMELO, NELBO of /DBOUT/
C   --   Uses ICOBEG, ICOEND of /DBXVAR/
C   --   Uses NDBIN of /DBASE/
C   --   Uses TITLEO of /DBTITL/
C   --   Uses NDIM, NUMNP, NUMEL, NSTEPS of /DBNUMS/
C   --   Uses ISZOOM of /ZOOM/

      include 'exodusII.inc'
      include 'ag_namlen.blk'
      include 'ag_numeqn.blk'
      include 'ag_progqa.blk'
      include 'ag_ent.blk'
      include 'ag_var.blk'
      include 'ag_dbtitl.blk'
      include 'ag_dbout.blk'
      include 'ag_dbxvar.blk'
      include 'ag_dbnums.blk'
      include 'ag_dbnumg.blk'
      include 'ag_dbnumq.blk'
      include 'ag_zoom.blk'
      include 'ag_filter.blk'
      include 'ag_remove.blk'
      include 'ag_dbws.blk'

      DIMENSION A(1)
      integer   ia(1)
      CHARACTER*1 C(1)
      INTEGER NPTIMS
      INTEGER NUMSTO
      INTEGER LTMENT

      LOGICAL STEP1, WSTEP1
      INTEGER MERR
      INTEGER CERR, NERR
      LOGICAL MEMBUG
      DATA MEMBUG /.FALSE./

      IF (MEMBUG) THEN
        CALL MLIST()
      END IF

      MERR  = 0
      IOERR = 0
      NWRIT = 0

C************************************************************
C     Reserve Dynamic arrays needed to process information
C************************************************************
C     Cumulative element counts for each element block
      CALL MDRSRV ('IXELB' , KIXELB, 1+NELBLK)
C     Cumulative element counts for each output block
      CALL MDRSRV ('IXELBO', KIXEBO, 1+NELBLK)

C************************************************************
C     Find dynamic array indices for element block information
C************************************************************
C     Find element block indices of arrays already reserved
      CALL MDFIND ('IDELB' , KIDELB, NELBLK)
      CALL MDFIND ('NUMELB', KNELB, NELBLK)
      CALL MDFIND ('NUMLNK', KNLNK, NELBLK)
      CALL MDFIND ('NUMATR', KNATR, NELBLK)
      CALL MCFIND ('BLKTYP', KNMLB, NELBLK * MXSTLN)
      CALL MDFIND ('LINK',   KLINK, IELNK)
      CALL MDFIND ('ATRIB',  KATRIB, IEATR)

      IF (NSTEPS .GT. 0) THEN
        CALL MDFIND ('IPTIMS', KPTIMS, NSTEPS)
        CALL MDFIND ('TIMES',  KTIMES, NSTEPS)
      END IF
      CALL MDFIND ('VISELB', KVISEB, NELBLK)
      IF (NVAREL .GT. 0) THEN
        CALL MDFIND ('ISEVOK', KIEVOK, NELBLK * NVAREL)
      END IF
      CALL MCFIND ('NAMECO', KNACOR, namlen*NDIM)

C     Check for errors
      CALL MDSTAT (NERR, MEM)
      CALL MCSTAT (CERR, MEM)
      IF ((NERR .GT. 0) .AND. (CERR .GT. 0)) THEN
        CALL MEMERR
        MERR = 1
        RETURN
      END IF

C ****************************************************************
C     Calculate the zoom mesh, if any
C     numeql - count number of occurrences of logical in a list
C ****************************************************************

      IF ((NUMEQL (.TRUE., NELBLK, A(KVISEB)) .LT. NELBLK)
     &  .OR. ISZOOM .OR. ISFILTER .OR. ISREMOVE) THEN
C        Read the coordinates (if needed for zoom)
        IF (ISZOOM) THEN
          CALL MDRSRV ('XN', KXN, NUMNP)
          CALL MDRSRV ('YN', KYN, NUMNP)
          IF (NDIM .GE. 3) THEN
            CALL MDRSRV ('ZN', KZN, NUMNP)
          ELSE
            KZN = 1
          END IF
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) THEN
            CALL MEMERR
            MERR = 1
            RETURN
          END IF
C           Read the coordinates of the nodes
          CALL EXGCOR(NDBIN, A(KXN), A(KYN), A(KZN), IERR)
        ELSE
          KXN = 1
          KYN = 1
          KZN = 1
        END IF

C        Indices of the zoomed node
        CALL MDRSRV ('IXNODE', KXNODE, NUMNP)
C        The zoom mexh index for each node
        CALL MDRSRV ('NODIX' , KNODIX, NUMNP)
C        Indices of the zoomed elements
        CALL MDRSRV ('IXELEM', KXELEM, NUMEL)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) THEN
          CALL MEMERR
          MERR = 1
          RETURN
        END IF

        CALL NUM2IX (NELBLK, A(KNELB), A(KIXELB))
C        A(KIXELB) - stores a cumulative count of the number of elements
C                    in each element block
C        A(KIXELB+0)= 0
C        A(KIXELB+1)= num_elem_blk1
C        A(KIXELB+2)= num_elem_blk1 + num_elem_blk2
C        A(KIXELB+3)= num_elem_blk1 + num_elem_blk2 + num_elem_blk3 ...
C        ZMFIXD finds the nodes within the zoomed coordinates and the elements
C        within the zoomed mesh (if any nodes are in).  The connectivity
C        array is fixed to point to the renumbered nodes.  The number of
C        element blocks is recalculated.
        if (isfilter) then
C ... Get space to hold the element variable used for filtering...
          call mdrsrv('ELVAR', KELVAR, numel)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) THEN
            CALL MEMERR
            MERR = 1
            RETURN
          END IF

C ... Find the step
          istep = locrea(timflt, nsteps, a(ktimes))

          IEL = 0
          DO IELB = 1, NELBLK
            IDEB = iarray(ia(kidelb), ielb)
            NELEM = iarray(ia(knelb), ielb)
            IF (i2array(ia(kievok), nelblk, ielb, idxflt) .ne. 0) THEN
              call exgev(ndbin, istep, idxflt, IDEB, NELEM,
     *          a(kelvar+iel), ierr)
            ELSE
C                  --Make sure values for undefined elements are zero
              DO I = IEL, IEL+iarray(ia(knelb), ielb) - 1
                a(kelvar+i) = 0.0
              end do
            END IF
            IEL = IEL + NELEM
          end do
          CALL FILTEL (A(KIXELB), A(KNLNK), A(KLINK), A(KVISEB),
     *      A(KXNODE), A(KIXEBO), A(KXELEM), A(KNODIX),
     *      ia(KIEVOK+(idxflt-1)*nelblk), a(kelvar))
          call MDDEL('ELVAR')
        else if (isremove) then
          CALL MDRSRV ('mapel', kmapel, NUMEL)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) THEN
            CALL MEMERR
            MERR = 1
            RETURN
          END IF
          call exgenm (ndbin, a(kmapel), ierr)
          CALL REMEL (A(KIXELB), A(KNLNK), A(KLINK), A(KVISEB),
     *      A(KXNODE), A(KIXEBO), A(KXELEM), A(KNODIX), A(KMAPEL))
          call mddel('mapel')
        else
          CALL ZMFIXD (A(KXN), A(KYN), A(KZN),A(KIXELB),
     &      A(KNLNK), A(KLINK), A(KVISEB), A(KXNODE),
     &      A(KIXEBO), A(KXELEM), A(KNODIX))
        end if
        IF ((NUMNPO .EQ. 0) .AND. (NUMELO .EQ. 0) .AND.
     &    (NELBO .EQ. 0)) THEN
          WRITE(*,500)
          WRITE(*,500)'ERROR - Output database will not be complete'
 500      FORMAT(A)
          IOERR = 1
          RETURN
        END IF

C        delete coordinate arrays - no longer needed
        IF (ISZOOM) THEN
          CALL MDDEL ('XN')
          CALL MDDEL ('YN')
          IF (NDIM .GE. 3) CALL MDDEL ('ZN')
        END IF

C        Resize the dynamic memory arrays
C        IXNODE - number of nodes to output - indices of zoomed nodes
C        IXELEM - number of elements to output - indices of zoomed elements
C        NUMNPO and NUMELO are computed in ZMFIXD
        CALL MDLONG ('IXNODE', KXNODE, NUMNPO)
        CALL MDLONG ('IXELEM', KXELEM, NUMELO)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) THEN
          CALL MEMERR
          MERR = 1
          RETURN
        END IF

C ****************************************************************
C        Read and munch the node sets
C ****************************************************************

C        output number of node sets
        NNPSO  = NUMNPS
C        output length of the concatenated node sets node list
        LNPSNO = LNPSNL

C        IF (num_nodes_in <> num_nodes_out) THEN
        IF ((NUMNP .NE. NUMNPO).AND. (numnps .gt. 0)) THEN
C           KIDNS - array containing the node set ID's for each node set
          CALL MDRSRV ('IDNPS',   KIDNS, NUMNPS)
C           KNNNS - array containing the number of nodes for each node set
          CALL MDRSRV ('NNNPS',   KNNNS, NUMNPS)
C           KNDNPS - array containing number of dist. fact for each node set
          CALL MDRSRV ('NDNPS',  KNDNPS, NUMNPS)
C           KIXNNS - array containing indices into the LTNNPS array which
C                    are the location of the 1st nodes for each set
          CALL MDRSRV ('IXNNPS', KIXNNS, NUMNPS)
C           KDISNS - array containing indices into the FACNPS array which
C                    are the location of the 1st dist factor for each set
          CALL MDRSRV ('IXDNPS', KDISNS, NUMNPS)
C           KLTNNS -  Returned array containing the nodes for all node sets
C                     Internal node IDs
          CALL MDRSRV ('LTNNPS', KLTNNS, LNPSNL)
C           KFACNS - Returned array containing the distribution factors
C                    for all node sets
          CALL MDRSRV ('FACNPS', KFACNS, LNPSNL)
          CALL MDSTAT(NERR,MEM)
          IF (NERR .GT. 0) THEN
            CALL MEMERR
            MERR = 1
            RETURN
          END IF

C           if (num_of_node_sets > 0) then
          IF (numnps .gt. 0) then
             if (LNPSNL .gt. 0) THEN
C     Read the concatenated node sets
                CALL EXGCNS(ndbin, ia(kidns), ia(knnns), ia(kndnps),
     &               ia(kixnns), ia(kdisns), ia(kltnns), a(kfacns),
     $               ierr)

C     Reserve scratch arrays for subroutine ZMNPS
                CALL MDRSRV ('NEWIX', KNEWIX, NUMNP)
                CALL MDRSRV ('IXNPS', KIXNPS, LNPSNO)
                CALL MDSTAT (NERR, MEM)
                IF (NERR .GT. 0) THEN
                   CALL MEMERR
                   MERR = 1
                   RETURN
                END IF
C     Compress the node set information by renumbering the
C     nodes and removing deleted nodes
                CALL ZMNPS (NUMNP, NUMNPO, A(KXNODE), NNPSO, LNPSNO,
     &               IA(KIDNS), IA(KNNNS), IA(KIXNNS), IA(KNDNPS),
     $               IA(KDISNS), IA(KLTNNS), A(KFACNS), IA(KNEWIX),
     $               IA(KIXNPS))
C     Remove scratch arrays
                CALL MDDEL ('NEWIX')
                CALL MDDEL ('IXNPS')
                CALL MDSTAT (NERR, MEM)
                IF (NERR .GT. 0) THEN
                   CALL MEMERR
                   MERR = 1
                   RETURN
                END IF
             ELSE
C ... Even though there are no nodes in the node set, we have told the output
C     database that there are (empty) node sets.  We need to read and write
C     the node set ids
                CALL EXGNSI (NDBIN, IA(KIDNS), IERR)
             END IF
          END IF
       END IF

C ****************************************************************
C        Read and munch the side sets
C ****************************************************************

C        Output number of side sets
        NESSO  = NUMESS
C        Output length of the side set element list
        LESSEO = LESSEL
C        Output length of the side set distribution list
        LESSDO = LESSDF

C        if (num_of_elements <> number of elements to output) then
        IF ((NUMEL .NE. NUMELO) .AND. (numess .gt. 0)) THEN
C           Reserve Dynamic arrays needed to process information
C           IDESS - array containing side set IDS
          CALL MDRSRV ('IDESS' , KIDSS , NUMESS)
C           NEESS - array containing the number of sides for each sets
          CALL MDRSRV ('NEESS' , KNESS , NUMESS)
C           KNDSS - Returned array containing the number of dist
C                   factors for each set
          CALL MDRSRV ('NDESS' , KNDSS , NUMESS )
C           IXEESS - returned array containing the indices into the
C                    LTEESS array which are the locations of the 1st
C                    element of each set
          CALL MDRSRV ('IXEESS', KIXESS, NUMESS)
C           IXDESS - Returned array containing the indices into the
C                    FACESS array which are the locations of the 1st
C                    distribution factor for each set.
          CALL MDRSRV ('IXDESS', KIDESS, NUMESS)
C           LTEESS - Returned array containing the elements for all side
C                    sets. Internal element IDS are used in this list
          CALL MDRSRV ('LTEESS', KLTESS, LESSEL)
C           LTSESS - Returned array containing the sides for all side sets
          CALL MDRSRV ('LTSESS', KLTSSS, LESSEL)
C           FACESS - Returned array containing dist factors for all side sets
          CALL MDRSRV ('FACESS', KFACSS, LESSDF)
C           LTNNN - array of number of nodes for each side in a side set
          CALL MDRSRV ('LTNNN', KLTNNN, LESSEL)
          CALL MDSTAT(NERR, MEM)
          IF (NERR .GT. 0) THEN
            CALL MEMERR
            MERR = 1
            RETURN
          END IF

C           if (number_of_side_sets > 0) then
          IF (numess .gt. 0 .and. lessel .gt. 0) THEN
C              Read concatenated side sets
            CAll EXGCSS(NDBIN, A(KIDSS), A(KNESS), A(KNDSS),
     &        A(KIXESS), A(KIDESS), A(KLTESS), A(KLTSSS),
     &        A(KFACSS), IERR)

            CALL EXGCSSC(NDBIN, A(KLTNNN), IERR)

C              Reserve scratch arrays for ZMESS
            CALL MDRSRV ('NEWIX', KNEWIX, MAX (NUMEL, NUMNP))
            CALL MDRSRV ('NEWSD', KNEWSD, LESSEL)
            CALL MDRSRV ('IXESS', KXESS, LESSEL)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) THEN
              CALL MEMERR
              MERR = 1
              RETURN
            END IF

            CALL ZMESS (NUMEL, NUMELO, A(KXELEM),
     &        NUMNP, NUMNPO, A(KXNODE), NESSO, LESSEO, LESSDO,
     &        A(KIDSS), A(KNESS), A(KIXESS), A(KLTESS), A(KLTSSS),
     &        A(KNDSS), A(KLTNNN), A(KFACSS), A(KNEWIX), A(KNEWSD),
     &        A(KXESS))

C              Delete scratch arrays
            CALL MDDEL ('LTNNN')
            CALL MDDEL ('IXESS')
            CALL MDDEL ('NEWSD')
            CALL MDDEL ('NEWIX')
            CALL MDSTAT(NERR, MEM)
            IF (NERR .GT. 0) THEN
              CALL MEMERR
              MERR = 1
              RETURN
            END IF
          END IF
        END IF

      ELSE

C        Initialize the number of output variables
        numnpo = numnp
        numelo = numel
        nelbo  = nelblk
        nnpso  = numnps
        lnpsno = lnpsnl
        nesso  = numess
        lesseo = lessel
        lessdo = lessdf
        KXNODE = 1
        KXELEM = 1
        KNODIX = 1
      END IF

C ***************************************************************
C     WRITE TO THE OUTPUT FILE
C ***************************************************************

C ****************************************************************
C     Write initialization parameters to the EXODUSII database
C ****************************************************************

      call expini (ndbout, titleo, ndim, numnpo, numelo, nelbo,
     &  nnpso, nesso, ierr)

C ****************************************************************
C     Write the QA and Information records
C ****************************************************************
C     Find dynamic array indices for QA and Information record information
      call MCFIND ('QAREC',  kqarec, 4*(MXSTLN)*(nqarec+1))
      call MCFIND ('INFREC', kinfo,  (MXLNLN)*ninfo)

      CALL DBOQA(NDBOUT, QAINFO, NQAREC, C(KQAREC),
     *  NINFO, C(KINFO))
C     Delete unused dynamic memrory
      CALL MCDEL ('QAREC')
      CALL MCDEL ('INFREC')

C ****************************************************************
C     Write element block information to the database
C     including element block connectivity and attributed
C ****************************************************************

      CALL WELB (NDBOUT, NELBLK, A(KVISEB),
     &  (NUMEL .EQ. NUMELO), C(KNMLB), A(KNLNK), A(KNATR),
     &  A(KLINK), A(KATRIB), A(KNELB), A(KIXELB), A(KIXEBO),
     &  A(KXELEM),(NUMNP .NE. NUMNPO), A(KNODIX),
     &  A(KIDELB), A, IA, C, MERR)

      IF (NUMNP .NE. NUMNPO) CALL MDDEL ('NODIX')
      CALL MDDEL ('LINK')
      CALL MDDEL ('ATRIB')
      CALL MDDEL ('NUMLNK')
      CALL MDDEL ('NUMATR')

C ****************************************************************
C     Write the database names
C ****************************************************************

C     Write the element block variable truth table
C     Reserve scratch arrays for writing truth table
      if ((nelbo .gt. 0) .and. (nelblk .ne. nelbo)) then
        CALL MDRSRV ('XEVOK', KEVOK, NELBO*NVAREO)
        CALL MDRSRV ('NEWID', KNEWID, NELBO)
      else
        CALL MDRSRV ('XEVOK', KEVOK, NELBLK*NVAREO)
        CALL MDRSRV ('NEWID', KNEWID, NELBLK)
      end if
      CALL MDSTAT(NERR, MEM)
      CALL MCSTAT(CERR, MEM)
      IF ((NERR .GT. 0) .OR. (CERR .GT. 0)) THEN
        CALL MEMERR
        MERR = 1
        RETURN
      END IF

      CALL WNAM (NDBOUT ,NDIM, NELBLK, NELBO, A(KVISEB),
     &  NVARGO, NVARNO, NVAREO, C(KNACOR), C(KNMLB),
     &  NAMVAR(JGVBEG), NAMVAR(JNVBEG), NAMVAR(JEVBEG),
     &  IEVVAR(JEVBEG), A(KIEVOK), A(KIDELB), A(KNEWID),
     &  A(KEVOK))

      CALL MDDEL ('XEVOK')
      CALL MDDEL ('NEWID')
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0)THEN
        CALL MEMERR
        MERR = 1
        RETURN
      END IF

C ****************************************************************
C     Read and write the node sets and the side sets
C ****************************************************************

C ****************************************************************
C     If number of nodes = number of nodes to output or
C     If have not read concatenated nodes sets, read now
C ****************************************************************
      IF (NUMNP .EQ. NUMNPO) THEN
C        if (number of node sets > 0) then
        IF (numnps .gt. 0) THEN
C           KIDNS - array containing the node set ID's for each node set
          CALL MDRSRV ('IDNPS', KIDNS, NUMNPS)
C           KNNNS - array containing the number of nodes for each node set
          CALL MDRSRV ('NNNPS', KNNNS, NUMNPS)
C           KNDNPS - array containing number of dist. fact for each node set
          CALL MDRSRV ('NDNPS',KNDNPS, NUMNPS)
C           KIXNNS - array containing indices into the LTNNPS array which
C                    are the location of the 1st nodes for each set
          CALL MDRSRV ('IXNNPS', KIXNNS, NUMNPS)
C           KDISNS - array containing indices into the FACNPS array which
C                    are the location of the 1st dist factor for each set
          CALL MDRSRV ('IXDNPS', KDISNS, NUMNPS)
C           KLTNNS -  Returned array containing the nodes for all node sets
C                     Internal node IDs
          CALL MDRSRV ('LTNNPS', KLTNNS, LNPSNL)
C           KFACNS - Returned array containing the distribution factors
C                    for all node sets
          CALL MDRSRV ('FACNPS', KFACNS, LNPSNL)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) THEN
            CALL MEMERR
            MERR = 1
            RETURN
          END IF

          CALL EXGCNS(ndbin, ia(kidns), ia(knnns), ia(kndnps),
     &      ia(kixnns), ia(kdisns), ia(kltnns), a(kfacns), ierr)
        END IF
      END IF
C     Write the node set information
      if ((numnps .gt. 0) .AND. (nnpso .gt. 0)) then
        CALL EXPCNS(ndbout, ia(kidns), ia(knnns), ia(kndnps),
     &    ia(kixnns), ia(kdisns), ia(kltnns), a(kfacns), ierr)

C       Delete unneeded dynamic memory
        CALL MDDEL ('IDNPS')
        CALL MDDEL ('NNNPS')
        CALL MDDEL ('NDNPS')
        CALL MDDEL ('IXNNPS')
        CALL MDDEL ('IXDNPS')
        CALL MDDEL ('LTNNPS')
        CALL MDDEL ('FACNPS')
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) THEN
          CALL MEMERR
          MERR = 1
          RETURN
        END IF
      end if

C ****************************************************************
C     If number of elements = number of elements to output or
C     If have not read the concatenated side set information yet
C ****************************************************************
      IF (NUMEL .EQ. NUMELO) THEN
        IF (numess .gt. 0) THEN
C          Reserve Dynamic arrays needed to process information
C          IDESS - array containing side set IDS
          CALL MDRSRV ('IDESS' , KIDSS , NUMESS)
C          NEESS - array containing the number of sides/elements for each sets
          CALL MDRSRV ('NEESS' , KNESS , NUMESS)
C          KNDSS - Returned array containing the number of dist
C             factors for each set
          CALL MDRSRV ('NDESS' , KNDSS , NUMESS )
C          IXEESS - returned array containing the indices into the
C              LTEESS array which are the locations of the 1st
C              element of each set
          CALL MDRSRV ('IXEESS', KIXESS, NUMESS)
C          IXDESS - Returned array containing the indices into the
C              FACESS array which are the locations of the 1st
C              distribution factor for each set.
          CALL MDRSRV ('IXDESS', KIDESS, NUMESS)
C          LTEESS - Returned array containing the elements for all side
C              sets. Internal element IDS are used in this list
          CALL MDRSRV ('LTEESS', KLTESS, LESSEL)
C          LTSESS - Returned array containing the sides for all side sets
          CALL MDRSRV ('LTSESS', KLTSSS, LESSEL)
C          FACESS - Returned array containing dist factors for all side sets
          CALL MDRSRV ('FACESS', KFACSS, LESSDF)
          CALL MDSTAT (NERR, MEM)
          IF (NERR .GT. 0) THEN
            CALL MEMERR
            MERR = 1
            RETURN
          END IF
          CAll EXGCSS(ndbin, a(kidss), a(kness), a(kndss), a(kixess),
     &      a(kidess), a(kltess), a(kltsss), a(kfacss), ierr)
        END IF
      END IF

C     Write the side set information
      IF ((numess .gt. 0) .AND. (nesso .gt. 0)) THEN
        CAll EXPCSS(ndbout, a(kidss), a(kness), a(kndss), a(kixess),
     &    a(kidess), a(kltess), a(kltsss), a(kfacss), ierr)

C        Delete dynamic memory
        CALL MDDEL ('IDESS')
        CALL MDDEL ('NEESS')
        CALL MDDEL ('NDESS')
        CALL MDDEL ('IXEESS')
        CALL MDDEL ('IXDESS')
        CALL MDDEL ('LTEESS')
        CALL MDDEL ('LTSESS')
        CALL MDDEL ('FACESS')
      END IF

C ****************************************************************
C     Read and write the element map
C ****************************************************************

      IF (NUMEL .NE. NUMELO) CALL MDRSRV ('SCREL', KSCREL, NUMEL)
      CALL MDRSRV ('MAPEL', KMAPEL, NUMEL)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
        CALL MEMERR
        MERR = 1
        RETURN
      END IF

      CALL RWMAP (NDBIN, NDBOUT, NUMEL, NUMELO, A(KXELEM), A(KMAPEL),
     &  A(KSCREL))

      CALL MDDEL ('MAPEL')
      IF (NUMEL .NE. NUMELO) CALL MDDEL ('SCREL')

C ****************************************************************
C     Read and write the node number map
C ****************************************************************

      IF (NUMNP .NE. NUMNPO) CALL MDRSRV ('SCRND', KSCRND, NUMNP)
      CALL MDRSRV ('MAPND', KMAPND, NUMNP)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
        CALL MEMERR
        MERR = 1
        RETURN
      END IF

      CALL RWNMAP (NDBIN, NDBOUT, NUMNP, NUMNPO, A(KXNODE), A(KMAPND),
     &  A(KSCRND))

      CALL MDDEL ('MAPND')
      IF (NUMNP .NE. NUMNPO) CALL MDDEL ('SCRND')

C ****************************************************************
C     Read and write the coordinates (may be used later)
C ****************************************************************

      if (numnp .ne. numnpo) then
        CALL MDRSRV('CRDSCR', KCRDSC, NDIM*NUMNPO)
      else
        kcrdsc = 1
      end if
      CALL MDRSRV ('CORD',  KCORD, NDIM*NUMNP)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
        CALL MEMERR
        MERR = 1
        RETURN
      END IF

      CALL RWXYZ (NDBIN, NDBOUT, NDIM, NUMNP, NUMNPO, A(KXNODE),
     &  A(KCORD), A(KCRDSC))

      if (numnp .ne. numnpo) then
        call mddel('CRDSCR')
      end if
      IF (ICOBEG .GT. ICOEND) THEN
        CALL MDDEL ('CORD')
      END IF

      CALL MDCOMP ()
      CALL MCCOMP ()
C     Cumulative element counts for each element block
      CALL MDFIND ('IXELB' , KIXELB, 1+NELBLK)
C     Cumulative element counts for each output block
      CALL MDFIND ('IXELBO', KIXEBO, 1+NELBLK)
      CALL MDFIND ('IDELB' , KIDELB, NELBLK)
      CALL MDFIND ('NUMELB', KNELB, NELBLK)
      CALL MCFIND ('NAMECO', KNACOR, namlen*NDIM)
      CALL MCFIND ('NAMES', KNAMES, namlen*(NVARGL+NVARNP+NVAREL))
      CALL MCFIND ('BLKTYP', KNMLB, NELBLK * MXSTLN)
      CALL MDFIND ('SELELB', KSELEB, NELBLK)

      IF (NSTEPS .GT. 0) THEN
        CALL MDFIND ('IPTIMS', KPTIMS, NSTEPS)
        CALL MDFIND ('TIMES',  KTIMES, NSTEPS)
      END IF
      CALL MDFIND ('VISELB', KVISEB, NELBLK)
      IF (NVAREL .GT. 0) THEN
        CALL MDFIND ('ISEVOK', KIEVOK, NELBLK * NVAREL)
      END IF
      if (iszoom .or. isfilter .or. (nelblk .ne. nelbo) .or.
     $     (numnpo .ne. numnp)) then
        call mdfind ('IXNODE', KXNODE, NUMNPO)
      else
        KXNODE = 1
      end if
      if (iszoom .or. isfilter .or. (nelblk .ne. nelbo) .or.
     $     (numelo .ne. numel)) then
        call mdfind ('IXELEM', KXELEM, NUMELO)
      else
        KXELEM = 1
      end if
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
        CALL MEMERR
        MERR = 1
        RETURN
      END IF

C ****************************************************************
C     Reserve space for variable stack
C ****************************************************************

      MAXNE = MAX (LTMENT, NUMNP, NUMEL)
C      --MAXNE - maximum dimension for the variable stack and the variable
C      --   arrays; an entry must be large enough to hold the time /
C      --   history/global entry, the coordinates, a nodal variable
C      --   or an element variable

      CALL MDRSRV ('VARVAL', KVVAL, MAXNE * NUMSTO)
C      --VARVAL(1..MAXNE,1..NUMSTO) - the needed input variables and
C      --   the output variables

C   --Save any coordinates needed for equation evaluation

      IF (ICOBEG .LE. ICOEND) THEN
        call mdfind ('CORD', KCORD, NDIM*NUMNP)
        CALL SVCORD (A(KCORD), MAXNE, A(KVVAL))
        CALL MDDEL ('CORD')
      END IF

C   --Reserve space for the equation evaluation stack

      CALL MDRSRV ('STACK', KSTACK, MAXNE * MAXSTK)
C      --STACK(1.MAXNE,1..MAXSTK) - an area for evaluating equations
      CALL MDRSRV ('GVSCR', KGVSCR, NVARGO)
      CALL MDRSRV ('VARSCR', KVARSC, MAXNE)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) THEN
        CALL MEMERR
        MERR = 1
        RETURN
      END IF

C   --Initialize variable read routine
      CALL DBIVIN (.TRUE.)

c      NWHOL = 0

      WRITE (*, *)

C     Loop from 1 to the number of selected times
      if (nptims .eq. 0) then
         if (numeqn .gt. 0) then
            istep = 1
            STEP1 = (NPT .EQ. 1)

C     --Evaluate the equations to get the output variables

            DO NEQN = 1, NUMEQN
               CALL EVAL (STEP1, WSTEP1, MAXNE, MAXSTK,
     &            A(KXNODE), A(KIXELB), A(KIXEBO), A(KXELEM), A(KIEVOK),
     &            NEQN, NUMENT(NEQN), NAMENT(1,NEQN), TYPENT(1,NEQN),
     &            INXENT(1,NEQN), VALENT(1,NEQN),
     &            ITMENT(1,NEQN), IEVENT(1,NEQN), VSZENT(1,NEQN),
     &            A(KVVAL), A(KSTACK), MERR)
               IF (MERR .EQ. 1) RETURN
            END DO

C     --Write the variables for the time step

            CALL WRSTEP (NDBOUT, 1, MAXNE, A(KVVAL), A(KVISEB),
     *           A(KXNODE), A(KIXELB), A(KIXEBO), A(KXELEM),
     *           A(KIDELB), A(KIEVOK), A(KGVSCR), A(KVARSC), MERR)
            IF (MERR .EQ. 1) RETURN

            WRITE (*, 10000) 1
         end if
      else
        DO 110 NPT = 1, NPTIMS
          istep = iarray(iA(KPTIMS), NPT)
          STEP1 = (NPT .EQ. 1)

c         WSTEP1 = (WHOTIM(ISTEP) .AND. (NWHOL .LE. 0))
C        Flag was present when distinguishing between
C        whole and history time steps
          WSTEP1 = .TRUE.

C      --Read variables for one time step

          CALL RDSTEP (ISTEP, IARRAY(iA(KTIMES), ISTEP), A(KNELB),
     &      A(KIDELB), A(KIEVOK), A(KVISEB),MAXNE, A(KVVAL), MERR)
          IF (MERR .EQ. 1) RETURN

C      --Fix up the first time step for last and first time variables

          IF (STEP1) THEN
            CALL FIXONE (MAXNE, A(KVVAL))
          END IF

C      --Evaluate the equations to get the output variables

          DO 100 NEQN = 1, NUMEQN
            CALL EVAL (STEP1, WSTEP1, MAXNE, MAXSTK,
     &        A(KXNODE), A(KIXELB), A(KIXEBO), A(KXELEM), A(KIEVOK),
     &        NEQN, NUMENT(NEQN), NAMENT(1,NEQN), TYPENT(1,NEQN),
     &        INXENT(1,NEQN), VALENT(1,NEQN),
     &        ITMENT(1,NEQN), IEVENT(1,NEQN), VSZENT(1,NEQN),
     &        A(KVVAL), A(KSTACK), MERR)
            IF (MERR .EQ. 1) RETURN
 100      CONTINUE

C      --Write the variables for the time step

          CALL WRSTEP (NDBOUT, NPT, MAXNE, A(KVVAL), A(KVISEB),
     *      A(KXNODE), A(KIXELB), A(KIXEBO), A(KXELEM),
     *      A(KIDELB), A(KIEVOK), A(KGVSCR), A(KVARSC), MERR)
          IF (MERR .EQ. 1) RETURN

C      --Move the values for the current time step into locations for the
C      --last time step

          CALL NXTTIM

          NWRIT = NWRIT + 1

          WRITE (*, 10000) NWRIT
10000     FORMAT (' ', I8, ' time steps processed')
c         CALL NCSNC (NDBOUT, IERR)
 110    CONTINUE
      end if

      CALL MDDEL ('VARSCR')
      CALL MDDEL ('STACK')
      CALL MDDEL ('GVSCR')
      CALL MDDEL ('VARVAL')
      CALL MDDEL ('IXELB')
      CALL MDDEL ('IXELBO')
      if (iszoom .or. isfilter .or. (nelblk .ne. nelbo)
     $        .or. isremove) then
        call mddel ('IXNODE')
        CALL MDdel ('IXELEM')
      end if
      CALL MDSTAT(NERR, MEM)
      IF (NERR .GT. 0) THEN
        CALL MEMERR
        MERR = 1
      END IF

      RETURN
      END

      integer function iarray(intarr, ipos)

      integer intarr(*)
      integer ipos

      iarray = intarr(ipos)
      return
      end

      integer function i2array(intarr, irows, irow, icol)

      integer intarr(irows,*)
      integer irow, icol

      i2array = intarr(irow, icol)
      return
      end
