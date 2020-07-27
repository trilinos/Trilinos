C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE FACE3D (A, LENE, NLNKE, LINKE,
     &   NFACES, LENLNK, LENSC, NLNKSC, LINKSC, IF2ESC, NAMELB)
C=======================================================================

C   --*** FACE3D *** (MESH) Break 3D elements into faces
C   --   Written by Amy Gilkey - revised 03/29/88
C   --              Sam Key, 04/85
C   --
C   --FACE3D finds all the faces in the 3D mesh.
C   --
C   --Dynamic memory is manipulated by the routine and should be checked
C   --for errors after call.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LINKE - IN - the original connectivity; connectivity all zero
C   --      if element is undefined
C   --   NFACES - OUT - returned number of non-duplicate faces
C   --   LENLNK - OUT - the length of the face connectivity array without
C   --      duplicates
C   --   LENSC - OUT - the cumulative face counts by element block
C   --   NLNKSC - OUT - the number of nodes per face
C   --   LINKSC - OUT - the face connectivity
C   --   IF2ESC - OUT - the element number(s) of each face in LINKSC;
C   --      size = 2*NFACES
C   --      IF2ESC(2,x) = 0 iff surface face
C   --      IF2ESC(2,x) > 0 iff interior face
C   --      IF2ESC(2,x) < 0 iff duplicate interior face
C   --
C   --Common Variables:
C   --   Uses NUMNP, NELBLK of /DBNUMS/

      common /debugc/ cdebug
      common /debugn/ idebug
      character*8 cdebug

      include 'params.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'minmax.blk'

      DIMENSION A(*)
      INTEGER LENE(0:NELBLK), LINKE(*)
      INTEGER NLNKE(NELBLK)
      INTEGER LENSC(0:NELBLK), LINKSC(*)
      INTEGER NLNKSC(NELBLK)
      INTEGER IF2ESC(2,*)
      CHARACTER*(MXSTLN) NAMELB(*)

      logical hastet

      INTEGER NPFRAT(8)
      SAVE NPFRAT

C      --NPFRAT(i) is the maximum number of faces that a node can be
C        in if the number of elements that the node is in = i
C        NOTE: Will fail for collapsed hexes...
C        For i <= 4, assumes possible strange connections, for
C        i > 4, assumes a somewhat regular connectivity
      DATA NPFRAT / 3, 6, 9, 12, 12, 12, 12, 12 /
C   --Reserve storage for node-to-face pointers

C   --MAXNPF is the maximum number of faces with a given node
C   --that are unmatched as the face matching routine is executed
      CALL MDRSRV ('NFPN', KNFPN, 1+NUMNP)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 150
      CALL INIINT (1+NUMNP, 0, A(KNFPN))
      MAXEPN = 0
      HASTET = .FALSE.
      DO 100 IELB = 1, NELBLK
         IXL = IDBLNK (IELB, 0, LENE, NLNKE)
         NUME = LENE(IELB) - LENE(IELB-1)
         CALL MXEPN (NUME, NLNKE(IELB), LINKE(IXL), A(KNFPN))
         if (namelb(ielb)(:3) .eq. 'TET') hastet = .true.
  100 CONTINUE
      MAXEPN = MXEPNMX(A(KNFPN))

C ... NOTE: I think this needs to be fixed for tet elements.
      IF (MAXEPN .LE. 8) THEN
         MAXNPF = NPFRAT(MAXEPN)
      ELSE
         if (hastet) then
C ... Not sure whether this is the correct number or not...
            MAXNPF = MAXEPN * 3
         else
            MAXNPF = (MAXEPN * 3) / 2
         end if
      END IF
      CALL MDDEL ('NFPN')

      CALL MDRSRV ('NPFS', KNPFS, (1+MAXNPF) * NUMNP)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 150

C   --Find all faces, matching faces in same element block

      LENSC(0) = 0
      IXL = 1
      minnod = 1
      maxnod = numnp
      DO 110 IELB = 1, NELBLK
         IXE = IDBLNK (IELB, 0, LENE, NLNKE)
         NUME = LENE(IELB) - LENE(IELB-1)
         CALL FACELB (IELB, LENE, NLNKE(IELB), LINKE(IXE),
     &      NLNKSC(IELB), LINKSC(IXL), IF2ESC(1,LENSC(IELB-1)+1),
     &      MAXNPF, A(KNPFS), NFELB, NAMELB(IELB))
         LENSC(IELB) = LENSC(IELB-1) + NFELB
         IXL = IXL + NLNKSC(IELB) * NFELB
  110 CONTINUE

      NFACES = LENSC(NELBLK)
      LENLNK = IXL - 1

C   --Match faces in different element blocks

      NOVER = 0
      minnod = 1
      maxnod = numnp

      DO 140 IELBLK = 1, NELBLK
        if (nlnksc(ielblk) .ge. 4) then
          CALL ININPF (minnod, maxnod, MAXNPF, A(KNPFS))
          minnod = numnp
          maxnod = 1
          DO 130 IELB = IELBLK, NELBLK
            IF ((NLNKSC(IELB) .EQ. NLNKSC(IELBLK))
     &        .AND. (NAMELB(IELB)(:3) .NE. 'SHE')
     $        .AND. (NAMELB(IELB) .EQ. NAMELB(IELBLK))) THEN
              IXL = IDBLNK (IELB, 0, LENSC, NLNKSC)
              DO 120 IFAC = LENSC(IELB-1)+1, LENSC(IELB)
                IF (IF2ESC(2,IFAC) .EQ. 0) THEN
                  IF (IELB .NE. IELBLK) THEN
                    if (NAMELB(IELB)(:3) .eq. 'TET') THEN
                      IMATCH = MATFAT
     $                  (LINKSC(IXL), MAXNPF, A(KNPFS),
     $                  IF2ESC(1,IFAC), numnp, IERR)
                    ELSE
                      IMATCH = MATFAC
     &                  (LINKSC(IXL), MAXNPF, A(KNPFS),
     *                  IF2ESC(1,IFAC), numnp, IERR)
                    END IF
                  ELSE
                    IMATCH = 0
                  END IF

                  IF (IMATCH .GT. 0) THEN
                    NFACES = NFACES - 1
                    LENLNK = LENLNK - NLNKSC(IELB)
                    IF2ESC(2,IMATCH) = IF2ESC(1,IFAC)
                    IF2ESC(1,IFAC) = -999
                    IF2ESC(2,IFAC) = -999
                  ELSE
                    CALL FILNPF (NLNKSC(IELB), LINKSC(IXL), IFAC,
     &                MAXNPF, A(KNPFS), NOVER, numnp)
                  END IF
                END IF
                IXL = IXL + NLNKSC(IELB)
 120          CONTINUE
            END IF
 130      CONTINUE
        end if
 140  CONTINUE

      CALL MDDEL ('NPFS')
  150 CONTINUE
      RETURN
      END
