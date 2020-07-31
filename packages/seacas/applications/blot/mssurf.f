C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MSSURF (A, IA, LENE, NLNKE, LINKE,
     &   LENF, NLNKF, KLINKF, KIF2EL, NAMELB)
C=======================================================================

C   --*** MSSURF *** (MESH) Convert elements into faces
C   --   Written by Amy Gilkey - revised 03/29/88
C   --
C   --MSSURF converts 3D elements into faces and packs the faces by
C   --element block.
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   LENE - IN - the cumulative element counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   LINKE - IN - the original connectivity; connectivity all zero
C   --      if element is undefined
C   --   LENF - OUT - the cumulative face counts by element block
C   --   NLNKF - OUT - the number of nodes per face
C   --   KLINKF - OUT - the dynamic memory index of the connectivity
C   --      for all faces
C   --   KIF2EL - OUT - the dynamic memory index of the element number
C   --      of each face
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      include 'params.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      DIMENSION A(*)
      INTEGER IA(*)
      INTEGER LENE(0:*)
      INTEGER NLNKE(*)
      INTEGER LINKE(*)
      INTEGER LENF(0:*)
      INTEGER NLNKF(*)
      CHARACTER*(MXSTLN) NAMELB(*)

      IF (.NOT. IS3DIM) THEN

         CALL CNTLNK (NELBLK, LENE, NLNKE, LENLNK, MAXFAC)
         CALL MDLONG ('LINKF', KLINKF, LENLNK)
         CALL MDLONG ('IF2EL', KIF2EL, MAXFAC)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL SURF2D (LENE, NLNKE, LINKE,
     &      LENF, NLNKF, A(KLINKF), A(KIF2EL),
     &      MAXFAC, LENLNK)

      ELSE

         CALL MDRSRV ('NLKSCR', KNLNKS, NELBLK)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL CNTLK3 (NELBLK, LENE, NLNKE, LSCLNK, NSCFAC, NAMELB,
     $        A(KNLNKS))

         CALL MDRSRV ('LENSCR', KLENS, 1 + NELBLK)
         CALL MDRSRV ('LNKSCR', KLINKS, LSCLNK)
         CALL MDRSRV ('IFSCR', KF2ES, 2 * NSCFAC)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL FACE3D (A, LENE, NLNKE, LINKE, MAXFAC, LENLNK,
     $        A(KLENS), A(KNLNKS), A(KLINKS), A(KF2ES), NAMELB)

         CALL CNTLNK (NELBLK, IA(KLENS), IA(KNLNKS), LSCLNK, NSCFAC)
         CALL MDLONG ('LNKSCR', KLINKS, LSCLNK)
         CALL MDLONG ('IFSCR', KF2ES, 2 * NSCFAC)

         CALL MDLONG ('LINKF', KLINKF, LENLNK)
         CALL MDLONG ('IF2EL', KIF2EL, MAXFAC)
C      --MSSTEP and MSGEOM uses MDFIND to find IF2EL2
         CALL MDRSRV ('IF2EL2', KIF2E2, MAXFAC)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL SURF3D (A(KLENS), A(KNLNKS), A(KLINKS), A(KF2ES),
     &      LENF, NLNKF, A(KLINKF), A(KIF2EL), A(KIF2E2),
     &      MAXFAC, LENLNK)

         CALL MDDEL ('LENSCR')
         CALL MDDEL ('NLKSCR')
         CALL MDDEL ('LNKSCR')
         CALL MDDEL ('IFSCR')
      END IF

      CALL MDLONG ('LINKF', KLINKF, LENLNK)
      CALL MDLONG ('IF2EL', KIF2EL, MAXFAC)
      IF (IS3DIM) CALL MDLONG ('IF2EL2', KIF2E2, MAXFAC)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

  100 CONTINUE
      RETURN
      END
