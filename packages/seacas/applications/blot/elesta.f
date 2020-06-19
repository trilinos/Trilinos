C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: elesta.f,v $
C Revision 1.3  2009/03/25 12:36:43  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:21  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:00:16  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:49:46  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE ELESTA (ISYTYP, VARFAC, LENF, NLNKF, LINKF, HIDEF,
     &   XN, YN, ZN, ISVOK, *)
C=======================================================================

C   --*** ELESTA *** (DETOUR) Plot element symbol for state
C   --   Written by Amy Gilkey, revised 03/14/88
C   --   D. P. Flanagan, 12/08/83
C   --
C   --ELESTA drives the state symbol interface for element variables.
C   --It processes each element by element block and plots the element
C   --state.  Only elements in selected element blocks are drawn.
C   --
C   --Parameters:
C   --   ISYTYP - IN - the symbol type (as in MODTYP of /DETOPT/)
C   --   VARFAC - IN - the face variable values
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   ISVOK - IN - ISVOK(i) is true iff the symbol variable is defined
C   --      for element block i
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      CHARACTER*(*) ISYTYP
      REAL VARFAC(*)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      LOGICAL HIDEF(*)
      REAL XN(*), YN(*), ZN(*)
      LOGICAL ISVOK(NELBLK)

      LOGICAL GRABRT

      DO 110 IELB = 1, NELBLK
         IF (ISVOK(IELB)) THEN

            DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IF (IS3DIM) THEN
                  IF (HIDEF(IFAC)) GOTO 100
               END IF

               IVAR = NINT (VARFAC(IFAC))
               IF (IVAR .GT. 0) THEN
                  CALL GRCOLR (IVAR)

C               --Draw a box around the element

                  IF (GRABRT ()) RETURN 1
                  IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)
                  CALL ELESTF (NLNKF(IELB), LINKF(IXL), XN, YN, ZN)
               END IF
  100       CONTINUE

            CALL PLTFLU
         END IF
  110 CONTINUE

      RETURN
      END
