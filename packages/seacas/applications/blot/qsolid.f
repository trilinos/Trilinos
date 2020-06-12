C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: qsolid.f,v $
C Revision 1.3  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:22  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:08:52  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.3  1991/06/25  16:09:52  gdsjaar
c Fixed? problem with calls to ugrcol -- changed
c call ugrcol(idelb(ielb),...) to call ugrcol(ielb,...)
c
c Revision 1.2  1990/12/14  08:55:48  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE QSOLID (LENF, NLNKF, LINKF, HIDEF, XN, YN, ZN, IELBST,
     &   BLKCOL, IDELB, *)
C=======================================================================

C   --*** QSOLID *** (DETOUR) Paint solid mesh (quick)
C   --   Modified by John Glick - 11/29/88
C   --   Written by Amy Gilkey - revised 10/27/87
C   --
C   --QSOLID paints the mesh in the color of each element's element block.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   HIDEF(i) - IN - true iff face i is hidden (3D only)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NUMEL, NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      LOGICAL HIDEF(*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER IELBST(NELBLK)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)

      LOGICAL GRABRT

      DO 110 IELB = 1, NELBLK

         CALL UGRCOL (IELB, BLKCOL)
C         CALL UGRCOL (IDELB(IELB), BLKCOL)
         IF ( (.NOT. IS3DIM)  .AND.  (NLNKF(IELB) .EQ. 9)) THEN
            NNPF = 8
         ELSE
            NNPF = NLNKF(IELB)
         ENDIF

         DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IF (IS3DIM) THEN
               IF (HIDEF(IFAC)) GOTO 100
            END IF

            IF (GRABRT ()) RETURN 1
            IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)
            CALL SOLIDF (NNPF, LINKF(IXL), XN, YN, ZN)
  100    CONTINUE

         CALL PLTFLU
  110 CONTINUE

      RETURN
      END
