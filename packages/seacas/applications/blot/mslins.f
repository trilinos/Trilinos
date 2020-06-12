C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: mslins.f,v $
C Revision 1.3  2009/03/25 12:36:46  gdsjaar
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
C Revision 1.1  1994/04/07 20:05:43  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:53:54  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE MSLINS (A, LENF, NLNKF, LINKF, IF2EL, LENL, KLNSET)
C=======================================================================

C   --*** MSLINS *** (MESH) Create line set from surface faces
C   --   Written by Amy Gilkey - revised 03/29/88
C   --
C   --MSLINS creates a line set of all lines making up the surface faces.
C   --The lines are classified and sorted into the following parts:
C   --   -1) Edge
C   --    0) Element block boundary
C   --    n) Interior within element block 'n'
C   --
C   --Parameters:
C   --   A - IN - the dynamic memory base array
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all surface faces
C   --   IF2EL - IN - the element number of each face
C   --   LENL - OUT - the cumulative line counts by element block
C   --   KLNSET - OUT - the dynamic memory index of sorted line set,
C   --      for surface faces only
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      DIMENSION A(*)
      INTEGER LENF(0:*)
      INTEGER NLNKF(*)
      INTEGER LINKF(*)
      INTEGER IF2EL(*)
      INTEGER LENL(-2:NELBLK)

C   --Invalidate the linset array before reassigning to prevent the
C   --movement of worthless data
      CALL MDLONG ('LINSET', KLNSET, 0)

      IF (.NOT. IS3DIM) THEN
         LLNSET = 2
         CALL CNTLNK (NELBLK, LENF, NLNKF, MAXLIN, IDUM)
         CALL MDLONG ('LINSET', KLNSET, LLNSET * MAXLIN)
         CALL MDRSRV ('IEBSET', KEBSET, MAXLIN)
         CALL MDRSRV ('LINDEF', KLNDEF, (1+LLNSET) * MAXLIN)
         CALL MDRSRV ('NREF', KNREF, NUMNPF)
         CALL MDRSRV ('LREF', KLREF, NUMNPF)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL GEOM2D (LENF, NLNKF, LINKF, IF2EL,
     &      LENL, A(KLNSET), A(KEBSET), A(KLNDEF),
     &      A(KNREF), A(KLREF))
         MAXLIN = LENL(NELBLK)

         CALL MDDEL ('IEBSET')
         CALL MDDEL ('LINDEF')
         CALL MDDEL ('NREF')
         CALL MDDEL ('LREF')

      ELSE
         LLNSET = 3
         CALL CNTLNK (NELBLK, LENF, NLNKF, MAXLIN, IDUM)
         CALL MDLONG ('LINSET', KLNSET, LLNSET * MAXLIN)
         CALL MDRSRV ('LINDEF', KLNDEF, 6 * MAXLIN)
         CALL MDRSRV ('LREF', KLREF, NUMNPF)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 100

         CALL GEOM3D (LENF, NLNKF, LINKF, IF2EL,
     &      LENL, A(KLNSET), A(KLNDEF), A(KLREF))
         MAXLIN = LENL(NELBLK)

         CALL MDDEL ('LREF')
         CALL MDDEL ('LINDEF')
      END IF

      CALL MDLONG ('LINSET', KLNSET, LLNSET * MAXLIN)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 100

  100 CONTINUE
      RETURN
      END
