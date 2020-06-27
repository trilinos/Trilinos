C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: fixal2.f,v $
C Revision 1.3  2009/03/25 12:36:44  gdsjaar
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
C Revision 1.1  1994/04/07 20:00:53  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:50:36  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE FIXAL2 (ALIVE, LENF, IF2EL, IE2ELB, NEWELB)
C=======================================================================

C   --*** FIXAL2 *** (MESH) Adjust for element birth/death (2D)
C   --   Written by Amy Gilkey - revised 03/10/88
C   --
C   --FIXAL2 adjusts the face array to reflect the new element states.
C   --Faces that are dead are moved to LENF(NELBLK+2).
C   --
C   --Parameters:
C   --   ALIVE - IN - true iff element i is alive
C   --   LENF - IN - the cumulative face counts by element block
C   --   IF2EL - IN - the element number of each face
C   --   IE2ELB - IN - the element block for each element
C   --   NEWELB - OUT - size = LENF(NELBLK+1)
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/

      include 'dbnums.blk'
      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      LOGICAL ALIVE(NUMEL)
      INTEGER LENF(0:NELBLK+2)
      INTEGER IF2EL(*)
      INTEGER IE2ELB(NUMEL)
      INTEGER NEWELB(*)

      DO 110 IELB = 1, NELBLK+2
         DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
            NEWELB(IFAC) = IELB
  100    CONTINUE
  110 CONTINUE

C   --Check each face by its defining elements, and move face if necessary

      DO 130 IELB = 1, NELBLK+2
         DO 120 IFAC = LENF(IELB-1)+1, LENF(IELB)

C         --Determine if the element for this face is alive

            IEL = IF2EL(IFAC)

            IF (ALIVE(IEL)) THEN

C            --If alive, change to an ALIVE face
               NEWELB(IFAC) = IE2ELB(IEL)

            ELSE

C            --If not alive, change to a DEAD face
               NEWELB(IFAC) = NELBLK+2

            END IF

  120    CONTINUE
  130 CONTINUE

      RETURN
      END
