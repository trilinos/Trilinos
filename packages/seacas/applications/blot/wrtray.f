C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE wrtray (LENF, NLNKF, LINKF, NXFAC, IXFAC,
     &   XN, YN, ZN, IELBST, BLKCOL, IDELB, *)
C=======================================================================

C   --*** WRTRAY *** (DETOUR) Write polygons to rayshade input file
C   --
C   --WRTRAY writes the polygons to a rayshade input file format.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   NXFAC - IN - the number of ordered faces (if DOIXF)
C   --   IXFAC - IN - the indices of the ordered faces (if DOIXF)
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

      include 'dbname.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IXFAC(*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER IELBST(NELBLK)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)

      REAL XPTS(0:20), YPTS(0:20), ZPTS(0:20)
      CHARACTER*20 STRING
      CHARACTER*2048 FILNAM
      SAVE NFIL
      DATA NFIL   /0/

C ... Open the output file.
      filnam = basenam(:lenstr(basenam)) // '.ray'
      lnam = lenstr(filnam)
      NFIL = NFIL + 1
      IF (NFIL .GT. 1) then
        call intstr (1, -1, NFIL, STRING, LSTR)
        FILNAM(LNAM+1:LNAM+LSTR) = STRING(:LSTR)
        LNAM = LNAM + LSTR
      END IF
      open (unit=95, file=filnam(:lnam), form='formatted',
     *  status='unknown')

      PI = ATAN2(0.0, -1.0)

      xmin = 1.0e38
      ymin = 1.0e38
      zmin = 1.0e38
      xmax =-1.0e38
      ymax =-1.0e38
      zmax =-1.0e38

      LSTBLK = -999
      DO 100 IX = 1, NXFAC
         IFAC = IXFAC(IX)
         IELB = 0
         IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)

         if (lstblk .ne. ielb) then
           write (string, 900) idelb(ielb)
           call pckstr (1, string)
           lstr = lenstr(string)
           lstblk = ielb
         end if

         NNPF = NLNKF(IELB)
         do 90 ilink = 0, nnpf-1
           XPTS(ILINK) = XN(LINKF(IXL+ILINK))
           YPTS(ILINK) = YN(LINKF(IXL+ILINK))
           ZPTS(ILINK) = ZN(LINKF(IXL+ILINK))

           xmin = min(xmin, xpts(ilink))
           ymin = min(ymin, ypts(ilink))
           zmin = min(zmin, zpts(ilink))

           xmax = max(xmax, xpts(ilink))
           ymax = max(ymax, ypts(ilink))
           zmax = max(zmax, zpts(ilink))
 90      continue

      write (95, 910) string(:lstr),
     *  (xpts(i), ypts(i), zpts(i),i=0, nnpf-1)

  100 CONTINUE

      write (*,*) 'RAYSHADE polygon file written to ',filnam(:lnam)
      write (*,*)
      write (*,920) 'X', xmin, xmax, xmax-xmin, (xmin+xmax)/2
      write (*,920) 'Y', ymin, ymax, ymax-ymin, (ymin+ymax)/2
      write (*,920) 'Z', zmin, zmax, zmax-zmin, (zmin+zmax)/2
      write (*,*)

      dmax  = max(xmax-xmin, ymax-ymin) / 2.0
      do 110 i=1, 9
        fov = i * 5.0
        angle = FOV/2.0 * PI/180.
        zdist = dmax / tan(angle)
        eyep  = zmax + zdist
        write (*,930) eyep, FOV
 110  continue
        write (*,*)
        close (95)

      RETURN

 900  FORMAT ('id',i10)
 910  format ('poly ',A, 20(1pe13.5))
 920  format (A1,': Min =',1pe11.3,' Max =',1pe11.3,' Range =',
     *  1pe11.3,' Center =',1pe11.3)
 930  FORMAT ('Z-coordinate of eye point =',1pe11.3,' for a ',
     *  0pf4.1,' degree field of view')
      END
