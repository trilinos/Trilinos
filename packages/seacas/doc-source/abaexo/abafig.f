      program test

      logical ldum, pltgtt, pltstt, pltstv
      character*5 oldnod, newnod

      call pltint
      call pltbgn

      penlin = 160.0 !normal line width
      ldum = pltstv (2, penlin) !normal line width
      ldum = pltgtt (11, pentxt) !normal pen width
      ldum = pltgtt (2, x) !character size
      vcs = .75*x
      ldum = pltstt (2, vcs) !character size

      oldnod = '\CDO '
      newnod = '\CCI '
      xnode = .003
      ynode = -(vcs+.003)
      xelem = 0.0
      yelem = -vcs

      rng2d = .160
      rng3d = rng2d * (8./7.)
      xsft3d = rng3d/6
      ysft3d = rng3d/8

      r3d = rng3d + xsft3d + xnode+2*vcs*(5./7.)
      xbord = (1.0 - 2 * (rng2d + r3d)) / 3.6
      xcen2d = xbord + rng2d
      xcen3d = 1.0 - (xbord + r3d)
      ybord = xbord
      ycen3d = .75/2
      ycen2d = ycen3d + 2*ysft3d

      yfig = 0.0 + ybord
      ytitl = .75 - ybord - 1.25*vcs

d     call pltvct (1, 0.0, 0.0, 1.0, 0.0)
d     call pltvct (1, 1.0, 0.0, 1.0, .75)
d     call pltvct (1, 1.0, .75, 0.0, .75)
d     call pltvct (1, 0.0, .75, 0.0, 0.0)

C *** 2D 8-node element

      xl = xcen2d - rng2d
      xh = xcen2d + rng2d
      yl = ycen2d - rng2d
      yh = ycen2d + rng2d
      xm = (xl+xh) / 2
      ym = (yl+yh) / 2

C   --Plot legend

      ldum = pltgtt (6, splin) !line spacing
      ymid = (yl - vcs) / 2 + (yfig + 1.25*vcs) / 2
      y = ymid + (3*vcs*splin) / 2
      call pltxsl ('= Original element', xlen)
      xlins = xcen2d - (xlen + 6.0*vcs) / 2
      xleg = xlins + 6.0*vcs
      xline = xleg - 1.0*vcs
      xch = xleg - 1.5*vcs
      call option ('old spots', .true., vcs, penlin, pentxt)
      call pltxtc (xch, y, '\DO ')
      call option ('old spots', .false., vcs, penlin, pentxt)
      call pltxts (xleg, y, '= Original node')
      call pltxsn (x, y)
      call option ('new spots', .true., vcs, penlin, pentxt)
      call pltxtc (xch, y, '\CI ')
      call option ('new spots', .false., vcs, penlin, pentxt)
      call pltxts (xleg, y, '= Generated node')
      call pltxsn (x, y)
      call option ('old lines', .true., vcs, penlin, pentxt)
      call pltvct (1, xlins, y+vcs/2, xline, y+vcs/2)
      call option ('old lines', .false., vcs, penlin, pentxt)
      call pltxts (xleg, y, '= Original element')
      call pltxsn (x, y)
      call option ('new lines', .true., vcs, penlin, pentxt)
      call pltvct (1, xlins, y+vcs/2, xline, y+vcs/2)
      call option ('new lines', .false., vcs, penlin, pentxt)
      call pltxts (xleg, y, '= New element')

C   --Label plot

      ldum = pltstt (2, 1.25*vcs) !big letters
      call pltxtc (xcen2d, yfig, 'Figure 1a.  2D Transformation')

      call pltxtc (xcen2d, ytitl, '8-Node 2D Element =>')
      call pltxsn (x, y)
      call pltxtc (xcen2d, y, 'Four 4-Node Elements')
      ldum = pltstt (2, 1.00*vcs) !normal letters

C   --Draw original element

      call option ('old lines', .true., vcs, penlin, pentxt)
      call pltvct (1, xl, yl, xh, yl)
      call pltvct (1, xh, yl, xh, yh)
      call pltvct (1, xh, yh, xl, yh)
      call pltvct (1, xl, yh, xl, yl)
      call option ('old lines', .false., vcs, penlin, pentxt)

C   --Number original nodes

      call option ('old nodes', .true., vcs, penlin, pentxt)
      call pltxts (xl+xnode, yl+ynode, '1')
      call pltxts (xh+xnode, yl+ynode, '2')
      call pltxts (xh+xnode, yh+ynode, '3')
      call pltxts (xl+xnode, yh+ynode, '4')
      call pltxts (xm+xnode, yl+ynode, '5')
      call pltxts (xh+xnode, ym+ynode, '6')
      call pltxts (xm+xnode, yh+ynode, '7')
      call pltxts (xl+xnode, ym+ynode, '8')
      call option ('old nodes', .false., vcs, penlin, pentxt)

      call option ('old spots', .true., vcs, penlin, pentxt)
      call pltxts (xl, yl, oldnod)
      call pltxts (xh, yl, oldnod)
      call pltxts (xh, yh, oldnod)
      call pltxts (xl, yh, oldnod)
      call pltxts (xm, yl, oldnod)
      call pltxts (xh, ym, oldnod)
      call pltxts (xm, yh, oldnod)
      call pltxts (xl, ym, oldnod)
      call option ('old spots', .false., vcs, penlin, pentxt)

C   --Draw new element

      call option ('new lines', .true., vcs, penlin, pentxt)
      call pltvct (1, xl, ym, xh, ym)
      call pltvct (1, xm, yl, xm, yh)
      call option ('new lines', .false., vcs, penlin, pentxt)

C   --Number new nodes

      call option ('new nodes', .true., vcs, penlin, pentxt)
      call pltxts (xm+xnode, ym+ynode, '9')
      call option ('new nodes', .false., vcs, penlin, pentxt)

      call option ('new spots', .true., vcs, penlin, pentxt)
      call pltxts (xm, ym, newnod)
      call option ('new spots', .false., vcs, penlin, pentxt)

C   --Number new elements

      call option ('elements', .true., vcs, penlin, pentxt)
      xlm = (xl+xm) / 2
      xhm = (xh+xm) / 2
      ylm = (yl+ym) / 2
      yhm = (yh+ym) / 2
      call pltxtc (xlm+xelem, ylm+yelem, '1')
      call pltxtc (xhm+xelem, ylm+yelem, '2')
      call pltxtc (xhm+xelem, yhm+yelem, '3')
      call pltxtc (xlm+xelem, yhm+yelem, '4')
      call option ('elements', .false., vcs, penlin, pentxt)

C *** 3D 8-node element

      xl0 = xcen3d - rng3d
      xh0 = xcen3d + rng3d
      yl0 = ycen3d - rng3d
      yh0 = ycen3d + rng3d
      xl1 = xl0 - xsft3d
      xh1 = xh0 - xsft3d
      yl1 = yl0 - ysft3d
      yh1 = yh0 - ysft3d
      xl2 = xl0 + xsft3d
      xh2 = xh0 + xsft3d
      yl2 = yl0 + ysft3d
      yh2 = yh0 + ysft3d
      xm1 = (xl1+xh1) / 2
      ym1 = (yl1+yh1) / 2
      xm2 = (xl2+xh2) / 2
      ym2 = (yl2+yh2) / 2
      xm0 = (xl0+xh0) / 2
      ym0 = (yl0+yh0) / 2

C   --Label plot

      ldum = pltstt (2, 1.25*vcs) !big letters
      call pltxtc (xcen3d, yfig, 'Figure 1b.  3D Transformation')

      call pltxtc (xcen3d, ytitl, '20-Node 3D Element =>')
      call pltxsn (x, y)
      call pltxtc (xcen3d, y, 'Eight 8-Node Elements')
      ldum = pltstt (2, 1.00*vcs) !normal letters

C   --Draw original element

      call option ('old lines', .true., vcs, penlin, pentxt)
      call pltvct (1, xl1, yl1, xh1, yl1)
      call pltvct (1, xh1, yl1, xh1, yh1)
      call pltvct (1, xh1, yh1, xl1, yh1)
      call pltvct (1, xl1, yh1, xl1, yl1)
      call pltvct (1, xh2, yl2, xh2, yh2)
      call pltvct (1, xh2, yh2, xl2, yh2)
      call pltvct (1, xh1, yl1, xh2, yl2)
      call pltvct (1, xh1, yh1, xh2, yh2)
      call pltvct (1, xl1, yh1, xl2, yh2)
      call option ('old lines', .false., vcs, penlin, pentxt)

      call option ('hidden lines', .true., vcs, penlin, pentxt)
      call pltvct (1, xl2, yl2, xh2, yl2)
      call pltvct (1, xl2, yh2, xl2, yl2)
      call pltvct (1, xl1, yl1, xl2, yl2)
      call option ('hidden lines', .false., vcs, penlin, pentxt)

C   --Number original nodes

      call option ('old nodes', .true., vcs, penlin, pentxt)
      call pltxts (xl1+xnode, yl1+ynode, '1')
      call pltxts (xh1+xnode, yl1+ynode, '2')
      call pltxts (xh1+xnode, yh1+ynode, '6')
      call pltxts (xl1+xnode, yh1+ynode, '5')
      call pltxts (xl2+xnode, yl2+ynode, '4')
      call pltxts (xh2+xnode, yl2+ynode, '3')
      call pltxts (xh2+xnode, yh2+ynode, '7')
      call pltxts (xl2+xnode, yh2+ynode, '8')
      call pltxts (xm1+xnode, yl1+ynode, '9')
      call pltxts (xh1+xnode, ym1+ynode, '18')
      call pltxts (xm1+xnode, yh1+ynode, '13')
      call pltxts (xl1+xnode, ym1+ynode, '17')
      call pltxts (xm2+xnode, yl2+ynode, '11')
      call pltxts (xh2+xnode, ym2+ynode, '19')
      call pltxts (xm2+xnode, yh2+ynode, '15')
      call pltxts (xl2+xnode, ym2+ynode, '20')
      call pltxts (xl0+xnode, yl0+ynode, '12')
      call pltxts (xh0+xnode, yl0+ynode, '10')
      call pltxts (xh0+xnode, yh0+ynode, '14')
      call pltxts (xl0+xnode, yh0+ynode, '16')
      call option ('old nodes', .false., vcs, penlin, pentxt)

      call option ('old spots', .true., vcs, penlin, pentxt)
      call pltxts (xl1, yl1, oldnod)
      call pltxts (xh1, yl1, oldnod)
      call pltxts (xh1, yh1, oldnod)
      call pltxts (xl1, yh1, oldnod)
      call pltxts (xl2, yl2, oldnod)
      call pltxts (xh2, yl2, oldnod)
      call pltxts (xh2, yh2, oldnod)
      call pltxts (xl2, yh2, oldnod)
      call pltxts (xm1, yl1, oldnod)
      call pltxts (xh1, ym1, oldnod)
      call pltxts (xm1, yh1, oldnod)
      call pltxts (xl1, ym1, oldnod)
      call pltxts (xm2, yl2, oldnod)
      call pltxts (xh2, ym2, oldnod)
      call pltxts (xm2, yh2, oldnod)
      call pltxts (xl2, ym2, oldnod)
      call pltxts (xl0, yl0, oldnod)
      call pltxts (xh0, yl0, oldnod)
      call pltxts (xh0, yh0, oldnod)
      call pltxts (xl0, yh0, oldnod)
      call option ('old spots', .false., vcs, penlin, pentxt)

C   --Draw new element

      call option ('new lines', .true., vcs, penlin, pentxt)
      call pltvct (1, xh0, yl0, xh0, yh0)
      call pltvct (1, xh0, yh0, xl0, yh0)
      call pltvct (1, xl1, ym1, xh1, ym1)
      call pltvct (1, xm1, yl1, xm1, yh1)
      call pltvct (1, xh1, ym1, xh2, ym2)
      call pltvct (1, xm1, yh1, xm2, yh2)
      call option ('new lines', .false., vcs, penlin, pentxt)

      call option ('new hidden lines', .true., vcs, penlin, pentxt)
      call pltvct (1, xl0, yl0, xh0, yl0)
      call pltvct (1, xl0, yh0, xl0, yl0)
      call pltvct (1, xl2, ym2, xh2, ym2)
      call pltvct (1, xm2, yl2, xm2, yh2)
      call pltvct (1, xl0, ym0, xh0, ym0)
      call pltvct (1, xm0, yl0, xm0, yh0)
      call pltvct (1, xm1, yl1, xm2, yl2)
      call pltvct (1, xl1, ym1, xl2, ym2)
      call pltvct (1, xm1, ym1, xm2, ym2)
      call option ('new hidden lines', .false., vcs, penlin, pentxt)

C   --Number new nodes

      call option ('new nodes', .true., vcs, penlin, pentxt)
      call pltxts (xm1+xnode, ym1+ynode, '22')
      call pltxts (xm2+xnode, ym2+ynode, '24')
      call pltxts (xm0+xnode, ym0+ynode, '27')
      call pltxts (xm0+xnode, yl0+ynode, '26')
      call pltxts (xh0+xnode, ym0+ynode, '23')
      call pltxts (xm0+xnode, yh0+ynode, '25')
      call pltxts (xl0+xnode, ym0+ynode, '21')
      call option ('new nodes', .false., vcs, penlin, pentxt)

      call option ('new spots', .true., vcs, penlin, pentxt)
      call pltxts (xm1, ym1, newnod)
      call pltxts (xm2, ym2, newnod)
      call pltxts (xm0, ym0, newnod)
      call pltxts (xm0, yl0, newnod)
      call pltxts (xh0, ym0, newnod)
      call pltxts (xm0, yh0, newnod)
      call pltxts (xl0, ym0, newnod)
      call option ('new spots', .false., vcs, penlin, pentxt)

C   --Number new elements

      call option ('elements', .true., vcs, penlin, pentxt)
      x0 = (xl0-xl1) / 2
      y0 = (yl0-yl1) / 2
      xlm1 = (xl1+xm1) / 2 + x0
      xhm1 = (xh1+xm1) / 2 + x0
      ylm1 = (yl1+ym1) / 2 + y0
      yhm1 = (yh1+ym1) / 2 + y0
      xlm2 = (xl2+xm2) / 2 - x0
      xhm2 = (xh2+xm2) / 2 - x0
      ylm2 = (yl2+ym2) / 2 - y0
      yhm2 = (yh2+ym2) / 2 - y0
      call pltxtc (xlm1+xelem, ylm1+yelem, '1')
      call pltxtc (xhm1+xelem, ylm1+yelem, '2')
      call pltxtc (xhm1+xelem, yhm1+yelem, '6')
      call pltxtc (xlm1+xelem, yhm1+yelem, '5')
      call pltxtc (xlm2+xelem, ylm2+yelem, '4')
      call pltxtc (xhm2+xelem, ylm2+yelem, '3')
      call pltxtc (xhm2+xelem, yhm2+yelem, '7')
      call pltxtc (xlm2+xelem, yhm2+yelem, '8')
      call option ('elements', .false., vcs, penlin, pentxt)

      call pltend
      end

      subroutine option (what, ison, vcs, penlin, pentxt)

      character*(*) what
      logical ison

      logical ldum, pltgtt, pltstt, pltstv, pltstd

      call pltflu

      if (ison) then
         if (what .eq. 'old lines') then
            call pltstd (1, 7.0)

         else if (what .eq. 'hidden lines') then
            call pltstd (1, 3.0)
            ldum = pltstv (2, .5*penlin) !thin line

         else if (what .eq. 'new lines') then
            call pltstd (1, 1.0)
            ldum = pltstv (1, 4.0) !short dash line

         else if (what .eq. 'new hidden lines') then
            call pltstd (1, 3.0)
            ldum = pltstv (1, 4.0) !short dash line
            ldum = pltstv (2, .5*penlin) !thin line

         else if (what .eq. 'old nodes') then
            call pltstd (1, 7.0)

         else if (what .eq. 'old spots') then
            call pltstd (1, 7.0)
            ldum = pltstt (2, 1.5*vcs) !big letters

         else if (what .eq. 'new nodes') then
            call pltstd (1, 1.0)

         else if (what .eq. 'new spots') then
            call pltstd (1, 1.0)
            ldum = pltstt (2, 0.75*vcs) !small letters

         else if (what .eq. 'elements') then
            call pltstd (1, 1.0)
            ldum = pltstt (2, 2*vcs) !big letters

         else
            call pltflu
            print *, 'unknown option ', what
         end if

      else
         call pltstd (1, 7.0)
         ldum = pltstt (2, vcs) !normal letters
         ldum = pltstv (1, 1.0) !solid line
         ldum = pltstv (2, penlin) !normal line
         ldum = pltstt (11, pentxt) !normal pen
      end if

      return
      end
