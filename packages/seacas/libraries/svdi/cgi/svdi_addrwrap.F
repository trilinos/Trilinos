C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

        subroutine cgidum
        end

#ifdef  Build64

c       This wrapper module handles two kinds of interface problems:
c       the routines starting with 'v' arbitrate between Fortran caller
c       and C callee; the routines; starting with 'c' arbitrate between
c       C caller and Fortran callee.  The only exception is cgtxx1,
c       which arbitrates between Fortran caller and C callee.

        subroutine vdinit (aspect,justif)
        real*8 aspect
        integer*8 justif
        real*4 aspect4
        integer*4 justif4
        aspect4 = aspect
        justif4 = justif
        call vdinit4 (aspect4,justif4)
c       aspect = aspect4
c       justif = justif4
        end
c-----------------------------------------------------------------------
        subroutine vdiqdc (indexa,value)
        integer*8 indexa
        real*8 value
        integer*4 indexa4
        real*4 value4
        indexa4 = indexa
C ... 'value' is writeonly,
        call vdiqdc4 (indexa4,value4)
        indexa = indexa4
        value = value4
        end
c-----------------------------------------------------------------------
        subroutine vdstco (num,index_array,color_array,color_mod)
        integer*8 num
        integer*8 index_array(num)
        real*8 color_array(3,num)
        integer*8 color_mod
        integer*4 num4
        integer*4 index_array4(num)
        real*4 color_array4(3,num)
        integer*4 color_mod4
        INTEGER J1X, J2X
        num4 = num
        DO J1X = 1, NUM
          INDEX_ARRAY4(J1X) = INDEX_ARRAY(J1X)
        END DO
        DO J2X = 1, 3*NUM
          COLOR_ARRAY4(J2X,1) = COLOR_ARRAY(J2X,1)
        END DO
        color_mod4 = color_mod
        call vdstco4 (num4,index_array4,color_array4,color_mod4)
c       num = num4
c       DO J1X = 1, NUM
c         INDEX_ARRAY(J1X) = INDEX_ARRAY4(J1X)
c       END DO
c       DO J2X = 1, 3*NUM
c         COLOR_ARRAY(J2X,1) = COLOR_ARRAY4(J2X,1)
c       END DO
c       color_mod = color_mod4
        end
c-----------------------------------------------------------------------
        subroutine vdiqco (num,index_array,color_array,color_mod)
        integer*8 num
        integer*8 index_array(num)
        real*8 color_array(3,num)
        integer*8 color_mod
        integer*4 num4
        integer*4 index_array4(num)
        real*4 color_array4(3,num)
        integer*4 color_mod4
        INTEGER J1X, J2X
        num4 = num
        DO J1X = 1, NUM
          INDEX_ARRAY4(J1X) = INDEX_ARRAY(J1X)
        END DO
C ... 'color_array' is writeonly.
        color_mod4 = color_mod
        call vdiqco4 (num4,index_array4,color_array4,color_mod4)
c       num = num4
c       DO J1X = 1, NUM
c         INDEX_ARRAY(J1X) = INDEX_ARRAY4(J1X)
c       END DO
        DO J2X = 1, 3*NUM
          COLOR_ARRAY(J2X,1) = COLOR_ARRAY4(J2X,1)
        END DO
c       color_mod = color_mod4
        end
c-----------------------------------------------------------------------
        subroutine vdescp (escape_code,n,args)
        integer*8 escape_code
        integer*8 n
        real*8 args(n)
        integer*4 escape_code4
        integer*4 n4
        real*4 args4(n)
        INTEGER J1X
        escape_code4 = escape_code
        n4 = n
        DO J1X = 1, N
          ARGS4(J1X) = ARGS(J1X)
        END DO
        call vdescp4 (escape_code4,n4,args4)
c       escape_code = escape_code4
c       n = n4
        DO J1X = 1, N
          ARGS(J1X) = ARGS4(J1X)
        END DO
        end
c-----------------------------------------------------------------------
        subroutine vdiqes (escape_code,support)
        integer*8 escape_code
        integer*8 support
        integer*4 escape_code4
        integer*4 support4
C... escape_code unused, support is write only.
        call vdiqes4 (escape_code4,support4)
        support = support4
        end
c-----------------------------------------------------------------------
        subroutine vdiqnd (x_ndc,y_ndc)
        real*8 x_ndc
        real*8 y_ndc
        real*4 x_ndc4
        real*4 y_ndc4
C ... x_ndc, y_ndc write-only arguments
        call vdiqnd4 (x_ndc4,y_ndc4)
        x_ndc = x_ndc4
        y_ndc = y_ndc4
        end
c-----------------------------------------------------------------------
        subroutine vdmova (x,y)
        real*8 x
        real*8 y
        real*4 x4
        real*4 y4
        x4 = x
        y4 = y
        call vdmova4 (x4,y4)
        x = x4
        y = y4
        end
c-----------------------------------------------------------------------
        subroutine vdlina (x,y)
        real*8 x
        real*8 y
        real*4 x4
        real*4 y4
        x4 = x
        y4 = y
        call vdlina4 (x4,y4)
        x = x4
        y = y4
        end
c-----------------------------------------------------------------------
        subroutine vdpnta (x,y)
        real*8 x
        real*8 y
        real*4 x4
        real*4 y4
        x4 = x
        y4 = y
        call vdpnta4 (x4,y4)
        x = x4
        y = y4
        end
c-----------------------------------------------------------------------
        subroutine vdtext (length,char_array)
        integer*8 length
        integer*8 char_array(length)
        integer*4 length4
        integer*4 char_array4(length)
        INTEGER J1X
        length4 = length
        DO J1X = 1, LENGTH
          CHAR_ARRAY4(J1X) = CHAR_ARRAY(J1X)
        END DO
        call vdtext4 (length4,char_array4)
        length = length4
        DO J1X = 1, LENGTH
          CHAR_ARRAY(J1X) = CHAR_ARRAY4(J1X)
        END DO
        end
c-----------------------------------------------------------------------
        subroutine vdpoly (x_array,y_array,npts)
        integer*8 npts
        real*8 x_array(npts)
        real*8 y_array(npts)
        real*4 x_array4(npts)
        real*4 y_array4(npts)
        integer*4 npts4
        INTEGER J1X
        DO J1X = 1, NPTS
          X_ARRAY4(J1X) = X_ARRAY(J1X)
          Y_ARRAY4(J1X) = Y_ARRAY(J1X)
        END DO
        npts4 = npts
        call vdpoly4 (x_array4,y_array4,npts4)
        DO J1X = 1, NPTS
          X_ARRAY(J1X) = X_ARRAY4(J1X)
          Y_ARRAY(J1X) = Y_ARRAY4(J1X)
        END DO
        npts = npts4
        end
c-----------------------------------------------------------------------
        subroutine vdiqcp (x,y)
        real*8 x
        real*8 y
        real*4 x4
        real*4 y4
C... 'x' and 'y' are write-only arguments
        call vdiqcp4 (x4,y4)
        x = x4
        y = y4
        end
c-----------------------------------------------------------------------
        subroutine vdiqos (attr_array)
        real*8 attr_array(7)
        real*4 attr_array4(7)
        INTEGER J1X
C ... attr_array4 argument is write only, don't initialize
        call vdiqos4 (attr_array4)
        DO J1X = 1, 7
          ATTR_ARRAY(J1X) = ATTR_ARRAY4(J1X)
        END DO
        end
c-----------------------------------------------------------------------
        subroutine vdstfc (color_index)
        integer*8 color_index
        integer*4 color_index4
        color_index4 = color_index
        call vdstfc4 (color_index4)
c       color_index = color_index4
        end
c-----------------------------------------------------------------------
        subroutine vdstbc (color_index)
        integer*8 color_index
        integer*4 color_index4
        color_index4 = color_index
        call vdstbc4 (color_index4)
c       color_index = color_index4
        end
c-----------------------------------------------------------------------
        subroutine vdstin (intensity)
        real*8 intensity
        real*4 intensity4
        intensity4 = intensity
        call vdstin4 (intensity4)
c       intensity = intensity4
        end
c-----------------------------------------------------------------------
        subroutine vdstls (line_style)
        integer*8 line_style
        integer*4 line_style4
        line_style4 = line_style
        call vdstls4 (line_style4)
c       line_style = line_style4
        end
c-----------------------------------------------------------------------
        subroutine vdstlw (line_width)
        real*8 line_width
        real*4 line_width4
        line_width4 = line_width
        call vdstlw4 (line_width4)
c       line_width = line_width4
        end
c-----------------------------------------------------------------------
        subroutine vdstcs (y_size)
        real*8 y_size
        real*4 y_size4
        y_size4 = y_size
        call vdstcs4 (y_size4)
c       y_size = y_size4
        end
c-----------------------------------------------------------------------
        subroutine vdaabu (button)
        integer*8 button
        integer*4 button4
C ... 'button' is write-only argument
        call vdaabu4 (button4)
        button = button4
        end
c-----------------------------------------------------------------------
        subroutine vdaloc (x,y)
        real*8 x
        real*8 y
        real*4 x4
        real*4 y4
C ... 'x' and 'y' are write-only arguments
        call vdaloc4 (x4,y4)
        x = x4
        y = y4
        end
c-----------------------------------------------------------------------
        subroutine vdabgl (button,x,y)
        integer*8 button
        real*8 x
        real*8 y
        integer*4 button4
        real*4 x4
        real*4 y4
C... 'button', 'x', and 'y' are write-only
        call vdabgl4 (button4,x4,y4)
        button = button4
        x = x4
        y = y4
        end
c-----------------------------------------------------------------------
        subroutine vdakgl (charac,x,y)
        integer*8 charac
        real*8 x
        real*8 y
        integer*4 charac4
        real*4 x4
        real*4 y4
C ... arguments are write-only
        charac4 = 0
        x4 = 0.0
        y4 = 0.0
        call vdakgl4 (charac4,x4,y4)
        charac = charac4
        x = x4
        y = y4
        end
c-----------------------------------------------------------------------
        subroutine vdstla (x,y)
        real*8 x
        real*8 y
        real*4 x4
        real*4 y4
C ... arguments are write-only
        x4 = 0.0
        y4 = 0.0
        call vdstla4 (x4,y4)
        x = x4
        y = y4
        end
c-----------------------------------------------------------------------
        subroutine vdstos (attr_array)
        real*8 attr_array(6)
        real*4 attr_array4(6)
        INTEGER J1X
        DO J1X = 1, 6
          ATTR_ARRAY4(J1X) = ATTR_ARRAY(J1X)
        END DO
        call vdstos4 (attr_array4)
        DO J1X = 1, 6
          ATTR_ARRAY(J1X) = ATTR_ARRAY4(J1X)
        END DO
        end
c-----------------------------------------------------------------------
        subroutine vdfram (itype)
        integer*8 itype
        integer*4 itype4
        itype4 = itype
        call vdfram4 (itype4)
c       itype = itype4
        end
c-----------------------------------------------------------------------
        subroutine cesc24 (i1, i2, r1)
        integer*4 i1, i2
        real*4 r1
        integer*8 l1, l2
        real*8 d1
        l1 = i1
        l2 = i2
        d1 = r1
        call cesc2 (l1, l2, d1)
        i1 = l1
        i2 = l2
        r1 = d1
        end
c-----------------------------------------------------------------------
        subroutine ctx24 (r1, r2, i1, i2)
        integer*4 i1(150), i2
        real*4 r1, r2
        integer*8 l1(150), l2
        real*8 d1, d2
        INTEGER J1X
        DO J1X = 1, 150
          L1(J1X) = I1(J1X)
        END DO
        l2 = i2
        d1 = r1
        d2 = r2
        call ctx2 (d1, d2, l1, l2)
        DO J1X = 1, 150
          I1(J1X) = L1(J1X)
        END DO
        i2 = l2
        r1 = d1
        r2 = d2
        end
c-----------------------------------------------------------------------
        subroutine cgtxx24 (r1, r2, i1, i2, r3, r4, r5, r6, r7,
     &          r8, r9, r10, r11, r12)
        integer*4 i1, i2
        real*4 r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12
        integer*8 l1, l2
        real*8 d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12
        l1 = i1
        l2 = i2
        d1 = r1
        d2 = r2
        d3 = r3
        d4 = r4
        d5 = r5
        d6 = r6
        d7 = r7
        d8 = r8
        d9 = r9
        d10 = r10
        d11 = r11
        d12 = r12
        call cgtxx2 (d1, d2, l1, l2, d3, d4, d5, d6, d7,
     &          d8, d9, d10, d11, d12)
        i1 = l1
        i2 = l2
        r1 = d1
        r2 = d2
        r3 = d3
        r4 = d4
        r5 = d5
        r6 = d6
        r7 = d7
        r8 = d8
        r9 = d9
        r10 = d10
        r11 = d11
        r12 = d12
        end
c-----------------------------------------------------------------------
        subroutine cqchh24 (i1, i2, i3, i4, i5, i6, i7)
        integer*4 i1, i2, i3, i4, i5, i6, i7
        integer*8 l1, l2, l3, l4, l5, l6, l7
        l1 = i1
        l2 = i2
        l3 = i3
        l4 = i4
        l5 = i5
        l6 = i6
        l7 = i7
        call cqchh2 (l1, l2, l3, l4, l5, l6, l7)
c       i1 = l1
c       i2 = l2
c       i3 = l3
        i4 = l4
        i5 = l5
        i6 = l6
        i7 = l7
        end
c-----------------------------------------------------------------------
        subroutine cesc1 (l1, l2, s, l3)
        integer*8 l1, l2, l3
        character*(*) s
        integer*4 i1, i2, i3
        i1 = l1
        i2 = l2
        i3 = l3
        call cesc14 (i1, i2, s, i3)
        l1 = i1
        l2 = i2
        l3 = i3
        end
c-----------------------------------------------------------------------
        subroutine ctx1 (d1, d2, l1, s, l2)
        integer*8 l1, l2
        real*8 d1, d2
        character*(*) s
        integer*4 i1, i2
        real*4 r1, r2
        i1 = l1
        i2 = l2
        r1 = d1
        r2 = d2
        call ctx14( r1, r2, i1, s, i2)
c       l1 = i1
        l2 = i2
c       d1 = r1
c       d2 = r2
        end
c-----------------------------------------------------------------------
        subroutine cgtxx1 (d1, d2, s, l1, l2, d3, d4, d5, d6, d7,
     &          d8, d9, d10, d11, d12)
        integer*8 l1, l2
        real*8 d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11, d12
        character*(*) s
        integer*4 i1, i2, i3
        real*4 r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12
        i1 = l1
        i2 = l2
        r1 = d1
        r2 = d2
        r3 = d3
        r4 = d4
        r5 = d5
        r6 = d6
        r7 = d7
        r8 = d8
        r9 = d9
        r10 = d10
        r11 = d11
        r12 = d12
        i3 = len(s)
        call cgtxx14 (r1, r2, s, i1, i2, r3, r4, r5, r6, r7,
     &          r8, r9, r10, r11, r12, i3)
        l1 = i1
        l2 = i2
c       d1 = r1
c       d2 = r2
        d3 = r3
        d4 = r4
        d5 = r5
        d6 = r6
        d7 = r7
        d8 = r8
        d9 = r9
        d10 = r10
        d11 = r11
        d12 = r12
        end
c-----------------------------------------------------------------------
        subroutine cqchh1 (s, l1, l2, l3, l4, l5, l6, l7, l8)
        integer*8 l1, l2, l3, l4, l5, l6, l7, l8
        character*(*) s
        integer*4 i1, i2, i3, i4, i5, i6, i7, i8
        i1 = l1
        i2 = l2
        i3 = l3
        i4 = l4
        i5 = l5
        i6 = l6
        i7 = l7
        i8 = l8
        call cqchh14 (s, i1, i2, i3, i4, i5, i6, i7, i8)
C       l1 = i1
C       l2 = i2
C       l3 = i3
        l4 = i4
        l5 = i5
        l6 = i6
        l7 = i7
        l8 = i8
        end
c-----------------------------------------------------------------------
#endif
