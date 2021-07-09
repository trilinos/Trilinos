C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      block data cdrabc
      integer*4 machin(3),maclen
      integer*4 kidsiz,kjobid(4),kusrsz,kusrid(4),kszrou
      integer*4 kjrout(4),ksecur,kjtime(4),kjdate(4)
      common /vcjob/ kidsiz,kjobid,kusrsz,kusrid,kszrou,
     1               kjrout,ksecur,kjtime,kjdate,machin,maclen

      integer*4 kwrtfl,krdfl,koutfl,kinfl,kwrdsz,kbytel,kcpw,kbaud,
     1kcomtp
      common /cdrcom/ kwrtfl,krdfl,koutfl,kinfl,kwrdsz,kbytel,kcpw,
     1 kbaud,kcomtp

      common /cdrcm2 / KGNAME
      character*80 KGNAME

      integer*4 kuntfd(0:999)
      common /cdrunx/ kuntfd

      logical*4 onode
      common /vconod/ onode

#if defined(CRAY)
      data kwrdsz /64/
      data kcpw /8/
#else
      data kwrdsz /32/
      data kcpw /4/
#endif

      data kwrtfl /6/
      data krdfl /5/
      data koutfl /77/
      data kinfl /66/
      data kbytel /8/
      data kbaud /9600/

c  COMPUTER TYPE WHERE DATA ORIGINATED
c  VAX NETWORK = 4
c  NOS         = 5
c  SCOPE 2.1   = 6
c  UNIVAC      = 7
c  cray        = 8
c  ALLIANT     = 9
c  APOLLO      = 10
c  SUN         = 11
c  HP          = 12
c  DEC/ULTRIX  = 13

      data kcomtp /11/

c initialize unit to file descriptor table (cdrunx)
      data kuntfd/1000 * 0/

C initialize graphics output filename to blank (cdrcm2)
      data kgname/' '/

C following initializations done just to force common block symbol
C to object file on some compilers...
      data onode /.true./
      data kidsiz /0/
      end

