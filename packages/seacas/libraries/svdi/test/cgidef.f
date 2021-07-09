C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

c  CGI enumerated type definitions for FORTRAN programs
c  8 Sep 1989, last date modified
c  Pat McGee, jpm@lanl.gov
c  Tricia Crotty, plcrott@sandia.gov

c  SRCP escapes
      integer XEMFNM, XEMXCL, XEPCTL, XEAGMD, XEPCCL, XESVDI
      parameter (XEMFNM= -28372, XEMXCL= -19281, XEPCTL= -190,
     *           XEAGMD= -23671, XEPCCL= -12048, XESVDI= -1001)

c  SRCP definitions
c  maximum error class
      integer XMXERR
      parameter (XMXERR = 9)

c  maximum function identifier
      integer XMXFCL
      parameter (XMXFCL = 61)

c  CGI definitions
c  (used in many places)
      integer CNO, CYES
      parameter (CNO= 0, CYES= 1)

c  force clear viewsurface (argument of CPDS)
      integer CFORCC, CCONDC
      parameter (CFORCC= 0, CCONDC= 1)

c  view surface state (arg of most CQxxxx routines)
      integer CDIRTY, CCLEAN
      parameter (CDIRTY=0, CCLEAN=1)

c  clip indicator
      integer COFF, CON
      parameter (COFF=0, CON=1)

c  drawing surface clip indicator
      integer CDCOFF, CDCREC, CVPORT
      parameter (CDCOFF=0, CDCREC=1, CVPORT=2)

c  error handling flag (arg of CERHCT)
      integer CEHON, CEHROF, CEHDOF
      parameter (CEHON=0, CEHROF=1, CEHDOF=2)

c  device class (arg of CQID)
      integer COUTPT, CINPUT, COUTIN
      parameter (COUTPT=0, CINPUT=1, COUTIN=2)

c  hard/soft copy flag
      integer CHARD, CSOFT
      parameter (CHARD= 0, CSOFT= 1)

c  display type
      integer CVECT, CRAST, COTHER
      parameter (CVECT= 0, CRAST= 1, COTHER= 2)

c  dynamic modification
      integer CIRG, CCBS, CIMM
      parameter (CIRG= 0, CCBS= 1, CIMM= 2)

c  pixel location relative to coordinates */
      integer CPXON, CPXBET
      parameter (CPXON=0, CPXBET=1)

c  support indicator
      integer CSNO, CSYES, CSUNR
      parameter (CSNO= 0, CSYES= 1, CSUNR=2)

c  text final flag
      integer CNOTFI, CFINAL
      parameter (CNOTFI=0, CFINAL=1)

c  clipping mode
      integer CLOCUS, CSHAPE, CLOCSH
      parameter(CLOCUS=0, CSHAPE=1, CLOCSH=2)

c  text precision
      integer CSTRNG, CCHAR, CSTROK
      parameter(CSTRNG=0, CCHAR=1, CSTROK=2)

c  text path
      integer CTPRIT, CTPLFT, CTPUP, CTPDWN
      parameter(CTPRIT=0, CTPLFT=1, CTPUP=2, CTPDWN=3)

c  text horizontal alignment
      integer CTALFT, CTACTR, CTARIT, CTANH, CTACOH
      parameter(CTALFT=0, CTACTR=1, CTARIT=2, CTANH=3, CTACOH=4)

c  text vertical alignment
      integer CTATOP, CTACAP, CTAHAF, CTABAS
      integer CTABOT, CTANV, CTACOV
      parameter(CTATOP=0, CTACAP=1, CTAHAF=2, CTABAS=3)
      parameter(CTABOT=4, CTANV=5, CTACOV=6)

c  interior style
      integer CHOLLO, CSOLID
      parameter (CHOLLO= 0, CSOLID= 1)

c  color selection mode (arg of CCSM, CQTXA)
      integer CDRECT, CINDEX
      parameter (CDRECT= 0, CINDEX= 1)

c  line, edge width specification mode
c  marker specification mode

c  cell array fill capability
      integer COUTLN, CFILLD
      parameter (COUTLN=0, CFILLD=1)

c  cell array alignment
      integer CAXIS, CSKEW
      parameter (CAXIS=0, CSKEW=1)

c  compound text capability
c  and closed figure capability
      integer CCNONE, CGLOBL, CLOCAL
      parameter (CCNONE=0, CGLOBL=1, CLOCAL=2)

c  pattern transformation support

c  color selection mode availability
      integer CCLRI, CCLRID
      parameter (CCLRI=0, CCLRID=1)

c  color overwrite capability

c  response validity
      integer CINVAL, CVAL
      parameter (CINVAL= 0, CVAL= 1)

c  input class

c  request status

c  input device state

c  direction
      integer CINCR, CDECR
      parameter (CINCR= 0, CDECR= 1)

c  action required flag
      integer CNOACT, CACT
      parameter (CNOACT= 0, CACT= 1)

c  pixel validity flag
      integer CVNONE, CVALL, CVSOME
      parameter (CVNONE= 0, CVALL= 1, CVSOME= 2)
