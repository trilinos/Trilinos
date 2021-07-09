C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

c vdicgi_char - FORTRAN shell for handling strings

C  CESC - Escape
      SUBROUTINE CESC2 (FUNCID, N, ARGS)
      INTEGER FUNCID
      REAL ARGS(N+1)
      INTEGER LDR
      CHARACTER*80 DATA
c  CGI enumerated type definitions for FORTRAN programs
c  8 Sep 1989, last date modified
c  Pat McGee, jpm@lanl.gov

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

      LDR = 1
      DATA = ' '
      IF (FUNCID .EQ. XEAGMD) THEN
        IF (ARGS(1) .EQ. 0) THEN
            DATA(1:5) = 'ALPHA'
        ELSE
            DATA(1:8) = 'GRAPHICS'
        ENDIF
      ELSE
c*** I know this won't work for real, real numbers
        DO 10 I = 1,N
          DATA(I:I) = CHAR(INT(ARGS(I)))
  10    CONTINUE
      ENDIF

      CALL CESC( FUNCID, LDR, DATA )
      RETURN
      END

C  CTX - Text
      SUBROUTINE CTX2(X, Y, TEXT1, LENGTH )
      REAL X,Y
      INTEGER FLAG
      INTEGER TEXT1(LENGTH)
      CHARACTER*136 TEXT2

      FLAG = 1

      DO 20 I=1,LENGTH
        TEXT2(I:I) = CHAR(TEXT1(I))
  20  CONTINUE

      CALL CTX( X, Y, FLAG, TEXT2(1:LENGTH))
      RETURN
      END

C  CGTXX - Get Text Extent
      SUBROUTINE CGTXX2 ( X, Y, VSTAT, VCONC, XCONC, YCONC,
     1                    X1, Y1, X2, Y2, X3, Y3, X4, Y4)
      REAL X,Y
      CHARACTER*1 STRING
      INTEGER VSTAT, VCONC
      REAL XCONC, YCONC
      REAL X1, Y1, X2, Y2, X3, Y3, X4, Y4

      STRING = 'M'
      CALL CGTXX( X, Y, STRING, VSTAT, VCONC, XCONC, YCONC,
     1            X1, Y1, X2, Y2, X3, Y3, X4, Y4 )
      RETURN
      END

C  CQCHH - Inquire Character Height
      SUBROUTINE CQCHH2 (txp, nreq, first, vstat, ntotal,
     1                   nlist, chhit)
      CHARACTER*1 FONT
      INTEGER TXP, NREQ, FIRST, VSTAT, NTOTAL, NLIST
      INTEGER CHHIT(*)

      FONT = char(0)
      CALL CQCHH( FONT, TXP, NREQ, FIRST, VSTAT, NTOTAL,
     1             NLIST, CHHIT)
      RETURN
      END

C  vdgnam -
      subroutine vdgnam(name)
      character *(*) name
c  CGI enumerated type definitions for FORTRAN programs
c  8 Sep 1989, last date modified
c  Pat McGee, jpm@lanl.gov

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

      call cesc(XEMFNM,1,name)
      return
      end

