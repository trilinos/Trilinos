/*
 * Copyright (C) 2009-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
/* cgidef.h   CGI enumerated type definitions for C programs
 * 8 Sep 1989, last date modified
 * Pat McGee, jpm@lanl.gov
 */

/* escape enumerations */
#define XEMFNM -28372 /* metafile name */
#define XEAGMD -23671 /* alpha/graphics mode */
#define XESVDI -1001  /* svdi escape  */
#define XEPGSZ -1002  /* page size */
#define XEBELL -1003  /* ring the bell */

/* used in many places */
#define CNO 0
#define CYES 1

/* force clear viewsurface (argument of CPDS) */
#define CFORCC 0
#define CCONDC 1

/* clip indicator */
#define COFF 0
#define CON 1

/* drawing surface clip indicator */
#define CDCOFF 0
#define CDCREC 1
#define CVPORT 2

/* error handling flag (arg of CERHCT) */
#define CEHON 0
#define CEHROF 1
#define CEHDOF 2

/* device class (arg of CQID) */
#define COUTPT 0
#define CINPUT 1
#define COUTIN 2

/* hard/soft copy flag (arg of CQD) */
#define CHARD 0
#define CSOFT 1

/* display type (arg of CQD) */
#define CVECT 0
#define CRAST 1
#define COTHER 2

/* dynamic modification (arg of CQD, etc) */
#define CIRG 0
#define CCBS 1
#define CIMM 2

/* pixel location relative to coordinates (arg of CQD) */
#define CPXON 0
#define CPXBET 1

/* support indicator */
#define CSNO 0
#define CSYES 1
#define CSUNR 2

/* view surface state */
#define CDIRTY 0
#define CCLEAN 1

/* text final flag */
#define CNOTFI 0
#define CFINAL 1

/* clipping mode */
#define CLOCUS 0
#define CSHAPE 1
#define CLOCSH 2

/* text precision */
#define CSTRNG 0
#define CCHAR 1
#define CSTROK 2

/* text path */
#define CTPRIT 0
#define CTPLFT 1
#define CTPUP 2
#define CTPDWN 3

/* text horizonal alignment */
#define CTALFT 0
#define CTACTR 1
#define CTARIT 2
#define CTANH 3
#define CTACOH 4

/* text vertical alignment */
#define CTATOP 0
#define CTACAP 1
#define CTAHAF 2
#define CTABAS 3
#define CTABOT 4
#define CTANV 5
#define CTACOV 6

/* interior style */
#define CHOLLO 0
#define CSOLID 1
#define CPATRN 2
#define CHATCH 3
#define CEMPTY 4
#define CBITMP 5

/* color selection mode (arg of CCSM, CQTXA) */
#define CDRECT 0
#define CINDEX 1

/* line, edge width specification mode */
/* marker specification mode */
#define CVDC 0
#define CSCA 1

/* cell array fill capability (arg of CQPRL) */
#define COUTLN 0
#define CFILLD 1

/* cell array alignment (arg of CQPRL) */
#define CAXIS 0
#define CSKEW 1

/* compound text capability */
/* and closed figure capability */
#define CCNONE 0
#define CGLOBL 1
#define CLOCAL 2

/* pattern transformation support */
#define CPTNO 0
#define CPTUNS 1
#define CPTFUL 2

/* color selection mode availability (arg of CQC) */
#define CCLRI 0
#define CCLRID 1

/* color overwrite capability (arg of CQC) */
#define CCSUB 0
#define CCADD 1
#define CCREP 2
#define CCDMOD 3

/* response validity (arg of most CQxxxx routines) */
#define CINVAL 0
#define CVAL 1

/* input class */
#define CLOCAT 0
/* #define CSTROK 1 */
#define CVALUA 2
#define CCHOIC 3
#define CPICK 4
/* #define CSTRNG 5 */
#define CRASTR 6
#define CGEN 7

/* request status */
#define CTRIGR 0
#define CBREAK 1
#define CTIME 2
#define CMEASC 3

/* input device state */
#define CRELEA 0
#define CREADY 1
#define CREQPD 2
#define CERQON 3
#define CERQPD 4
#define CERQCM 5
#define CEVTON 6

/* direction (arg of CPXA) */
#define CINCR 0
#define CDECR 1

/* action required flag */
#define CNOACT 0
#define CACT 1

/* pixel validity flag (arg of CGPXA) */
#define CVNONE 0
#define CVALL 1
#define CVSOME 2

/* end cgidef.h */
