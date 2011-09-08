C Copyright (c) 2008 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++ Copyright 1988, Sandia Corporation. The United States Government
C+++ retains a limited license in this software as prescribed in AL 88-1
C+++ and AL 91-7. Export of this program may require a license from
C+++ the United States Government.
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C -*- Mode: fortran -*-
C=======================================================================
C $Id: gjoin2.f,v 1.10 2008/07/31 20:15:56 gdsjaar Exp $
c 
      PROGRAM GJOIN2
C=======================================================================

C                         *** GJOIN2 ***

C   --*** GJOIN2 *** (GJOIN) GENESIS database combination program
C   --   Written by Amy Gilkey - revised 03/04/88
C   --
C   --GJOIN2 combines two or more GENESIS II databases into a single database.
C   --
C   --Expected input:
C   --   o Responses from the user on the standard input device.
C   --   o The input databases (name requested)
C   --
C   --Output:
C   --   o Prompts on the standard output device.
C   --   o The output database (name requested)
C
C   --Developed at Sandia National Laboratories.
C   --
C   --Current author and code sponsor: Greg Sjaardema
C   --
C   --Source is in FORTRAN 77
C   --
C   --External software used:
C   --   SUPES package (dynamic memory, free-field reader, FORTRAN extensions)
C   --
C   --Documentation:
C   --   "User's Manual for GJOIN"

      include 'exodusII.inc'

      include 'params.blk'
      INCLUDE 'progqa.blk'
      INCLUDE 'titles.blk'
      INCLUDE 'dbvars.blk'
      INCLUDE 'filnum.blk'
      INCLUDE 'xyzrot.blk'

      DIMENSION A(1), IA(1)
C      --A - the dynamic numeric memory base array
      EQUIVALENCE (A(1),IA(1))
      character*1 c(1)

      LOGICAL USESDF, NONQUD
      LOGICAL RENNP, RENEL, REN, DELNP, DELEL, BATCH, CLOSE, MATMAT
      LOGICAL FIRST, DONE, MDEBUG
      character*(256) filnam, string

C... String containing name of common element topology in model
C    or 'MULTIPLE_TOPOLOGIES' if not common topology.
      character*(MXSTLN) comtop

      character*(MXSTLN) qarec(4,MAXQA)
      character*(MXLNLN) infrec(MAXINF)
C      --QAREC - the QA records
C      --INFREC - the information records

      INCLUDE 'qainfo.blk'

      CALL STRTUP (QAINFO)

      CALL BANNER (0, QAINFO,
     &   'A GENESIS DATABASE COMBINATION PROGRAM',
     &   ' ', ' ')

      call cpyrgt (0, '1988')

      call exinq (netid, EXLBVR, idummy, exlibversion, name, nerr)
      write(*,'(A,F6.3)')'ExodusII Library version ',
     &          exlibversion

C     --Open LOG file if running interactively
      IF (.NOT. BATCH()) THEN
         KLOG = 99
         OPEN (UNIT=KLOG, FILE='gjoin.log', FORM='formatted',
     &        STATUS='unknown', IOSTAT=IERR)
         IF (IERR .NE. 0) THEN
            CALL PRTERR ('ERROR', 'Could not open log file')
            GOTO 150
         END IF
      ELSE
         KLOG = 0
      END IF

C   --Initialize dynamic memory

      NQAREC = 0
      NINFO  = 0
      CALL MDINIT (A)
      CALL MCINIT (C)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 140

      MDEBUG = .false.
C      if (MDEBUG) then
C         call mlist()
C      end if

      FIRST = .TRUE.
      USESDF = .FALSE.

 80   CONTINUE

      CALL INIGEN (A, FIRST,
     &   KXN, KYN, KZN, KMAPEL,
     &   KIDELB, KNELB, KNLNK, KNATR, KLINK, KATRIB,
     &   KIDNS, KNNNS, KIXNNS, KLTNNS, KFACNS, 
     &   KIDSS, KNESS, KNDSS, KIXESS, KIXDSS, KLTESS, kltsss,
     &   kltsnc, kfacss, KNMLB)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 140

C   --Open and read the database

      WRITE (*, *)

      IF (FIRST) THEN
 90      CONTINUE
         CALL GETINP (0, 0, 'First input file> ', FILNAM, IOSTAT)
         IF (IOSTAT .LT. 0) GOTO 150
         IF (FILNAM(1:1) .EQ. '$' .OR. LENSTR(FILNAM) .LE. 1)  GOTO 90
         CALL OUTLOG (KLOG, 1, 0, FILNAM, IDUM, RDUM)
      END IF

      CALL RDGEN (A, IA, C, FIRST, FILNAM,
     &  TITLE1, NDIM, NUMNP1, NUMEL1, NELBL1,
     &  NNPS1, LNPSN1, NESS1, LESSE1, LESSD1,
     &  KXN, KYN, KZN, KMAPEL,
     &  KIDELB, KNELB, KNLNK, KNATR, KLINK, KATRIB,
     &  KIDNS, KNNNS, KIXNNS, KLTNNS, KFACNS,
     &  KIDSS, KNESS, KNDSS, KIXESS, KIXDSS, KLTESS, kltsss,
     &  kltsnc, kfacss, NQAREC, QAREC, NINFO, INFREC, KNMLB,
     &  USESDF, *150)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 140

 100  CONTINUE
C   --Save the length of the LINK array for later
      CALL MDFIND ('LINK', IDUM, LLINK1)

C   --Open and read the second database

      WRITE (*, *)

  110 CONTINUE
      CALL GETINP (0, 0, 'Next input file> ', FILNAM, IOSTAT)
      IF (IOSTAT .LT. 0) GOTO 150
      CALL OUTLOG (KLOG, 1, 0, FILNAM, IDUM, RDUM)
      TWODB = (FILNAM .NE. ' ')


      IF (TWODB) THEN
        CALL RDGEN (A, IA, C, .TRUE., FILNAM,
     &    TITLE2, NDIM2, NUMNP2, NUMEL2, NELBL2,
     &    NNPS2, LNPSN2, NESS2, LESSE2, LESSD2, 
     &    KXN, KYN, KZN, KMAPEL,
     &    KIDELB, KNELB, KNLNK, KNATR, KLINK, KATRIB,
     &    KIDNS, KNNNS, KIXNNS, KLTNNS, KFACNS,
     &    KIDSS, KNESS, KNDSS, KIXESS, KIXDSS, KLTESS, kltsss,
     &    kltsnc, kfacss, NQAREC, QAREC, NINFO, INFREC, KNMLB,
     &    USESDF, *110)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 140
        
        IF (NDIM .NE. NDIM2) THEN
          CALL PRTERR ('FATAL', 'Number of dimensions must match')
          GOTO 110
        END IF

        call CHKTOP(NELBL2, C(KNMLB), COMTOP)
      ELSE
        TITLE2 = ' '
        NUMNP2 = 0
        numel2 = 0
        nelbl2 = 0
        nnps2  = 0
        lnpsn2 = 0
        ness2  = 0
        lesse2 = 0
        comtop = 'MULTIPLE_TOPOLOGIES'
      END IF

C   --Add an offset to the nodal point and element side set pointers
C   --for the second database

      IF (TWODB) THEN
         CALL RENIX (NNPS2, LNPSN1, -999, IA(KIXNNS+NNPS1))
         CALL RENIX (NESS2, LESSE1, -999, IA(KIXESS+NESS1))
         CALL RENIX (NESS2, LESSD1, -999, IA(KIXDSS+NESS1))
      END IF

C   --Combine the nodes

      NEWNP = NUMNP1 + NUMNP2

      IF (TWODB) THEN
         NMAT = 0

  120    CONTINUE
         xoff = 0.0
         yoff = 0.0
         zoff = 0.0

         xscl = 1.0
         yscl = 1.0
         zscl = 1.0
         rot3d = .false.

         CALL IRENNP (A, NNPS1, NNPS2, IA(KIDNS), IA(KNNNS),
     &        REN, MATNS1, MATNS2, TOLER, CLOSE, MATMAT,
     $        XSCL, YSCL, ZSCL, XOFF, YOFF, ZOFF, NDIM, IEXPCT)

         IF (ROT3D) call dorot(ndim, numnp2, a(kxn+numnp1),
     $        a(kyn+numnp1), a(kzn+numnp1),
     $        rotmat, rotcen)

         if (xoff .ne. 0.0 .or. xscl .ne. 1.0) then
            call offset( xoff, xscl, a(kxn+numnp1), numnp2)
            write (string, 9000) 'X', xscl, 'X', xoff
            call sqzstr(string, lstr)
            call prterr('CMDSPEC', string(:lstr))
         end if

         if (yoff .ne. 0.0 .or. yscl .ne. 1.0) then
            call offset( yoff, yscl, a(kyn+numnp1), numnp2)
            write (string, 9000) 'Y', yscl, 'Y', yoff
            call sqzstr(string, lstr)
            call prterr('CMDSPEC', string(:lstr))
         end if

         if ((zoff .ne. 0.0 .or. zscl .ne. 1.0) .and. ndim .eq. 3) then
            call offset( zoff, zscl, a(kzn+numnp1), numnp2)
            write (string, 9000) 'Z', zscl, 'Z', zoff
            call sqzstr(string, lstr)
            call prterr('CMDSPEC', string(:lstr))
         else
            zoff = 0.0
            zscl = 1.0
         end if


         IF (XSCL * YSCL * ZSCL .LT. 0.0) THEN
            kidel2 = kidelb + nelbl1
            knelb2 = knelb  + nelbl1
            knlnk2 = knlnk  + nelbl1
            klink2 = klink  + llink1
            knmlb2 = knmlb  + MXSTLN*nelbl1
            
            CALL DBMIRR (1, NELBL2, IA(KIDEL2), IA(KNELB2), IA(KNLNK2),
     *        IA(KLINK2), C(KNMLB2), NDIM, NONQUD)

C ... Note that at this point, the index arrays have already been offset
C    for the second database, so the arrays containing lists of
C    nodes/dist-fact are passed in with no offset.
            CALL MIRSS (NESS2, LESSE2, LESSD2, 
     *        IA(KIDSS+NESS1), IA(KNESS+NESS1), IA(KNDSS+NESS1),
     *        IA(KIXESS+NESS1), IA(KIXDSS+NESS1), IA(KLTESS),
     *        IA(KLTSSS), IA(KLTSNC), A(KFACSS),
     *        USESDF, NONQUD, COMTOP)
         END IF

         IF (REN) THEN
           IF (NMAT .LE. 0) THEN
             CALL MDRSRV ('IXNP2', KIXNP2, NUMNP2)
             CALL MDSTAT (NERR, MEM)
             IF (NERR .GT. 0) GOTO 140
           END IF
           
           CALL MDRSRV ('IX1', KIX1, NUMNP1)
           CALL MDRSRV ('IX2', KIX2, NUMNP2)
           CALL MDSTAT (NERR, MEM)
           IF (NERR .GT. 0) GOTO 140
           
           if (.not. MATMAT) then
             CALL MATXYZ (NDIM,
     &         MATNS1, MATNS2, IA(KNNNS), IA(KIXNNS), IA(KLTNNS),
     &         NUMNP1, A(KXN),        A(KYN),        A(KZN),
     &         NUMNP2, A(KXN+NUMNP1), A(KYN+NUMNP1), A(KZN+NUMNP1),
     &         IA(KIX1), IA(KIX2), IA(KIXNP2), NMAT, TOLER, CLOSE,
     $         IEXPCT)
           else
             kidel2 = kidelb + nelbl1
             knelb2 = knelb  + nelbl1
             knlnk2 = knlnk  + nelbl1
             klink2 = klink  + llink1
             
             CALL MDRSRV ('IX3', KIX3, NUMNP1)
             CALL MDRSRV ('IX4', KIX4, NUMNP2)
             CALL MDSTAT (NERR, MEM)
             CALL EXPXYZ (NDIM,
     &         MATNS1, MATNS2, IA(KNNNS), IA(KIXNNS), IA(KLTNNS),
     &         NUMNP1, A(KXN),        A(KYN),        A(KZN),
     &         NUMNP2, A(KXN+NUMNP1), A(KYN+NUMNP1), A(KZN+NUMNP1),
     &         IA(KIX1), IA(KIX2), IA(KIX3), IA(KIX4), IA(KIXNP2),
     $         NELBL1, IA(KIDELB), IA(KNELB),  IA(KNLNK),  IA(KLINK),
     $         NELBL2, IA(KIDEL2), IA(KNELB2), IA(KNLNK2), A(KLINK2),
     $         NMAT, TOLER, CLOSE, MATMAT)
             CALL MDDEL ('IX3')
             CALL MDDEL ('IX4')
           end if
           
           IF (NMAT .LE. 0) CALL MDDEL ('IXNP2')
           CALL MDDEL ('IX1')
           CALL MDDEL ('IX2')
           
           GOTO 120
         END IF
         
         RENNP = (NMAT .GT. 0)
         IF (.NOT. RENNP) KIXNP2 = 1
         
C      --"Munch" the coordinates (second set)
         
         IF (RENNP) THEN
           NEWJNP = NUMNP2
           CALL MUNXYZ (NDIM, NEWJNP, NUMNP1, IA(KIXNP2),
     &       A(KXN+NUMNP1), A(KYN+NUMNP1), A(KZN+NUMNP1))
           
           NEWNP = NUMNP1 + NEWJNP
           
           CALL MDLONG ('XN', KXN, NEWNP)
           CALL MDLONG ('YN', KYN, NEWNP)
           CALL MDLONG ('ZN', KZN, NEWNP)
         END IF
         
         IF (RENNP) THEN
           IOFFNP = -999
         ELSE
           IOFFNP = NUMNP1
         END IF
         
C      --Renumber the element block nodes (second set)
         
         CALL RENELB (NELBL2, IOFFNP, IA(KIXNP2),
     &     IA(KNELB+NELBL1), IA(KNLNK+NELBL1), IA(KLINK+LLINK1))
         
C      --Renumber the nodal point set nodes (second set)
         
         CALL RENIX (LNPSN2, IOFFNP, IA(KIXNP2), IA(KLTNNS+LNPSN1))
         
         IF (RENNP) CALL MDDEL ('IXNP2')
      END IF
C     End of IF(TWODB)      
C   --Initialize items for output database

      TITLE = TITLE1

      NEWELB = NELBL1 + NELBL2
      NEWEL = NUMEL1 + NUMEL2
      
      NEWNPS = NNPS1 + NNPS2
      NEWNNL = INTADD (NNPS1+NNPS2, IA(KNNNS))

      NEWESS = NESS1 + NESS2
      NEWSEL = INTADD (NESS1+NESS2, IA(KNESS))
      NEWSDL = INTADD (NESS1+NESS2, IA(KNDSS))

C   --Set up status arrays for user manipulation of element blocks and sets

      CALL MDRSRV ('IELBST', KIELBS, NEWELB)
      CALL MDRSRV ('INPSST', KINPSS, NEWNPS)
      CALL MDRSRV ('IESSST', KIESSS, NEWESS)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 140

C   --Allow user to change element blocks and sets

      CALL COMAND (A,
     &   IA(KIELBS), IA(KIDELB), IA(KNELB), IA(KNLNK), IA(KNATR),
     &   C(KNMLB), IA(KINPSS), IA(KIDNS), IA(KNNNS),
     &   IA(KIESSS), IA(KIDSS), IA(KNESS), DONE, *150)

C   --"Munch" the element blocks

      I = INTCNT (0, IA(KIELBS), NEWELB)
      RENEL = (I .LT. NEWELB)

      IF (RENEL) THEN
         CALL MDRSRV ('IXEL', KIXEL, NEWEL)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 140

         CALL MDFIND ('LINK', IDUM, LLNK)
         CALL MDRSRV ('LINKO', KLINKO, LLNK)
         CALL MDFIND ('ATRIB', IDUM, LATR)
         CALL MDRSRV ('ATRIBO', KATRO, LATR)
         CALL MDRSRV ('IXELB', KIXELB, NEWELB)
         CALL MDRSRV ('JNELB', KJNELB, NEWELB)
         CALL MDRSRV ('ISCR', KISCR, NEWELB)

         CALL MCFIND ('NAMELB', IDUM, LNAM)
         CALL MCRSRV ('NAMSCR', KNMSC, LNAM)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 140

         CALL MUNELB (NEWELB, IA(KIELBS), NEWEL,
     &      IA(KIDELB), IA(KNELB), IA(KNLNK), IA(KNATR),
     &      IA(KLINK), A(KATRIB), IA(KLINKO), A(KATRO),
     &      IA(KIXEL), IA(KIXELB), IA(KJNELB), IA(KISCR),
     $      C(KNMLB),  C(KNMSC), LLINK, LATRIB)

         CALL MDDEL ('LINKO')
         CALL MDDEL ('ATRIBO')
         CALL MDDEL ('IXELB')
         CALL MDDEL ('JNELB')
         CALL MDDEL ('ISCR')
         CALL MCDEL ('NAMSCR')
         CALL MDLONG('LINK',  KLINK,  LLINK) 
         CALL MDLONG('ATRIB', KATRIB, LATRIB) 
      END IF

      CALL MDDEL ('IELBST')

C   --Mark if any elements are deleted

      DELEL = NEWEL .LT. (NUMEL1 + NUMEL2)

      IF (DELEL) THEN

C      --Make up an index of nodes in the existing element blocks

         CALL MDRSRV ('IXNP', KIXNP, NEWNP)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 140

         N = NEWNP
         CALL ZMFIXD (NEWELB, IA(KNELB), IA(KNLNK), IA(KLINK),
     &        N, IA(KIXNP))

         DELNP = (N .LT. NEWNP)

         IF (.NOT. DELNP) THEN
            CALL MDDEL ('IXNP')
         END IF
      ELSE
         DELNP = .FALSE.
      END IF

C   --Squeeze the coordinates

      IF (DELNP) THEN
         CALL ZMXYZ (NDIM, NEWNP, IA(KIXNP),
     &      A(KXN), A(KYN), A(KZN))

         CALL MDLONG ('XN', KXN, NEWNP)
         CALL MDLONG ('YN', KYN, NEWNP)
         CALL MDLONG ('ZN', KZN, NEWNP)
      END IF

C   --Renumber the element map
c$$$
c$$$      IF (TWODB) THEN
c$$$         CALL RENIX (NUMEL2, NUMEL1, IDUM, IA(KMAPEL+NUMEL1))
c$$$      END IF
c$$$
c$$$      IF (RENEL) THEN
c$$$         CALL RENIX (NUMEL1+NUMEL2, -999, IA(KIXEL), IA(KMAPEL))
c$$$      END IF
c$$$
c$$$C   --Squeeze the element map
c$$$
c$$$      IF (DELEL) THEN
c$$$         NEW = NUMEL1+NUMEL2
c$$$         CALL ZMMAP (NEW, IA(KMAPEL))
c$$$
c$$$         CALL MDLONG ('MAPEL', KMAPEL, NEW)
c$$$      END IF
C ... The above code assumes that the element map is a permutation of
C     the sequence (1..numel). For example, if an optimizer has been 
C     run on the input databases.  It will fail if the map contains
C     values >numel which can happen in some instances.  Since
C     an optimization would have to be redone for the combined mesh and
C     there is no good way to ensure that combining 2 or more arbitrary
C     maps will give unique ids, we just punt and create a map which is
C     1..numel.  Since we don't need it until output, we allocate and
C     create it at the expmap call.
      
C   --Renumber the element block nodes

      IF (DELNP) THEN
         CALL RENELB (NEWELB, -999, IA(KIXNP),
     &      IA(KNELB), IA(KNLNK), IA(KLINK))
      END IF

C   --Renumber the nodal point set nodes

      IF (DELNP) THEN
         CALL RENIX (LNPSN1+LNPSN2, -999, IA(KIXNP), IA(KLTNNS))
      END IF

C   --Renumber the element side set elements

      IF (TWODB) THEN
         CALL RENIX (LESSE2, NUMEL1, IDUM, IA(KLTESS+LESSE1))
      END IF

      IF (RENEL) THEN
         CALL RENIX (LESSE1+LESSE2, -999, IA(KIXEL), IA(KLTESS))
      END IF

      IF (RENEL) THEN
         CALL MDDEL ('IXEL')
      END IF

C   --"Munch" the nodal point sets

      I = INTCNT (0, IA(KINPSS), NEWNPS)

      IF ((I .LT. NEWNPS) .OR. DELNP) THEN
         CALL MDLONG ('LTNNPS', KLTNNS, NEWNNL)
         CALL MDLONG ('FACNPS', KFACNS, NEWNNL)

         CALL MDRSRV ('LTNNPO', KLTNNO, NEWNNL)
         CALL MDRSRV ('FACNPO', KFACNO, NEWNNL)
         CALL MDRSRV ('IXNNPO', KIXNNO, NEWNPS)
         CALL MDRSRV ('NNNPO',  KNNNO,  NEWNPS)
         CALL MDRSRV ('ISCR',   KISCR,  NEWNPS)
         call mdrsrv ('nodscr', kndscr, newnp)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 140

         CALL MUNNPS (NEWNPS, IA(KINPSS), NEWNNL,
     &      IA(KIDNS), IA(KNNNS), IA(KIXNNS), IA(KLTNNS), A(KFACNS),
     &      IA(KLTNNO), A(KFACNO), IA(KIXNNO), IA(KNNNO), IA(KISCR),
     *      IA(KNDSCR), NEWNP)

         CALL MDDEL ('LTNNPO')
         CALL MDDEL ('FACNPO')
         CALL MDDEL ('IXNNPO')
         CALL MDDEL ('NNNPO')
         CALL MDDEL ('ISCR')
         call mddel ('nodscr')

C      --Squeeze the nodal point sets

         IF (DELNP) THEN
            CALL ZMNPS (NEWNPS, NEWNNL,
     &         IA(KIDNS), IA(KNNNS), IA(KIXNNS), IA(KLTNNS), A(KFACNS))
         END IF

         CALL MDLONG ('IDNPS', KIDNS, NEWNPS)
         CALL MDLONG ('NNNPS', KNNNS, NEWNPS)
         CALL MDLONG ('IXNNPS', KIXNNS, NEWNPS)
         CALL MDLONG ('LTNNPS', KLTNNS, NEWNNL)
         CALL MDLONG ('FACNPS', KFACNS, NEWNNL)
      END IF

      CALL MDDEL ('INPSST')

C   --"Munch" the element side sets

      I = INTCNT (0, IA(KIESSS), NEWESS)

      IF ((I .LT. NEWESS) .OR. DELEL) THEN
         CALL MDLONG ('LTEESS', KLTESS, NEWSEL)

         CALL MDRSRV ('LTEESO', KLTESO, NEWSEL)
         CALL MDRSRV ('LTSSO',  KLTSSO, NEWSEL)
         CALL MDRSRV ('FACS0',  KFACS0, NEWSDL)
         CALL MDRSRV ('IXEESO', KIXESO, NEWESS)
         CALL MDRSRV ('IXEDS0', KIXDS0, NEWESS)
         CALL MDRSRV ('NEESO',  KNESO,  NEWESS)
         CALL MDRSRV ('NEDS0',  KNDS0,  NEWESS)
         CALL MDRSRV ('ISCR',   KISCR,  NEWESS)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 140

         CALL MUNESS (NEWESS, IA(KIESSS), NEWSEL, NEWSDL, 
     &     IA(KIDSS), IA(KNESS), IA(KNDSS), IA(KIXESS), IA(KIXDSS),
     &     IA(KLTESS), IA(KLTSSS), A(KFACSS), 
     &     IA(KLTESO), IA(KLTSSO), A(KFACS0), IA(KIXESO), IA(KIXDS0),
     &     IA(KNESO), IA(KNDS0), IA(KISCR))

         CALL MDDEL ('LTEESO')
         CALL MDDEL ('LTSSO')
         CALL MDDEL ('FACS0')
         CALL MDDEL ('IXEESO')
         CALL MDDEL ('IXEDS0')
         CALL MDDEL ('NEESO')
         CALL MDDEL ('NEDS0')
         CALL MDDEL ('ISCR')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 140

C      --Squeeze the element side sets

         IF (DELEL) THEN
           CALL ZMESS (NEWESS, NEWSEL, NEWSDL, 
     &       IA(KIDSS), IA(KNESS), IA(KNDSS), IA(KIXESS),
     *       IA(KIXDSS), IA(KLTESS), IA(KLTSSS), IA(KLTSNC), A(KFACSS))
         END IF

         CALL MDLONG ('IDESS', KIDSS, NEWESS)
         CALL MDLONG ('NEESS', KNESS, NEWESS)
         CALL MDLONG ('IXEESS', KIXESS, NEWESS)
         CALL MDLONG ('LTEESS', KLTESS, NEWSEL)
         CALL MDLONG ('LTSESS', KLTSSS, NEWSEL)
         CALL MDLONG ('FACESS', KFACSS, NEWSDL)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 140
      END IF

      CALL MDDEL ('IESSST')

      IF (DELNP) THEN
         CALL MDDEL ('IXNP')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 140
      END IF

C   --Find out if done processing (MOVED TO COMAND -- EXIT/ADD/FINISH)
      IF (.NOT. DONE) THEN
        numnp1 = newnp
        numel1 = newel
        nelbl1 = newelb
        nnps1  = newnps
        lnpsn1 = newnnl
        ness1  = newess
        lesse1 = newsel
        lessd1 = newsdl

C ... Reset array sizes to match current combined database
        call mdlong('XN', KXN, newnp)
        call mdlong('YN', KYN, newnp)
        call mdlong('ZN', KZN, newnp)
C        call mdlong('MAPEL', KMAPEL, newel)
        call mdlong('IDELB', KIDELB, newelb)
        call mdlong('IDNPS', KIDNS, newnps)
        call mdlong('LTNNPS', KLTNNS, newnnl)
        call mdlong('IDESS', KIDSS, newess)
        call mdlong('LTEESS', KLTESS, newsel)
        CALL MDSTAT (NERR, MEM)
        IF (NERR .GT. 0) GOTO 140
         GOTO 100
       END IF
       
C   --Write out the new database

  130 CONTINUE
      IF (DONE) THEN
         WRITE (*, *)
         CALL GETINP (0, 0, 'Output file> ', FILNAM, IOSTAT)
         CALL OUTLOG (KLOG, 1, 0, FILNAM, IDUM, RDUM)

         IF (IOSTAT .LT. 0) GOTO 150
         IF (FILNAM .EQ. ' ') GOTO 130
C   --Write the QA records

         IF (NQAREC .LT. MAXQA) THEN
            NQAREC = NQAREC + 1
            QAREC(1,NQAREC) = QAINFO(1)
            QAREC(2,NQAREC) = QAINFO(3)
            QAREC(3,NQAREC) = QAINFO(5)
            QAREC(4,NQAREC) = QAINFO(6)
         END IF
      ELSE
         FILNAM = '%gjoin'
      END IF


      CALL WRGEN (A, A, FILNAM, TITLE, NDIM, NEWNP, NEWEL, NEWELB,
     &   NEWNPS, NEWNNL, NEWESS, NEWSEL, NEWSDL, 
     &   KXN, KYN, KZN, KMAPEL,
     &   KIDELB, KNELB, KNLNK, KNATR, KLINK, KATRIB,
     &   KIDNS, KNNNS, KIXNNS, KLTNNS, KFACNS,
     &   KIDSS, KNESS, KNDSS, KIXESS, KIXDSS, KLTESS, KFACSS,
     &   kltsss, NQAREC, QAREC, NINFO, INFREC, C(KNMLB), *140)

      FIRST = .FALSE.

      IF (.NOT. DONE) GOTO 80
      GOTO 150

  140 CONTINUE
      CALL MEMERR
      GOTO 150

  150 CONTINUE
      CALL WRAPUP (QAINFO(1))
      call addlog (QAINFO(1))
      OPEN (UNIT=9, FILE='%gjoin', FORM='unformatted',
     &   STATUS='old', IOSTAT=IERR)
      IF (IERR .EQ. 0) THEN
         CLOSE (9, STATUS='DELETE')
      END IF
 9000 format (A1,'new = ',1pe10.3,' * ',A1,'old + ',1pe10.3)
      END

C...Check whether model contains elements of a single topology.
C   This is currently used in the sideset mirroring code
      subroutine chktop(nelblk, namelb, comtop)

      include 'exodusII.inc'
      integer nelblk
      character*(MXSTLN) namelb(nelblk)
      character*(MXSTLN) comtop
      
      comtop = namelb(1)
      do 10 i=2, nelblk
         if (namelb(i) .ne. comtop) then
            comtop = 'MULTIPLE_TOPOLOGIES'
            return
         end if
 10   continue
      return
      end
