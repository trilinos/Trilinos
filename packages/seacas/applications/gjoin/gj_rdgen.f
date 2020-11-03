C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RDGEN (A, IA, C, FIRST, FILNAM,
     &   TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSDL,
     &   KXN, KYN, KZN, KMAPEL,
     &   KIDELB, KNELB, KNLNK, KNATR, KLINK, KATRIB,
     &   KIDNS, KNNNS, KIXNNS, KLTNNS, KFACNS,
     &   KIDSS, KNESS, KNDSS, KIXESS, KIXDSS, KLTESS, kltsss,
     &   kltsnc, kfacss, NQAREC, QAREC, NINFO, INFREC, KNMLB, USESDF, *)
C=======================================================================

C   --*** RDGEN *** (GJOIN) Read the GENESIS database
C   --   Written by Amy Gilkey - revised 12/04/87
C   --
C   --RDGEN expands the memory for the GENESIS database then reads the
C   --database.  Note that the memory is reserved (with a length of zero)
C   --in INIGEN.  File may be deleted.
C   --
C   --Parameters:
C   --   A - IN/OUT - the dynamic memory base array
C   --   FIRST - IN - true iff file is NOT to be deleted
C   --   FILNAM - IN - the database filename
C   --   TITLE - OUT - the database title
C   --   NDIM - OUT - the number of coordinates per node
C   --   NUMNP - OUT - the number of nodes
C   --   NUMEL - OUT - the number of elements
C   --   NELBLK - OUT - the number of element blocks
C   --   NUMNPS - OUT - the number of nodal point sets
C   --   LNPSNL - OUT - the length of the nodal point sets node list
C   --   NUMESS - OUT - the number of side sets
C   --   LESSEL - OUT - the length of the element side sets element list
C   --   KXN, KYN, KZN - OUT - index of XN, YN, ZN; nodal coordinates
C   --   KMAPEL - OUT - index of MAPEL; the element order map
C   --   KIDELB - OUT - index of IDELB; the element block IDs for each block
C   --   KNELB - OUT - index of NUMELB; the number of elements in each block
C   --   KNLNK - OUT - index of NUMLNK; the number of nodes per element
C   --      in each block
C   --   KNATR - OUT - index of NUMATR; the number of attributes in each block
C   --   KLINK - OUT - index of LINK; the connectivity for each block
C   --   KATRIB - OUT - index of ATRIB; the attributes for each block
C   --   KIDNS - OUT - index of IDNPS; the nodal point set ID for each set
C   --   KNNNS - OUT - index of NNNPS; the number of nodes for each set
C   --   KIXNNS - OUT - index of IXNNPS; the index of the first node
C   --      for each set
C   --   KLTNNS - OUT - index of LTNNPS; the nodes for all sets
C   --   KFACNS - OUT - index of FACNPS; the distribution factors for all sets
C ---------sidesets----------
C   --   KIDSS - OUT - index of IDESS; the element side set ID for each set
C   --   KNESS - OUT - index of NEESS; the number of elements for each set
C   --   KNDSS - OUT - index of NDESS; the number of dist-fact for each set
C   --   KIXESS - OUT - index of IXEESS; the index of the first element
C   --      for each set
C   --   KIXDSS - OUT - index of IXDESS; the index of the first dist-fact for each set
C   --   KLTESS - OUT - index of LTEESS; the elements for all sets
C   --   KFACSS - OUT - index of FACESS; the distribution factors for all sets
C   --   kltsss - OUT - index of LTSESS; the sides for all sets
C   --   kltsnc - OUT - index of LTSSNC; the df count for each element/face in the list
C   --   NAMELB - OUT - names of the element blocks

      include 'exodusII.inc'
      include 'gj_params.blk'

      DIMENSION A(*), IA(*)
      CHARACTER*1 C(*)

      LOGICAL FIRST, TEMP, USESDF

      CHARACTER*(*) FILNAM
      character*(MXLNLN) title, tmpstr
      character*(MXSTLN) qarec(4,MAXQA)
      character*(MXLNLN) infrec(MAXINF)
      character*(MXSTLN) blname
      character*132 errmsg, name
C      --QAREC - the QA records
C      --INFREC - the information records
      CHARACTER*(MXSTLN) CDUMMY
      CHARACTER*80 TMPFIL

      character*20 stra, strb, strc

      DATA TMPFIL /'%gjoin'/
      DATA NPART /0/

C   --Initialize

      TITLE = ' '
      NDIM = 0
      NUMNP = 0
      NUMEL = 0
      NELBLK = 0
      NUMNPS = 0
      LNPSNL = 0
      NUMESS = 0
      LESSEL = 0

      CALL SQZSTR(FILNAM, LNAM)
      LNAM = LENSTR (FILNAM)
      IF (FILNAM .EQ. TMPFIL) THEN
         TEMP = .TRUE.
      ELSE
         TEMP = .FALSE.
      ENDIF

C     Make netCDF and exodus errors not show up

      call exopts (0, ierr)

C   --Open the netcdf file
      icpuws = 0
      iows   = 0
      netid = exopen(filnam(:lnam), EXREAD, icpuws, iows, vers, nerr)
      if (nerr .lt. 0) then
         write(errmsg,10) filnam(:lnam), nerr
 10      format("Could not open exodusII file '",A,"', error = ",i3)
         call exerr ('gjoin2', errmsg, exlmsg)
         goto 960
      endif

      call exinq (netid, EXVERS, idummy, versi, cdummy, nerr)
      write(*,'(A,F6.3)')
     & 'This database was created by ExodusII version ', versi

C   --Read global information from the database

      call exgini (netid, title, ndim, numnp, numel,
     &             nelblk, numnps, numess, nerr)
      if (nerr .lt. 0) then
         call exerr ('gjoin2', 'Error from exgini', exlmsg)
         goto 960
      endif

C     Get the length of the node sets node list

      if (numnps .gt. 0) then
         call exinq (netid, EXNSNL, lnpsnl, dummy, cdummy, nerr)
         if (nerr .lt. 0) then
            call exerr ('gjoin2', 'Error from exqini', exlmsg)
            goto 960
         endif
      else
         lnpsnl = 0
      endif

      if (numess .gt. 0) then

C        Get the length of the side sets distribution factor list

         call exinq (netid, EXSSDF, lessdl, dummy, cdummy, nerr)
         if (nerr .lt. 0) then
            call exerr ('gjoin2', 'Error from exqini', exlmsg)
            goto 960
         endif

C        Get the length of the side sets element list

         call exinq (netid, EXSSEL, lessel, dummy, cdummy, nerr)
         if (nerr .lt. 0) then
            call exerr ('gjoin2', 'Error from exqini', exlmsg)
            goto 960
         endif
      else
         lessel = 0
         lessdl = 0
      endif

      IF (.NOT. TEMP) THEN
         CALL DBPINI ('NTIS', NDB, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &                NUMNPS, LNPSNL, lnpsnl, NUMESS, LESSEL,
     &                lessdl, IDUM, IDUM, IDUM, FILNAM(:LNAM))
      ENDIF

C   --Read the coordinates

      CALL MDFIND ('XN', KXN, LOLD)
      CALL MDLONG ('XN', KXN, LOLD+NUMNP)
      CALL MDLONG ('YN', KYN, LOLD+NUMNP)
      IF (NDIM .GE. 3) CALL MDLONG ('ZN', KZN, LOLD+NUMNP)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 950

      if (ndim .ge. 3) then
         call exgcor(netid, a(kxn+lold), a(kyn+lold), a(kzn+lold),
     &               nerr)
      else
         call exgcor(netid, a(kxn+lold), a(kyn+lold), dummy, nerr)
      endif
      if (nerr .lt. 0) then
         call exerr ('gjoin2', 'Error from exgcor', exlmsg)
         goto 960
      endif

C   --Read the element order map
C ... See note in gjoin.f about element map
c$$$      CALL MDFIND ('MAPEL', KMAPEL, LOLD)
c$$$      CALL MDLONG ('MAPEL', KMAPEL, LOLD+NUMEL)
c$$$      CALL MDSTAT (NERR, MEM)
c$$$      IF (NERR .GT. 0) GOTO 950
c$$$
c$$$      call exgmap (netid, a(kmapel+lold), nerr)
c$$$      if (nerr .ne. 0) then
c$$$         if (nerr .eq. 17) then
c$$$
c$$$C           -- no element order map in the EXODUS II file
c$$$C           -- create a dummy one
c$$$            do 30 i=1,numel
c$$$               ia(kmapel+lold+i-1) = i
c$$$ 30         continue
c$$$         else
c$$$            goto 950
c$$$         endif
c$$$      endif

C   --Read in the element block ID array

      CALL MDFIND ('IDELB', KIDELB, LOLD)
      LOLDBL = LOLD
      CALL MDLONG ('IDELB', KIDELB, LOLD+NELBLK)
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 950

      call exgebi (netid, a(kidelb+lold), nerr)
      if (nerr .lt. 0) then
         call exerr ('gjoin2', 'Error from exgebi', exlmsg)
         goto 960
      endif

C   --Read the element blocks

      CALL MDLONG ('NUMELB', KNELB, LOLD+NELBLK)
      CALL MDLONG ('NUMLNK', KNLNK, LOLD+NELBLK)
      CALL MDLONG ('NUMATR', KNATR, LOLD+NELBLK)
      CALL MCLONG ('NAMELB', KNMLB, (LOLD+NELBLK)*MXSTLN)
      call mdfind ('LINK', klink, loldlk)
      call mdfind ('ATRIB', katrib, loldat)

      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 950
      loldlks = loldlk
      loldats = loldat

      do 40 ielb = 1, nelblk
         ioff = lold+ielb-1
         idelb = ia(kidelb+ioff)
         call exgelb (netid, idelb, blname, num, numlnk,
     &                numatr, nerr)
         if (nerr .lt. 0) then
            call exerr ('gjoin2', 'Error from exgelb', exlmsg)
            goto 960
         endif

C ... Wrapper to deal with character*1 'C' array
         call cpnam(blname, c(knmlb + mxstln*ioff))
         ia(knelb+ioff) = num
         ia(knlnk+ioff) = numlnk
         ia(knatr+ioff) = numatr

         lnewlk = loldlk+(num*numlnk)
         loldlk = lnewlk

         lnewat = loldat+(num*numatr)
         loldat = lnewat
 40   continue

      call mdlong ('LINK', klink, lnewlk)
      call mdlong ('ATRIB', katrib, lnewat)
      call mdstat (nerr, mem)
      if (nerr .gt. 0) goto 950

      loldlk = loldlks
      loldat = loldats
      nel = 0
      do 50 ielb = 1, nelblk
         ioff   = lold+ielb-1
         idelb  = ia(kidelb+ioff)
         numlnk = ia(knlnk+ioff)
         numatr = ia(knatr+ioff)
         num    = ia(knelb+ioff)

         if (numlnk .gt. 0) then
            lnewlk = loldlk+(num*numlnk)
            call exgelc (netid, idelb, a(klink+loldlk), nerr)
            if (nerr .lt. 0) then
               call exerr ('gjoin2', 'Error from exgelc', exlmsg)
               goto 960
            endif
            loldlk = lnewlk
         endif

         if (numatr .gt. 0) then
            lnewat = loldat+(num*numatr)
            call exgeat (netid, idelb, a(katrib+loldat), nerr)
            if (nerr .lt. 0) then
               call exerr ('gjoin2', 'Error from exgeat', exlmsg)
               goto 960
            endif
            loldat = lnewat
         endif
 50   continue

C   --Read the node sets

      CALL MDFIND ('IDNPS', KIDNS, LOLD)
      CALL MDLONG ('IDNPS', KIDNS, LOLD+NUMNPS)
      CALL MDLONG ('NNNPS', KNNNS, LOLD+NUMNPS)
      call mdlong ('NDNPS', kndns, numnps)   ! Node set df count array
      CALL MDLONG ('IXNNPS', KIXNNS, LOLD+NUMNPS)
      call mdlong ('IXDNPS', kixdns, numnps) ! Node set df index array
      CALL MDFIND ('LTNNPS', KLTNNS, LOLD2)
      CALL MDLONG ('LTNNPS', KLTNNS, LOLD2+LNPSNL)
      CALL MDLONG ('FACNPS', KFACNS, LOLD2+LNPSNL)
      call mdlong ('CFACNP', kcfacn, lnpsnl) ! Compressed df list array
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 950

      if (numnps .gt. 0) then
         call exgcns (netid, a(kidns+lold), a(knnns+lold),
     &                a(kndns), a(kixnns+lold), a(kixdns),
     &                a(kltnns+lold2), a(kcfacn), nerr)
         if (nerr .lt. 0) then
            call exerr ('gjoin2', 'Error from exgcns', exlmsg)
            goto 960
         endif
      endif

C     Massage node sets distribution factors to include 'gj_1' for node sets
C     without DFs by walking KNDNS array, checking for 0, and filling where
C     necessary.

      do 80 i=0, numnps-1
         if (ia(kndns+i) .eq. 0) then
            do 60 ii=0, ia(knnns+lold+i)-1
               a(kfacns+lold2+ia(kixnns+lold+i)-1+ii) = 1.0
 60         continue
         else
            do 70 ii=0, ia(kndns+i)-1
               a(kfacns+lold2+ia(kixnns+lold+i)-1+ii) =
     &              a(kcfacn+ia(kixdns+i)-1+ii)
 70         continue
         endif
 80   continue

C   --Read the side sets

      CALL MDFIND ('IDESS', KIDSS, LOLD)
      CALL MDLONG ('IDESS', KIDSS, LOLD+NUMESS)
      CALL MDLONG ('NEESS', KNESS, LOLD+NUMESS)
      call mdlong ('NDESS', kndss, LOLD+numess)      ! number of dist factors array
      CALL MDLONG ('IXEESS', KIXESS, LOLD+NUMESS)
      call mdlong ('IXDESS', kixdss, LOLD+numess)    ! index into dist factors array
      CALL MDFIND ('LTEESS', KLTESS, LOLD2)
      CALL MDLONG ('LTEESS', KLTESS, LOLD2+LESSEL)
      call mdlong ('LTSESS', kltsss, lold2+lessel)    ! side list
      call mdlong ('LTSSNC', kltsnc, lold2+lessel)
      call mdfind ('FACESS', KFACSS, LOLD3)
      call mdlong ('FACESS', kfacss, lold3+lessdl)    ! Compressed dist factors list
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 950

      if (numess .gt. 0) then
         call exgcss (netid, a(kidss+lold), a(kness+lold),
     &    a(kndss+lold), a(kixess+lold), a(kixdss+lold),
     &    a(kltess+lold2), a(kltsss+lold2), a(kfacss+lold3), nerr)
         if (nerr .lt. 0) then
            call exerr ('gjoin2', 'Error from exgcss', exlmsg)
            goto 960
         endif

         call exgcssc(netid, ia(kltsnc+lold2), nerr)
         if (nerr .lt. 0) then
           call exerr ('gjoin2', 'Error from exgcssc', exlmsg)
           goto 960
         endif

         ioff = 0
         ntot = 0
         do 90 iess = 1, numess
           id = ia(kidss+lold+iess-1)
           call chksnc(ia(kltsnc+lold2+ioff), ia(kness+lold+iess-1),
     *       ncnt)
           if (ncnt .ne. ia(kndss+lold+iess-1) .and.
     *       ia(kndss+lold+iess-1) .gt. 0 ) then

              CALL INTSTR (1, 0, ncnt, STRA, LSTRA)
              CALL INTSTR (1, 0, id,   STRC, LSTRC)
              CALL INTSTR (1, 0, ia(kndss+lold+iess-1), STRB, LSTRB)
              CALL PRTERR ('ERROR',
     &             'For sideset ' // STRC(:LSTRC)
     &             // ' the distribution factor count of '
     $             // STRB(:LSTRB)
     $             // ' does not match the face node count of '
     $             // STRA(:LSTRA))
             stop 'internal dist-factor count error'
           end if
           ntot = ntot + ncnt
           ioff = ioff + ia(kness+lold+iess-1)
 90      continue
         if (ntot .ne. lessdl .and. lessdl .gt. 0) then
           stop 'internal nodecount error'
         end if
      endif

C   --Read the QA and information records

      IF (.NOT. TEMP) THEN

         kqarec = 0
         call exinq (netid, EXQA, kqarec, r, name, nerr)
         if (nerr .lt. 0) then
            call exerr ('gjoin2', 'Error from exinq', exlmsg)
            goto 960
         endif

         if (kqarec .gt. 0) then
            if ((nqarec + kqarec) .le. MAXQA) then
               call exgqa (netid, qarec(1, nqarec+1), nerr)
               if (nerr .lt. 0) then
                  call exerr ('gjoin2', 'Error from exgqa', exlmsg)
                  goto 960
               endif
               nqarec = nqarec + kqarec
            endif
         endif

         kinfo = 0
         call exinq (netid, EXINFO, kinfo, r, name, nerr)
         if (nerr .lt. 0) then
            call exerr ('gjoin2', 'Error from exinq', exlmsg)
            goto 960
         endif

         if (kinfo .gt. 0) then
            if ((ninfo + kinfo) .le. MAXINF) then
               call exginf (netid, infrec(ninfo+1), nerr)
               if (nerr .lt. 0) then
                  call exerr ('gjoin2', 'Error from exginf', exlmsg)
                  goto 960
               endif

               NPART  = NPART + 1
               DO 320 INFO = NINFO+1, (NINFO+KINFO)
                  TMPSTR = INFREC(INFO)
                  WRITE (INFREC(INFO), 310) NPART, TMPSTR(:77)
 310              FORMAT (I2.2,'/',A)
 320           CONTINUE
               ninfo = ninfo + kinfo
            endif
         endif

C        Add info block

         IF (NINFO .LT. MAXINF) THEN
            NINFO = NINFO + 1
            WRITE (INFREC(NINFO), 310) NPART, FILNAM(:77)
         ENDIF
      ENDIF

      IF (FIRST) THEN
         call exclos (netid, ierr)
      ELSE
         IF (FILNAM .NE. TMPFIL) THEN
            CALL PRTERR ('PROGRAM', 'in RDGEN')
         ENDIF
         call exclos (netid, ierr)
      ENDIF

      RETURN

 950  CONTINUE
      call memerr(6)

 960  CONTINUE
      RETURN 1
      END

      subroutine cpnam(namin, namout)
      character*32 namin, namout
      namout = namin
      return
      end

      subroutine chksnc(lsnc, len, ncnt)
C ... Sum the face counts in 'lsnc' so can check that the sum equals the
C     df count for this list
      integer lsnc(len)
      ncnt = 0
      do 10 i=1, len
        ncnt = ncnt + lsnc(i)
 10   continue
      return
      end

      subroutine fixdf(len, df)
      real df(*)
      return
      end
