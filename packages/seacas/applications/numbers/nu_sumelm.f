C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SUMELM (CRD, DISP, SVAR, MAT, NDIM, NUMNP, INDX,
     *     ELMSEL, NAME, TIME, ITMSEL, AVER, AXI, DOABS, DODENS,
     &     NUMEL, link, nnodes, nelblk, volume, ISEVOK, SUM,
     &     nsel)
C=======================================================================

      include 'exodusII.inc'
      include 'nu_varcnt.blk'

      REAL CRD(NUMNP,*), DISP(NUMNP,*), SVAR(*)
      INTEGER MAT(6,*)
      LOGICAL ISEVOK(nvarel,*)
      CHARACTER*(*) NAME
      REAL TIME(*)
      LOGICAL ELMSEL(*), ITMSEL(*), AXI, AVER, DOABS, DODENS
      integer link(nnodes, *)
      real volume(*)
      REAL SUM(*)
      integer nsel(*)

      CHARACTER*32 STRA
      CHARACTER*132 STRB
      include 'nu_io.blk'
      include 'nu_ptim.blk'
      include 'nu_ndisp.blk'

      NLAST = 0

      IF (AVER) THEN
         IF (DOABS) THEN
            STRA = 'Average absolute value of'
         ELSE
            STRA = 'Average'
         END IF
         RSEL = DBLE(NUMEQL (.TRUE., NUMEL, SELECT))
      ELSE
         IF (DOABS) THEN
            STRA = 'Absolute value total of'
         ELSE
            STRA = 'Total'
         END IF
         RSEL = 1.0
      END IF
      if (dodens) then
         strb = '(Multiplied by element volumes)'
      else
         strb = ' '
      end if
      DO 10 IO=IOMIN, IOMAX
         IF (NDIM .EQ. 2 .AND. AXI) THEN
            WRITE (IO,20) STRA(:LENSTR(STRA)), NAME(:LENSTR(NAME)),
     &       strb(:lenstr(strb))
         ELSE
            WRITE (IO,30) STRA(:LENSTR(STRA)), NAME(:LENSTR(NAME)),
     &       strb(:lenstr(strb))
         END IF
         WRITE (IO, 40) (MAT(1,iblk),iblk=1,nelblk)
 10   CONTINUE
 20   FORMAT (/4X,A,' ',A,' on selected elements ',
     *     '(2 PI Radius Multiplier) ',A)
 30   FORMAT (/4X,A,' ',A,' on selected elements ',A)
 40   FORMAT (/,4X,'Time        Total  ',32(I8,4x))

C ... Determine number of selected elements in each block
      DO 45 iblk = 1, nelblk
         if (mat(5, iblk) .eq. 1) then
            ielbeg = mat(3, iblk)
            ielend = mat(4, iblk)
            nsel(iblk) = 0
            do 44 iel = ielbeg, ielend
               if (elmsel(iel)) nsel(iblk) = nsel(iblk) + 1
 44         continue
         else
            nsel(iblk) = 0
         end if
 45   continue

 50   CONTINUE
      NLAST = NLAST + 1
      IF (NLAST .GT. LSTSEL) THEN
         NLAST = 0
         GO TO 120
      ELSE IF (ITMSEL(NLAST)) THEN

C     ... If by density, then calculate volume of elements
         IF (DODENS) THEN
C     ... READ THE STEP AND STORE DISPLACEMENTS
           do 54 i=1, ndim
             call exgnv(ndb, nlast, ndisp(i), numnp, disp(1,i), ierr)
 54        continue
           call elvol(ndim, crd, disp, link, numnp, numel, nnodes,
     &       volume)
         else
           call inirea (numel, 1.0, volume)
         end if

C     ... Zero the volume of all non-selected elements
         do 55 iel = 1, numel
            if (.not. elmsel(iel)) volume(iel) = 0.0
 55      continue

C     ... Read and store the element variable
         call inirea(numel, 0.0, svar)
         ioff = 1
         do 56 i=1, nelblk
           if (isevok(indx, i)) then
             call exgev(ndb, nlast, indx, mat(1,i), mat(2,i),
     &         svar(ioff), ierr)
           end if
 56      continue
         TREAD = TIME(NLAST)

         DO 60 I = 1, NELBLK+1
            SUM(I) = 0.0
 60      CONTINUE

         neb1 = nelblk + 1
         do 90 iblk = 1, nelblk
            if (mat(5, iblk) .ne. 1) goto 80
            ielbeg = mat(3, iblk)
            ielend = mat(4, iblk)
C ... Note that volume is set to zero for all non-selected elements
            if (doabs) then
               do 70 iel = ielbeg, ielend
                  SUM(iblk) = SUM(iblk) + ABS(SVAR(I) * volume(i))
 70            continue
            else
               do 75 iel = ielbeg, ielend
                  SUM(IBLK) = SUM(IBLK) + SVAR(I) * volume(i)
 75            continue
            end if
            sum(neb1) = sum(neb1) + sum(iblk)
            if (aver) then
               sum(iblk) = sum(iblk) / nsel(iblk)
            end if
 80         continue
 90      continue
         if (aver) then
            sum(neb1) = sum(neb1) / rsel
         end if

         DO 100 IO=IOMIN,IOMAX
            WRITE (IO, 110) TREAD, sum(nelblk+1), (SUM(i),i=1, nelblk)
 100     CONTINUE
 110     FORMAT (1X,32(1PE15.8,2X))

      END IF
      GO TO 50

 120  CONTINUE
      RETURN

      CALL PRTERR ('PROGRAM',
     *     'Internal code error, contact sponsor')
      STOP 'SUMIT'
      END
