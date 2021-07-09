C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C-----------------------------------------------------------------------
      subroutine prnshd (nelblk, idelb, ishdcl, shdcol, ielbst)
C-----------------------------------------------------------------------
      INTEGER IDELB(*)
      INTEGER ISHDCL(3,NELBLK)
      REAL SHDCOL(7,NELBLK)
      INTEGER IELBST(NELBLK)

      CHARACTER*132 ELBLIN

      do 110 iblk = 1, nelblk
        if (ishdcl(1, iblk) .eq. iblk) then
          if (ishdcl(2, iblk) .gt. 1) then
            elblin = ' Shades'
          else
            elblin = ' Shade'
          end if
          write (*, 920) (shdcol(i,iblk),i=1,6), ishdcl(2, iblk),
     *      elblin(:lenstr(elblin))
          ELBLIN = ' Block IDs: '
          N = 2
          do 100 i = iblk, nelblk
            if (ishdcl(1, i) .eq. ishdcl(1, iblk)) then
              n = n + 1
              WRITE (ELBLIN((N-1)*6+1:N*6), '(I6)', IOSTAT=IDUM)
     *          idelb(i)
              IF (N .GE. 12) THEN
                WRITE (*, 940) ELBLIN(1:LENSTR(ELBLIN))
                N = 0
                ELBLIN = ' '
              END IF
            end if
 100      continue
          IF (N .GT. 0) THEN
            WRITE (*, 940) ELBLIN(1:LENSTR(ELBLIN))
            write (*,*) ' '
          END IF
        end if
 110  continue

 920  FORMAT (' Block Shading: RGB =', 6(f6.3),', ',i3,A)
 940  FORMAT (A)

      RETURN
      end

