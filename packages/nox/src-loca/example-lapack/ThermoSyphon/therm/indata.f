C----------------------------------------------------------------------|
      subroutine read_indata(lu_in,m_max,n_max,m,prandtl,rayleigh,bval1,
     &    bval2,file1,ra_step,ra_final,pr_final,jobvL,jobvR,tol)
C----------------------------------------------------------------------|
      implicit none
        integer lu_in, m_max, n_max
        include 'indata.h'
        character*8 file1
        
C----------------------------------------------------------------------|
cBegin executables
	open(unit=lu_in,file=file1,status='old',form='formatted')
        read(lu_in,*) debug
        read(lu_in,*) m
        read(lu_in,*) n
        read(lu_in,*) prandtl
        read(lu_in,*) sigma
        read(lu_in,*) rayleigh
        read(lu_in,*) bval1 
        read(lu_in,*) bval2
        read(lu_in,*) bc_case
       read(lu_in,*) oldfile
       read(lu_in,*) tol 
       read(lu_in,*) ra_step 
       read(lu_in,*) ra_final 
       read(lu_in,*) pr_final 
       read(lu_in,*) jobvL
       read(lu_in,*) jobvR
        close(lu_in)
c  create name for savefile
      write(unit=savefile,fmt=41) prandtl,rayleigh
41    format('./save_',f9.7,'_',f9.7)

C----------------------------------------------------------------------|
c... check input data
        if (m .gt. m_max .or. m .lt. 0) then
          stop 'm is out of range'
        endif
        pi = 4.0d0*datan(1.0d0)

      end
