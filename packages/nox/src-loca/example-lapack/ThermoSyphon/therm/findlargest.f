c subroutine to find the largest entry in an array
      subroutine findlargest(array, size, index)
      implicit none
      integer size, index
      double precision array(size)
c  local variables
      integer tempindex
      integer i
      double precision tempmax
c begin executables
      tempmax=array(1)
      index=1
      do i=2, size
        if (array(i) .gt. tempmax) then
	  tempmax=array(i)
	  index=i
	endif
      enddo
      return
      end
