program main
  use ,intrinsic :: iso_c_binding ,only : c_int,c_double
  implicit none

!
! C procedure bindings
!

interface 
  integer(c_int) function FEpetra_Map_Create( numGlobalElements ) bind(c,name='FEpetra_Map_Create')
    import :: c_int
    integer(c_int) ,value :: numGlobalElements
  end function FEpetra_Map_Create

  subroutine FEpetra_Map_Destroy( mapID ) bind(c,name='FEpetra_Map_Destroy')
    import :: c_int
    integer(c_int) ,value :: mapID
  end subroutine FEpetra_Map_Destroy

  integer(c_int) function FEpetra_Map_NumGlobalElements( mapID ) bind(c,name='FEpetra_Map_NumGlobalElements')
    import :: c_int
    integer(c_int) ,value :: mapID
  end function FEpetra_Map_NumGlobalElements
end interface

interface 
  integer(c_int) function FEpetra_Vector_Create(  mapID )bind(c,name='FEpetra_Vector_Create')
    import :: c_int
    integer(c_int) ,value :: mapID
  end function FEpetra_Vector_Create

  subroutine FEpetra_Vector_Destroy( vectorID )bind(c,name='FEpetra_Vector_Destroy')
    import :: c_int
    integer(c_int) ,value :: vectorID
  end subroutine FEpetra_Vector_Destroy

  subroutine FEpetra_Vector_Random( vectorID )bind(c,name='FEpetra_Vector_Random')
    import :: c_int
    integer(c_int) ,value :: vectorID 
  end subroutine FEpetra_Vector_Random

  subroutine FEpetra_Vector_Update(vectorID, alpha, vector2ID, beta )bind(c,name='FEpetra_Vector_Update')
    import :: c_int,c_double
    integer(c_int) ,value :: vectorID ,vector2ID
    real(c_double) ,value ::  alpha ,beta
  end subroutine FEpetra_Vector_Update

  real(c_double) function FEpetra_Vector_Norm2( vectorID )bind(c,name='FEpetra_Vector_Norm2')
    import :: c_int,c_double
    integer(c_int) ,value :: vectorID 
  end function FEpetra_Vector_Norm2
end interface

 !
 ! Data declarations 
 !

  integer(c_int) :: numGlobalElements, numGlobalElements_rtn
  integer(c_int) :: mapID    ! TYPE(MapID)
  integer(c_int) :: xID, bID ! TYPE(VectorID)
  real(c_double) :: bnorm, xnorm

! /*
!  * Executable code
!  */
  
! /* Create a map */
  numGlobalElements = 4;
  mapID = FEpetra_Map_Create(numGlobalElements);

  numGlobalElements_rtn = FEpetra_Map_NumGlobalElements(mapID)
  print *,'NumGlobalElements = ', numGlobalElements_rtn
! assert( numGlobalElements == numGlobalElements_rtn )
  
! /* Create vectors */
  xID = FEpetra_Vector_Create(mapID)
  bID = FEpetra_Vector_Create(mapID)

! /* Do some vector operations */
  call FEpetra_Vector_Random(bID)
  call FEpetra_Vector_Update(xID,2.0_c_double,bID,0.0_c_double) ! /* x = 2*b */

  bnorm = FEpetra_Vector_Norm2(bID)
  xnorm = FEpetra_Vector_Norm2(xID)

  print *, "2 norm of x = ", xnorm 
  print *, "2 norm of b = ", bnorm 

! /* Clean up memory (in reverse order)! */
  call FEpetra_Vector_Destroy(bID)
  call FEpetra_Vector_Destroy(xID)
  call FEpetra_Map_Destroy(mapID)

! /* This should throw an exception and print an error message! */
! /* FEpetra_Map_NumGlobalElements(mapID); */

end program main
