!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Zoltan Library for Parallel Applications                                   !
! For more info, see the README file in the top-level Zoltan directory.      ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  CVS File Information :
!     $RCSfile$
!     $Author$
!     $Date$
!     $Revision$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This is a simple preprocessor used for conditional compilation in f90gl.
! This is NOT a general purpose preprocessor.  Do not try to use it for any
! other purpose -- it'll bite you!
! If you need a general purpose preprocessor for Fortran 90, get Michel
! Olagnon's f90ppr from
! http://www.ifremer.fr/ditigo/molagnon/fortran90/contenu.html

module sppr_mod
implicit none

public :: get_directive, get_operands, substitutions, is_equal, &
          is_defined,    add_macro,    init_macro_list

! defined constants

integer, parameter, public :: MAXCHAR = 32 ! maximum length for character strings
integer, parameter, public :: UNDEFINED_DIR = 0, & ! preprocessor directives
                      IF_DIR        = 1, &
                      IFDEF_DIR     = 2, &
                      IFNDEF_DIR    = 3, &
                      ELSE_DIR      = 4, &
                      ENDIF_DIR     = 5, &
                      INCLUDE_DIR   = 6, &
                      DEFINE_DIR    = 7
integer, parameter, public :: NO_OP      = 0, & ! operators in preprocessor statements
                      EQUALS     = 1, &
                      NOT_EQUALS = 2, &
                      ASSIGNMNT  = 3

! linked list type for defined macros

type, private :: macro
   character(len=MAXCHAR) :: name, value
   type (macro), pointer :: next
end type macro

! the list of macros

type(macro), pointer, private :: macro_list

contains

!        -------------
function get_directive(line) result(getdir)
!        -------------
character(len=*), intent(in) :: line
integer :: getdir

character(len=MAXCHAR) :: word
integer :: endofword

! look at the beginning of a line beginning with pound to see
! if it is one of the recognized keywords, and return the keyword

endofword = index(line," ")-1
if (endofword < 0) then
   word = line(2:)
else
   word = line(2:endofword)
end if

if (word == "if") then
   getdir = IF_DIR
else if (word == "ifdef") then
   getdir = IFDEF_DIR
else if (word == "ifndef") then
   getdir = IFNDEF_DIR
else if (word == "else") then
   getdir = ELSE_DIR
else if (word == "endif") then
   getdir = ENDIF_DIR
else if (word == "include") then
   getdir = INCLUDE_DIR
else if (word == "define") then
   getdir = DEFINE_DIR
else
   getdir = UNDEFINED_DIR
end if

return
end function get_directive

!          ------------
subroutine get_operands(line,leftop,op,rightop)
!          ------------
character(len=*), intent(in) :: line
character(len=*), intent(out) :: leftop,rightop
integer, intent(out) :: op

integer :: equal_loc, blank2loc, startop, endop

! returns the values of left_operand operator right_operand

! determine the operator as == or /= or assignment (no = and more than one
! word) or none (no = and one word)

equal_loc = index(line,"=")
if (equal_loc==0) then
   blank2loc = index(line(index(line," ")+1:)," ")
   if (blank2loc == 0) then
      op = NO_OP
   else
      op = ASSIGNMNT
      startop = blank2loc + index(line," ")
      endop = startop
   end if
else if (line(equal_loc-1:equal_loc-1) == "/") then
   op = NOT_EQUALS
   startop = equal_loc-1
   endop = equal_loc
else if (line(equal_loc+1:equal_loc+1) == "=") then
   op = EQUALS
   startop = equal_loc
   endop = equal_loc+1
else ! the = must be part of one of the operands
   blank2loc = index(line(index(line," ")+1:)," ")
   if (blank2loc == 0) then
      op = NO_OP
   else
      op = ASSIGNMNT
      startop = blank2loc + index(line," ")
      endop = startop
   end if
end if

! determine the operands

if (op == NO_OP) then
   leftop = trim(line(index(line," ")+1:))
   rightop = ""
else
   leftop = line(index(line," ")+1:startop-1)
   rightop = trim(line(endop+1:))
end if

return
end subroutine get_operands

!          -------------
subroutine substitutions(line)
!          -------------
character(len=*), intent(in out) :: line

type(macro), pointer :: current_macro
integer :: sub_loc
character(len=132) :: templine ! for PGI beta bug workaround

! looks for macro substitutions in line

current_macro => macro_list
do
   if (.not.associated(current_macro)) then
      exit
   end if
   do
      sub_loc = index(line,trim(current_macro%name))
      if (sub_loc == 0) then
         exit
      end if
! PGI beta bug workaround
! Would like to use
!      line = line(1:sub_loc-1)//trim(current_macro%value)// &
!             line(sub_loc+len_trim(current_macro%name):)
! but a bug in the beta version of the PGI compiler prevents this.  Use temp.
      templine = line(1:sub_loc-1)//trim(current_macro%value)// &
                 line(sub_loc+len_trim(current_macro%name):)
      line = templine
   end do
   current_macro => current_macro%next
end do

return
end subroutine substitutions

!          --------
subroutine is_equal(name,value,is_eq)
!          --------
character(len=*), intent(in) :: name,value
logical, intent(out) :: is_eq

type(macro), pointer :: current_macro

! determines if name is in the list of macro names with the given value

current_macro => macro_list
do
   if (.not.associated(current_macro)) then
      is_eq = .false.
      exit
   end if
   if (current_macro%name == name) then
      if (current_macro%value == value) then
         is_eq = .true.
      else
         is_eq = .false.
      end if
      exit
   end if
   current_macro => current_macro%next
end do

return
end subroutine is_equal

!          ----------
subroutine is_defined(name,is_def)
!          ----------
character(len=*), intent(in) :: name
logical, intent(out) :: is_def

type(macro), pointer :: current_macro

! determines if name is in the list of macro names

current_macro => macro_list
do
   if (.not.associated(current_macro)) then
      is_def = .false.
      exit
   end if
   if (current_macro%name == name) then
      is_def = .true.
      exit
   end if
   current_macro => current_macro%next
end do

return
end subroutine is_defined

!          ---------
subroutine add_macro(name,value)
!          ---------
character(len=*), intent(in) :: name,value
type(macro), pointer :: newmacro

! add name to the list of defined macros with the given value

allocate(newmacro)
newmacro%next => macro_list
newmacro%name = name
newmacro%value = value
macro_list => newmacro

return
end subroutine add_macro

!          ---------------
subroutine init_macro_list()
!          ---------------

! initialize the linked list of macros to empty

nullify(macro_list)

return
end subroutine init_macro_list

end module sppr_mod

!-----------
program sppr
!-----------

! simple preprocessor.  This is very simple, makes assumptions that I know I
! will satisfy, and is not robust.  It is not intended for general use.  For
! a general purpose Fortran 90 preprocessor, get Michel Olagnon's f90ppr
! from the Fortran Market at http://www.fortran.com/fortran/market.html

use sppr_mod
implicit none

character(len=132) :: line
integer :: ioflag, directive, op, orig_len, if_depth, if_true_depth
logical :: skipping, including, is_eq, is_def
character(len=MAXCHAR) :: leftop, rightop

! set initial state

skipping = .false.
including = .false.
if_depth = 0
if_true_depth = 0
call init_macro_list()

! inifinte loop to read and processor lines of input.  Exits on an error

input_loop: do

! get next input line

   if (including) then
! input from an include file
      read(unit=10,fmt="(a)",iostat=ioflag) line
! at error, go back to standard input
      if (ioflag /= 0) then
         close(unit=10)
         including = .false.
         cycle input_loop
      end if
   else
! input from standard input
      read(unit=*,fmt="(a)",iostat=ioflag) line
      if (ioflag /= 0) then
         exit input_loop
      end if
   end if

! if the line begins with pound, get the directive

   if (line(1:1)=="#") then
      directive = get_directive(line)

! perform action depending on the directive

   select case(directive)

    case(IF_DIR)

      if_depth = if_depth + 1
      if (.not.skipping) then
! get the two strings and operator between them
         call get_operands(line,leftop,op,rightop)
! check equality of the operands and set skipping accordingly
         call is_equal(leftop,rightop,is_eq)
         if (op == EQUALS) then
            if (.not. is_eq) then
               skipping = .true.
               if_true_depth = if_depth - 1
            end if
         else if (op == NOT_EQUALS) then
            if (is_eq) then
               skipping = .true.
               if_true_depth = if_depth - 1
            end if
         end if
      end if


    case(IFDEF_DIR)

      if_depth = if_depth + 1
      if (.not.skipping) then
! get the operand
         call get_operands(line,leftop,op,rightop)
! check existence of the operand and set skipping accordingly
         call is_defined(leftop,is_def)
         if (.not.is_def) then
            skipping = .true.
            if_true_depth = if_depth - 1
         end if
      end if

    case(IFNDEF_DIR)

      if_depth = if_depth + 1
      if (.not.skipping) then
! get the operand
         call get_operands(line,leftop,op,rightop)
! check existence of the operand and set skipping accordingly
         call is_defined(leftop,is_def)
         if (is_def) then
            skipping = .true.
            if_true_depth = if_depth - 1
         end if
      end if

    case(ELSE_DIR)

! reverse the condition of skipping
      if (skipping) then
         if (if_true_depth == if_depth - 1) then
            skipping = .false.
            if_true_depth = if_depth
         end if
      else
         skipping = .true.
         if_true_depth = if_depth - 1
      end if

    case(ENDIF_DIR)

! terminate skipping
      if_depth = if_depth - 1
      if (if_true_depth == if_depth) then
         skipping = .false.
      end if

    case(INCLUDE_DIR)

      if (.not.skipping) then
! get the operand (filename surrounded by double quotes)
         call get_operands(line,leftop,op,rightop)
! open the file for reading
         open(unit=10,action="READ",file=leftop(2:len_trim(leftop)-1), &
              iostat=ioflag,status="UNKNOWN")
! if no error on opening, change the input source
         if (ioflag == 0) then
            including = .true.
         end if
      end if

    case(DEFINE_DIR)

      if (.not.skipping) then
! get the operands
         call get_operands(line,leftop,op,rightop)
! add the macro to the macro list
         call add_macro(leftop,rightop)
      end if

    case default

! do nothing if the directive is not recognized

   end select

! if the line does not start with pound ...
   else

! if we are not currently skipping code
! make macro substitutions and output the line
      if (.not.skipping) then
         orig_len = len_trim(line)
         call substitutions(line)
! don't output lines that became empty from substitutions
         if (len_trim(line) /= 0 .or. orig_len == 0) then
            write(unit=*,fmt="(1x,a)") trim(line)
         end if
      end if

   end if ! line starts with pound

end do input_loop

! Don't let Microsoft execute the stop statement, because it will 
! write "Program Terminated" to standard out 
 
!MS$IF .NOT. DEFINED (_MSFORTRAN_) 
stop 
!MS$ENDIF
end program sppr
