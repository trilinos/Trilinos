program structsize
use lb_user_const
type(LB_GID) :: a(2)
call c_structsize(a(1),a(2))
end program structsize
