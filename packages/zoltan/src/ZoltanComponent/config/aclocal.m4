dnl A way to find MakeInc.CCA_Component given a list of directories
dnl   Usage AC_FIND_FILE_STRICT_JR( filename, "dir1 dir2 dir3 ..")
dnl   On finding the directory where the file is ac_find_file_dir is set to the directory
dnl   else ac_find_file_strict_dir is set to not_found.

AC_DEFUN([AC_FIND_FILE_STRICT_JR], 
[
   ac_find_file_strict_dir="not_found"
   for dir in "$2" ; do
	if test -d "$dir" ; then
	   if test -f "$dir"/"$1" ; then
  	      ac_find_file_strict_dir="$dir" ;  break 3 ;
	   fi
        fi
   done

   if test "$ac_find_file_strict_dir" != "not_found" ; then
  	AC_MSG_RESULT( [found it in $ac_find_file_strict_dir] ) 
   else
	AC_MSG_ERROR( [Did not find "$1". Exiting] )
   fi
])
    
dnl This is the same as before, but does not exit if the file is not found
AC_DEFUN([AC_FIND_FILE_NOT_STRICT_JR], 
[
   ac_find_file_not_strict_dir="not_found"
   for dir in "$2" ; do
	if test -d "$dir" ; then
	   if test -f "$dir"/"$1" ; then
  	      ac_find_file_not_strict_dir="$dir" ;  break 3 ;
	   fi
        fi
   done

   if test "$ac_find_file_not_strict_dir" != "not_found" ; then
  	AC_MSG_RESULT( [found it in $ac_find_file_not_strict_dir] ) 
   fi
])
    