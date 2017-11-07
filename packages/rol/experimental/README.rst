XROL - Experimental Rapid Optimization Library 
----------------------------------------------

This "subpackage" contains experimental implementations of ROL 
components aimed at reducing runtime overhead and leveraging
new features that have been incorporated in the C++11 and 
C++14 standard. There are no concrete plans yet as to when
these features may be incorporated into ROL. Planning any
project depending on components in the experimenal directory
is *NOT RECOMMENDED* and no support or reliable behavior should
be expected. 

Components in this directory utilize a number of features 
added in the C++14 standard, so if you want to try using 
them, likely one of the following compilers will be needed:

gcc 5.0+
clang 3.4+ 
MSVC 19.0+ 
Intel 17.0+

Additionally CMake 3.1 or greater is needed to deduce compiler
C++14 compatibility. To build experimental tests requires
setting the CMake variable ROL_ENABLE_Experimental:BOOL=ON
