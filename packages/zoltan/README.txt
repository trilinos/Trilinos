lriesen@sandia.gov
July 1, 2010

WHAT THIS IS:

This is the Zoltan branch "zoltan_gid_64".

It is Zoltan 3.3, with modifications to permit 64 bit global IDs, and counts of 
global objects that require 64 bit integers.

Bug fixes subsequent to Zoltan 3.3 have been applied to this branch.

To checkout this branch from the Trilinos git repository:

git clone software.sandia.gov:/space/git/Trilinos
cd Trilinos
git branch --track zoltan_gid_64 origin/zoltan_gid_64
git checkout zoltan_gid_64
cd packages/zoltan

HOW IT IS DIFFERENT FROM THE LATEST ZOLTAN RELEASE:

The Zoltan library is compiled from C-language code that uses the data type
"int" for the global number of objects being partitioned.  This limits the
number of objects to about 2 billion on most architectures.  In this version
of Zoltan, global counts are data type "ssize_t" (signed "size_t") which
will be a 64 bit integer on a 64-bit system and a 32-bit integer on a
32-bit system.

In addition, a global ID is a user-defined multiple of a ZOLTAN_ID_TYPE
which is an "unsigned int".  While the application could use global IDs composed 
of two ZOLTAN_ID_TYPEs to represent a space of over 2 billion objects, this
would be inconvenient.  So the data type for ZOLTAN_ID_TYPE is set at
compile time.  It defaults to "long", but can be set to "int" or "long long"
with a configure script option.

There is a new library function:

  int Zoltan_get_global_id_type(char **name)

which returns the number of bytes in a ZOLTAN_ID_TYPE.  If "name" is non-NULL,
it also returns the name of the data type in "name".  This would be one of
"int", "long long", or "long".  This function is called by the tests in the
directory tests/Large_Data to ensure that the tests are compiled with the same
ZOLTAN_ID_TYPE that the Zoltan library was compiled with.

HOW TO USE IT:

Zoltan is part of the Trilinos framework of linear algebra solvers and tools.  
Trilinos uses CMake to configure and build its packages.  Zoltan uses both
CMake and autoconf to configure and build.  Only the autoconf build should
be used with "zoltan_gid_64".

So to configure, in your build directory:

  {srcdir}/configure --enable-mpi --with-id-type={int|long|llong}

The default ID type is "long".  See "{srcdir}/configure --help" for more options
that may be useful.

Then, in your build directory:

  make

You will find libzoltan.a in your build directory's "src" directory.

You can test your build with:

  make check

which runs some tests.

The API is unchanged, with the exception of the addition of "Zoltan_get_global_id_type"
function, so you should be able to compile your existing code after changing, if 
necessary, the data type you are using for global IDs.

LIMITATIONS:

The library uses data type "int" or "unsigned int" for local IDs and for counts
of local objects, objects being sent, and objects being received.  This means
that at no time can a single process own more than 2^31 objects (about 2.1 billion).

Coloring has not been modified to work with 64 bit global IDs.  You may see 
compiler warnings because of this.

The work required to make this version of Zoltan work with 64 bit third party libraries
is incomplete.  Let us know if you want to use ParMetis or Scotch with 64 bit indexes.

INFORMATION FOR DEVELOPERS:

Data type:

The following definitions are created at preprocessor-time in include/zoltan_types.h and
at run-time in zz/zz_util.c:Zoltan_set_mpi_types() which is called from Zoltan_Initialize().

ZOLTAN_GNO_TYPE - the data type for global counts (ssize_t)
ZOLTAN_GNO_SPECIFIER - the scanf characters for a ZOLTAN_GNO_TYPE, you can also just use "zd"
ZOLTAN_GNO_MPI_TYPE - the MPI_Datatype to use when sending/receiving ZOLTAN_GNO_TYPEs

ZOLTAN_ID_TYPE - self explanatory
ZOLTAN_ID_SPECIFIER - the scanf characters for a ZOLTAN_ID_TYPE
ZOLTAN_ID_MPI_TYPE - the MPI_Datatype to use when sending/receiving ZOLTAN_ID_TYPEs

If you are not using "configure" to build Zoltan, you can instead use the compile-time
flags -DZOLTAN_ID_TYPE_INT, -DZOLTAN_ID_TYPE_LONG, or -DZOLTAN_ID_TYPE_LONG_LONG.

intptr_t - use this type for an integer that may hold a pointer

Zoltan_Map_Create():

In the Zoltan trunk, Zoltan_Map_Create() assumes it is handling keys that are
multiples of ZOLTAN_ID_TYPEs.  In this branch the function is provided with
the number of bytes in the key, not the number of ZOLTAN_ID_TYPEs.

Memory:

Because testing of this branch involves running large memory problems, I added the
function Zoltan_write_linux_meminfo() which will write out the contents of /proc/meminfo
on a Linux machine.  The new function Zoltan_Memory_Get_Debug() returns the debug level
set in mem.c by Zoltan_Memory_Debug().  zdrive has a new input option

  zoltan memory debug level = n

which will set the debug level.  Then after partitioning, zdrive checks the debug level
and if there was an error and it is running on a linux machine it will dump out /proc/meminfo.

I modified the configure script to define HOST_LINUX on a linux machine.

I wrote three tests in tests/Large_Data that test PHG, RCB and RIB with arbitrarily large
numbers of objects.  They have signal handlers that call Zoltan_write_linux_meminfo() on
a Linux machine.

Code changes:

The hypergraph coarsening and refinement code was extensively rewritten in places where
results were shared across rows and columns, because the code assumed that pointers, floats,
and integral data were all the same size.  Merging new work from the trunk into the code
will require careful inspection first.

INTERESTING CHART:

32 and 64 bit data models (ILP - integer/long/pointer):

type            LP32    ILP32   ILP64   LLP64   LP64

char            8       8       8       8       8
short           16      16      16      16      16
_int32                          32
int             16      32      64      32      32
long            32      32      64      32      64
long long                               64
pointer         32      32      64      64      64

ILP32 is most widely used.
LP64 is most widely used.

LLP64 is ILP32 with new 64 bit int added to it - used for Win64.
