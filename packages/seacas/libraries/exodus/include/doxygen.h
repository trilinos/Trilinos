/* clang-format off */

/*!  \mainpage Exodus API Documentation

\section intro Introduction

EXODUS is the successor of the widely used finite element (FE) data file format EXODUS
(henceforth referred to as EXODUS I) developed by Mills-Curran and Flanagan. It
continues the concept of a common database for multiple application codes (mesh generators,
analysis codes, visualization software, etc.) rather than code-specific utilities, affording
flexibility and robustness for both the application code developer and application code user.
By using the EXODUS data model, a user inherits the flexibility of using a large array of
application codes (including vendor-supplied codes) which access this common data file
directly or via translators.

The uses of the EXODUS data model include the following:
    - Problem definition -- mesh generation, specification of locations of boundary conditions and
load application, specification of material types.
    - Simulation -- model input and results output.
    - Visualization -- model verification, results postprocessing, data interrogation, and analysis
tracking.

\section avail Availability

The Exodus library source code is available on Github at
https://github.com/gsjaardema/seacas

For bug reports, documentation errors, and enhancement suggestions, contact:
- Gregory D. Sjaardema
- WEB:   https://github.com/gsjaardema/seacas/issues
- EMAIL: gdsjaar@sandia.gov
- EMAIL: gsjaardema@gmail.com
- PHONE: (505) 844-2701 (office)

\section license License
The EXODUS library is licensed under the BSD open source license.

     Copyright(C) 1999-2020 National Technology & Engineering Solutions
     of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
     NTESS, the U.S. Government retains certain rights in this software.

     See packages/seacas/LICENSE for details

     Redistribution and use in source and binary forms, with or without
     modification, are permitted provided that the following conditions are
     met:

     * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.

     * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.

     * Neither the name of NTESS nor the names of its
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.

     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
     "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
     LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
     A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
     OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
     SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
     LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
     DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
     THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
     (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
     OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

\section devel Development of EXODUS

The evolution of the EXODUS data model has been steered by FE application code developers
who desire the advantages of a common data format. The EXODUS model has been
designed to overcome deficiencies in the EXODUS I file format and meet the following
functional requirements as specified by these developers:
   - Random read/write access.
   - Application programming interface (API) -- provide routines callable from FORTRAN, C, and C++
application codes.
   - Extensible -- allow new data objects to be added without modifying the application programs
that use the file format.
   - Machine independent -- data should be independent of the machine which generated it.
   - Real-time access during analysis -- allow access to the data in a file while the file is
being created.

To address these requirements, the open source database library
NetCDF (http://www.unidata.ucar.edu/software/netcdf/) was selected to handle the low-level data
storage. The EXODUS
II library functions provide the mapping between FE data objects and
NetCDF dimensions, attributes, and variables. Thus, the code developer
interacts with the data model using the vocabulary of an FE analyst
(element connectivity, nodal coordinates, etc.) and is relieved of the
details of the data access mechanism.

Because an EXODUS file is a NetCDF file, an application program can
access data via the EXODUS API or the NetCDF API directly. Although
accessing the data directly via the NetCDF API requires more in-depth
understanding of NetCDF, this capability is a powerful feature that
allows the development of auxiliary libraries of special purpose
functions not offered in the standard EXODUS library. For example,
if an application required access to the coordinates of a single node
(the standard library function returns the coordinates for all of the
nodes in the model), a simple function could be written that calls
NetCDF routines directly to read the data of interest.

\section descrip Description of Data Objects

The data in EXODUS files can be divided into three primary
categories: initialization data, model, and results.

Initialization data includes sizing parameters (number of nodes,
number of elements, etc.), optional quality assurance information
(names of codes that have operated on the data), and optional
informational text.

The model is described by data which are static (do not change through
time). These data include nodal coordinates, element connectivity
(node lists for each element), element attributes, and node sets and
side sets (used to aid in applying loading conditions and boundary
constraints).

The results are optional and include five types of variables -- nodal,
element, nodeset, sideset, and global -- each of which is stored
through time. Nodal results are output (at each time step) for all the
nodes in the model. An example of a nodal variable is displacement in
the X direction. Element, nodeset, and sideset results are output (at
each time step) for all entities (elements, nodes, sides) in one or
more entity block. For example, stress may be an element
variable. Another use of element variables is to record element status
(a binary flag indicating whether each element is "alive" or "dead")
through time. Global results are output (at each time step) for a
single element or node, or for a single property. Linear momentum of a
structure and the acceleration at a particular point are both examples
of global variables.  Although these examples correspond to typical FE
applications, the data format is flexible enough to accommodate a
spectrum of uses.

A few conventions and limitations must be cited:

 - There are no restrictions on the frequency of results output except
 that the time value associated with each successive time step must
 increase monotonically.
 - To output results at different frequencies (i.e., variable A at
 every simulation time step, variable B at every other time step)
 multiple EXODUS files must be used.
 - There are no limits to the number of each type of results, but once
 declared, the number cannot change.
 - If the mesh geometry or topology changes in time (i.e., number of
 nodes increases, connectivity changes), then the new geometry must be
 output to a new EXODUS file.

\section int64 Integer Bulkdata Storage Details (32-bit and 64-bit integers)

The EXODUS database can store integer bulk data, entity map data, and
mesh entity (block/set) ids in either 32-bit or 64-bit integer format. The data
considered "bulk data" are:

 - element, face, and edge connectivity lists,
 - element, face, edge, and node set entity lists,

The entity map data is any data stored in one of the 'map' objects on
the exodus file.  This includes:
 - id maps
 - number maps
 - order maps
 - processor node maps
 - processor element maps.

A mesh entity id is the id of any block (element block, edge block,
...); set (node set, face set, ...), coordinate frame, and
communication map.

When an EXODUS file is created via the ex_create() function, the
'mode' argument provides the mechanism for specifying how integer data
will be passed as arguments to the API functions and also how the
integer data will be stored on the database. The ex_open() function
also provides a mechanism for specifying how integer data will be
passed as arguments.

The method uses the 'mode' argument to the ex_open() and
ex_create() functions.  The mode is a 32-bit integer in which certain
bits are turned on by or'ing certain predefined constants.

    exoid = ex_create( "test.exo",
                       EX_CLOBBER|EX_MAPS_INT64_DB|EX_MAPS_INT64_API,
                       &appWordSize, &diskWordSize );

The constants related to the integer size (32-bit or 64-bit)
specification are:

|   Constant Name    | Which data are 64-bit
---------------------|----------------------
| #EX_MAPS_INT64_DB   | entity map data
| #EX_IDS_INT64_DB    | mesh entity ids
| #EX_BULK_INT64_DB   | bulk data
| #EX_ALL_INT64_DB    | (the above 3 or'd together)
| #EX_MAPS_INT64_API  | entity map data
| #EX_IDS_INT64_API   | mesh entity ids
| #EX_BULK_INT64_API  | bulk data
| #EX_INQ_INT64_API   | integers passed to/from ex_inquire()
| #EX_ALL_INT64_API   | (the above 4 or'd together)

The constants that end with `_DB` specify that that particular integer
data is stored on the database as 64-bit integers; the constants that
end with `_API` specify that that particular integer data is passed
to/from API functions as 64-bit integers.

If the range of the data being transmitted is larger than the
permitted integer range (for example, if the data is stored on the
database as 64-bit ints and the application specifies passing data as
32-bit ints), the API function will return an error.

The three types of integer data whose storage can be specified are
- maps (`EX_MAPS_INT64_`),
- "bulk data" including connectivity lists and entity lists (`EX_BULK_INT64_`), and
- entity ids which are the ids of element, face, edge, and node sets
   and blocks; and map ids (`EX_IDS_INT64_`)

The function ex_int64_status()(exoid) is used to determine the integer
storage types being used for the EXODUS database `exoid`.  It returns
an integer which can be and'ed with the above flags to determine
either the storage type or function parameter type.

For example, if
`(#EX_MAPS_INT64_DB & ex_int64_status()(exoid))` is true, then map data is
being stored as 64-bit integers for that database.

It is not possible to determine the integer data size on a database
without opening the database via an ex_open() call. However, the
integer size specification for API functions can be changed at any
time via the ex_set_int64_status()(exoid, mode) function. The mode is
one or more of `#EX_MAPS_INT64_API`, `#EX_IDS_INT64_API`, or
`#EX_BULK_INT64_API`, or'd together.  Any exodus function calls after
that point will use the specified integer size. Note that a call to
ex_set_int64_status()(exoid, mode) overrides any previous setting for
the integer sizes used in the API.  The ex_create() function is the
only way to specify the integer sizes specification for database
integers.

\subsection int64_fortran_api Fortran API
The fortran API is uses the same mechanism as was described above for
the C API. If using the "8-byte real and 8-byte int" fortran mode
typically used by the SEACAS applications (the compiler automatically
promotes all integers and reals to 8-byte quantities), then the
fortran exodus library will automatically enable the `EX_*_INT64_API`
options; the client still needs to specify the `EX_*_INT64_DB` options.

\subsection int64_fortran_imp Fortran Implementation

The new capability to pass 64-bit integer data through the fortran and
C API functions simplifies the implementation of the "8-byte real
8-byte int" usage of the exodus library. Previously, the wrapper
routines in addrwrap.F were required to convert the 8-byte integer
data on the client side to/from 4-byte integers on the library
side. This required extra memory allocation and complications that are
now handled at the lowest level in the NetCDF library.  The
functions in the fortran API have all been converted to
pass 64-bit integers down to the C API which has removed some code and
simplified those functions.


\section db_options Database Options (Compression, Name Length, File Type)

The ex_set_option() function call is used to set various options on the
database.  Valid values for 'option' are:

|   Option Name          | Option Values |
-------------------------|---------------|
| #EX_OPT_MAX_NAME_LENGTH | Maximum length of names that will be returned/passed via API call. |
| #EX_OPT_COMPRESSION_TYPE | Not currently used; default is gzip |
| #EX_OPT_COMPRESSION_LEVEL | In the range [0..9]. A value of 0 indicates no compression |
| #EX_OPT_COMPRESSION_SHUFFLE | 1 if enabled, 0 if disabled |
| #EX_OPT_INTEGER_SIZE_API | 4 or 8 indicating byte size of integers used in API functions. |
| #EX_OPT_INTEGER_SIZE_DB  | Query only, returns 4 or 8 indicating byte size of integers stored on the database. |

The compression-related options are only available on NetCDF-4 files
since the underlying hdf5 compression functionality is used for the
implementation. The compression level indicates how much effort should
be expended in the compression and the computational expense increases
with higher levels; in many cases, a compression level of 1 is
sufficient.

\section names Variable, Attribute, and Entity Block/Set Names
The length of the Variables, Attributes, and Entity Block/Set names is
variable.  The default length is 32 characters to provide backward
compatibility.  This is the default on both read and write, so if
there is a database with longer names and the reader does not change
the length of names to be returned, any API call that returns a name
will truncate the name at 32 characters.

To avoid this, the reading application can all
~~~{.c}
  // Determine maximum length of names stored on database
  int max_name_length = ex_inquire_int(exoid, EX_INQ_DB_MAX_USED_NAME_LENGTH);

  // Tell the library to return names this length
  ex_set_max_name_length(exodusFilePtr, max_name_length);
~~~

On write, you can call:

~~~{.c}
   ex_set_option(exoid, EX_OPT_MAX_NAME_LENGTH, {max_name_length});

   // or equivalently
   ex_set_max_name_length(exoid, {max_name_length});
~~~

which tells the database that you will be using names of that length or shorter.

Following this call, you can define (i.e., read/write) names of any
size; if the names are longer than `{max_name_length}`, then they will be truncated otherwise they will pass through unchanged.

There are three queries that can be made to ex_inquire() or
ex_inquire_int():

  - #EX_INQ_DB_MAX_ALLOWED_NAME_LENGTH -- returns the value of the
    maximum size that can be specified for `max_name_length`
    (netcdf/hdf5 limitation)
  - #EX_INQ_DB_MAX_USED_NAME_LENGTH -- returns the size of the longest
    name on the database.
  - #EX_INQ_MAX_READ_NAME_LENGTH -- returns the maximum name length
    size that will be passed back to the client. 32 by default,
    set by the previously mentioned ex_set_option() or
    ex_set_max_name_length() call.

\note
  - The length of the QA records (ex_get_qa(), ex_put_qa()) is not affected by this setting and each entry in the QA record is still limited to 32 characters.
  - The length of the `entity_descrip` type passed and returnen in the
    ex_get_block() and ex_put_block() calls is still limited to 32 characters.
  - The length of the title is limited to 80 characters
    (ex_get_init(), ex_get_init_ext(), ex_put_init(), ex_put_init_ext()).
  - The length of the info records is limited to 80 characters
  (ex_put_info(), ex_get_info()).

\defgroup ResultsData Results Data
@{
This section describes data file utility functions for creating
opening a file, initializing a file with global parameters, reading
writing information text, inquiring on parameters stored in the data
file, and error reporting.

The results are optional and include an optional variable type for
each block and set type (node, edge, face, and element) in addition
there are global variables and sideset variables -- each of which is
stored through time. Nodal results are output (at each time step) for
all the nodes in the model. An example of a nodal variable is
displacement in the X direction. Global results are output (at each
time step) for a single element or node, or for a single
property. Linear momentum of a structure and the acceleration at a
particular point are both examples of global variables.  The other
results are output (at each time step) for all entities (elements,
faces, edges, nodes, or sides) in one or more entity blocks. For
example, stress may be an element variable. Another use of element
variables is to record element status (a binary flag indicating
whether each element is "alive" or "dead") through time. Although
these examples correspond to typical FE applications, the data format
is flexible enough to accommodate a spectrum of uses.

A few conventions and limitations must be cited:

+ There are no restrictions on the frequency of results output except
that the time value associated with each successive time step should
increase monotonically.

+ All variables are output at the same time frequency. To output
results at different frequencies (i.e., variable A at every simulation
time step, variable B at every other time step) multiple files must be
used.

+ There are no limits to the number of each type of results, but once
declared, the number cannot change.

+ If the mesh geometry changes in time (i.e., number of nodes
increases, connectivity changes), the new geometry must be output to a
new file.
@}

\defgroup Utilities Data File Utilities
  @{
This section describes data file utility functions for creating
opening a file, initializing a file with global parameters, reading
writing information text, inquiring on parameters stored in the data
file, and error reporting.
  @}

\defgroup ModelDescription Model Description
  @{
The routines in this section read and write information which
describe an exodus finite element model. This includes nodal
coordinates, element order map, element connectivity arrays,
element attributes, node sets, side sets, and object properties.
  @}

@example ../test/CreateEdgeFace.c
@example ../test/ReadEdgeFace.c
@example ../test/create_mesh.c
@example ../test/rd_wt_mesh.c
@example ../test/test-empty.c
@example ../test/test_nemesis.c
@example ../test/test_ts_errval.c
@example ../test/test_ts_files.c
@example ../test/test_ts_nvar.c
@example ../test/test_ts_nvar_rd.c
@example ../test/test_ts_partial_nvar.c
@example ../test/test_ts_partial_nvar_rd.c
@example ../test/testcp.c
@example ../test/testcp_nl.c
@example ../test/testcp_tran.c
@example ../test/testcpd.c
@example ../test/testrd-groups.c
@example ../test/testrd-long-name.c
@example ../test/testrd-nfaced.c
@example ../test/testrd-nsided.c
@example ../test/testrd.c
@example ../test/testrd1.c
@example ../test/testrd_nc.c
@example ../test/testrd_par.c
@example ../test/testrd_ss.c
@example ../test/testrdd.c
@example ../test/testrdwt.c
@example ../test/testwt-compress.c
@example ../test/testwt-groups.c
@example ../test/testwt-long-name.c
@example ../test/testwt-nface-nside.c
@example ../test/testwt-nfaced.c
@example ../test/testwt-nsided.c
@example ../test/testwt-one-attrib.c
@example ../test/testwt-oned.c
@example ../test/testwt-partial.c
@example ../test/testwt-results.c
@example ../test/testwt-zeroe.c
@example ../test/testwt-zeron.c
@example ../test/testwt.c
@example ../test/testwt1.c
@example ../test/testwt2.c
@example ../test/testwt_clb.c
@example ../test/testwt_nc.c
@example ../test/testwt_nossnsdf.c
@example ../test/testwt_ss.c
@example ../test/testwtd.c
@example ../test/testwtm.c
@example ../test/twod.c

@example ../exodus_for/test/test_nem.f
@example ../exodus_for/test/testcp.f
@example ../exodus_for/test/testcpd.f
@example ../exodus_for/test/testcpnl.f
@example ../exodus_for/test/testrd.f
@example ../exodus_for/test/testrd1.f
@example ../exodus_for/test/testrd_nsid.f
@example ../exodus_for/test/testrdd.f
@example ../exodus_for/test/testwt.f
@example ../exodus_for/test/testwt1.f
@example ../exodus_for/test/testwt2.f
@example ../exodus_for/test/testwt3.f
@example ../exodus_for/test/testwt_nsid.f
@example ../exodus_for/test/testwtd.f
@example ../exodus_for/test/testwtm.f
*/

/* clang-format on */
