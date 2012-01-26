/*!  \mainpage ExodusII API Documentation

\section intro Introduction

EXODUS II is the successor of the widely used finite element (FE) data file format EXODUS
(henceforth referred to as EXODUS I) developed by Mills-Curran and Flanagan. It
continues the concept of a common database for multiple application codes (mesh generators,
analysis codes, visualization software, etc.) rather than code-specific utilities, affording
flexibility and robustness for both the application code developer and application code user.
By using the EXODUS II data model, a user inherits the flexibility of using a large array of
application codes (including vendor-supplied codes) which access this common data file
directly or via translators.

The uses of the EXODUS II data model include the following:
    - Problem definition -- mesh generation, specification of locations of boundary conditions and load application, specification of material types.
    - Simulation -- model input and results output.
    - Visualization -- model verification, results postprocessing, data interrogation, and analysis tracking.

\section avail License and Availability
The EXODUS II library is licensed under the BSD open source license.

Copyright (c) 2005 Sandia Corporation. Under the terms of Contract DE-AC04-94AL85000
with Sandia Corporation, the U.S. Government retains certain rights in this software.

Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
  - Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
  - Redistributions in binary form must reproduce the above copyright notice, this list
of conditions and the following disclaimer in the documentation and/or other
materials provided with the distribution.
  -Neither the name of Sandia Corporation nor the names of its contributors may be
used to endorse or promote products derived from this software without specific
prior written permission.

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

The ExodusII library source code is available on Sourceforge at
http://sourceforge.net/projects/exodusii

For bug reports, documentation errors, and enhancement suggestions, contact:
- Gregory D. Sjaardema
- PHONE: (505) 844-2701
- EMAIL: gdsjaar@sandia.gov

\section devel Development of EXODUS II

The evolution of the EXODUS data model has been steered by FE application code developers
who desire the advantages of a common data format. The EXODUS II model has been
designed to overcome deficiencies in the EXODUS I file format and meet the following
functional requirements as specified by these developers:
   - Random read/write access.
   - Application programming interface (API) -- provide routines callable from FORTRAN, C, and C++ application codes.
   - Extensible -- allow new data objects to be added without modifying the application programs that use the file format.
   - Machine independent -- data should be independent of the machine which generated it.
   - Real-time access during analysis -- allow access to the data in a file while the file is
being created.

To address these requirements, the open source database library
etCDF (http://www.unidata.ucar.edu/software/netcdf/) was selected to handle the low-level data storage. The EXODUS
II library functions provide the mapping between FE data objects and
netCDF dimensions, attributes, and variables. Thus, the code developer
interacts with the data model using the vocabulary of an FE analyst
(element connectivity, nodal coordinates, etc.) and is relieved of the
details of the data access mechanism. 

Because an EXODUS II file is a netCDF file, an application program can
access data via the EXODUS II API or the netCDF API directly. Although
accessing the data directly via the netCDF API requires more in-depth
understanding of netCDF, this capability is a powerful feature that
allows the development of auxiliary libraries of special purpose
functions not offered in the standard EXODUS II library. For example,
if an application required access to the coordinates of a single node
(the standard library function returns the coordinates for all of the
nodes in the model), a simple function could be written that calls
netCDF routines directly to read the data of interest.

\section descrip Description of Data Objects

The data in EXODUS II files can be divided into three primary
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
 multiple EXODUS II files must be used.
 - There are no limits to the number of each type of results, but once
 declared, the number cannot change.
 - If the mesh geometry or topology changes in time (i.e., number of
 nodes increases, connectivity changes), then the new geometrymust be
 output to a new EXODUS II file.

\defgroup ResultsData Results Data
@{
 This section describes data file utility functions for creating /
 opening a file, initializing a file with global parameters, reading /
 writing information text, inquiring on parameters stored in the data
 file, and error reporting.
@}

\defgroup Utilities Data File Utilities
  @{
This section describes data file utility functions for creating /
opening a file, initializing a file with global parameters, reading /
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


*/
