\page exodus_formats Exodus Formats

\section large Large Model (64-bit offset)

This is the current default format for exodus files.  In this format,
the size of an *individual* dataset within the file is limited
to 4 GiB (2 GiB for older versions of netcdf library).

Typically, the two largest datasets are an individual component of the
nodal coordinates and the connectivity for an element block.

An individual component of the nodal coordinates stored in double
precision limits the model to about 2^29 (537 Million) nodes.

The other dataset size is the connectivity of a single element
block. For 8-node hexes, the maximum number of elements in an element
block is 2^27 (134 Million) elements; for 27-node hexes, this is
reduced to 39.7 million elements per element block.

Since a complex model typically has many element blocks, the
node count is typically the controlling dataset, but
connectivity datasets can also be the limiting factor.

Note that the model can typically split a large element block
into two or more element blocks to avoid the connectivity
limit, but there is no way to split the node coordinate
dataset, so 2^29 is the hard upper limit on the number of
nodes stored in a Large Model (64-bit offset) format exodus
file.

\section normal Normal, or Classic

The non-transient portion of the model definition is limited to 2 GiB.
This means that the size of the coordinates, connectivtity, nodeset,
sideset, maps, and other datasets must total less than 2 GiB.  This
typically limits the model to about 35 million elements.

This is no longer supported for writing by the exodus library, but old
files in this format can still be read.

\section nc4 Netcdf-4 Classic

The netcdf-4 classic format uses HDF5 as the underlying file format
on disk.  This format eliminates the dataset size limits.  However,
since it only supports the use of 32-bit signed integers, the model
is limited to 2^31 nodes and 2^31 elements.

\section nc4nc Netcdf-4 Non-Classic

This netcdf-4 non-classic format also uses HDF5 as the underlying
file format on disk.  There are no dataset size limits and it
supports datasets storing 64-bit integers.  With this format, the
maximum model size should be unlimited.  Models well in excess of 2
Billion nodes and elements have been defined with this format
verifying that the previous limits do not exist in this format.
Currently when exodus uses this format, the only extension over the
classic model is the use of 64-bit integers in some of the
datasets. In the future, uses of additional netcdf-4 functionality
(groups, compound datatypes) is planned.

\section cdf5 CDF5

A PNetCDF developed format which provides 64-bit integer capability
in a non-HDF5-based file format. Being evaluated at this time.
