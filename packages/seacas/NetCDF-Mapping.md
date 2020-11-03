# Copyright(C) 1999-2020 National Technology & Engineering Solutions
# of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
# NTESS, the U.S. Government retains certain rights in this software.
#
# See packages/seacas/LICENSE for details
# Mapping of exodus entities onto NetCDF

If you are using NetCDF-4.5.1 or later, then you can ignore the information in this file.

The distributed version of netcdf sets the following limits
on dimensions and variables:
 * `#define NC_MAX_DIMS 1024`
 * `#define NC_MAX_VARS 8192`

For use with Exodus, it is recommended that these be increased to:
 * `#define NC_MAX_DIMS  65536`
 * `#define NC_MAX_VARS 524288`

The reason for these increases is due to the mapping of Exodus onto
NetCDF. The sections below show the number of Dimensions (controlled
by NC_MAX_DIMS) and Variables (controlled by NC_MAX_VARS) that are
used in an Exodus file.

## Entities
 * A mesh-entity is an individual node, edge, face, or element.
 * An entity is a set or block consisting of a single mesh-entity type.
 * Each entity can have variables, maps, and attributes which contain an entry per mesh-entity.
 * Each entity has an optional name and a required id (32-bit or 64-bit )which is non-negative.
 * A mesh-entity can be in one and only one entity block,
 * A mesh-entity can be in zero or more entity sets.
 * Currently there is only a single implicit node block containing all nodes in the model.

## Dimensions: (NC_MAX_DIMS)
* There are about 10 standard dimensions in every file.
* plus one for each set plus one if any attributes
* plus two for each block plus one if any attributes
* plus one for each transient variable on an entity (node, node set, element block, element set, ...)

## Variables: (NC_MAX_VARS)
* There are a few standard dimensions
  * times
  * names of each entity type (block set)
  * ids of each entity type (block set)
  * status of each entity type (block set)
  * #ndim coordinates (1,2,3)

* Each block adds 1 + 2*#attr_on_block + #var_on_block

* Each set adds 2 + 2*#attr_on_set + #var_on_set

* Each sideset add 3 + 2*#attr_on_sset + #var_on_sset

* Each map adds 1

## Example
If we have an exodus file with:
 * Nodes

 * 5 Element blocks
   * 4 transient variables per element block
   * 2 attributes per element block

 * 4 Node Sets
   * Distribution Factors defined on each set
   * 3 transient variables

 * 3 Side Sets
   * Distribution Factors defined on each set
   * 2 transient variables

Then there would be about:
 `10 + 5*(2+1) + 4*(2) + 3*(2) + 1 + 1 + 1 = 42` Dimensions

There would be about:
 `5*(1+2*2+4) + 4*(2+3) + 3*(3+2) + 3*(5+4+3) + 3 + 1 = 120` Variables.

From this, you can see that a moderately complicated model would
quickly overflow the standard values for `NC_MAX_DIMS` and `NC_MAX_VARS`.
