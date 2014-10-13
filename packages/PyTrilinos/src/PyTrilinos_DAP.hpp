// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef PYTRILINOS_DAP_H
#define PYTRILINOS_DAP_H

// Include the NumPy and Python headers
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"

// Teuchos includes
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Tuple.hpp"

////////////////////////////////////////////////////////////////////////

namespace PyTrilinos
{

////////////////////////////////////////////////////////////////////////

enum DistributionType
{
  NONE,
  BLOCK,
  CYCLIC,
  UNSTRUCTURED
};

////////////////////////////////////////////////////////////////////////

typedef Teuchos::ArrayRCP< int > IndicesType;

typedef Teuchos::Tuple< int, 2 > PaddingType;

////////////////////////////////////////////////////////////////////////

class DimensionDictionary
{
public:

  DimensionDictionary(PyObject * dim_dict);

  ~DimensionDictionary();

  DistributionType dist_type;
  int              size;
  int              proc_grid_size;
  int              proc_grid_rank;
  int              start;
  int              stop;
  IndicesType      indices;
  PaddingType      padding;
  int              block_size;
  bool             periodic;
  bool             one_to_one;

};

////////////////////////////////////////////////////////////////////////

class DistArrayProtocol
{
public:

  DistArrayProtocol(PyObject * distarray);

  ~DistArrayProtocol();

  int num_dims() const;
  DistributionType dist_type(int axis) const;
  int size(int axis) const;
  int proc_grid_size(int axis) const;
  int proc_grid_rank(int axis) const;
  int start(int axis) const;
  int stop(int axis) const;
  IndicesType indices(int axis) const;
  PaddingType padding(int axis) const;
  int block_size(int axis) const;
  bool periodic(int axis) const;
  bool one_to_one(int axis) const;
  PyObject * buffer() const;
  std::string version() const;

private:
  PyObject * __distarray__;
  std::string __version__;
  PyObject * __buffer__;
  Teuchos::Array< DimensionDictionary > dim_data;

  void checkAxis(int axis) const;
};

////////////////////////////////////////////////////////////////////////

}

#endif
