// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef PYTRILINOS_DAP_HPP
#define PYTRILINOS_DAP_HPP

// Include the NumPy and Python headers
#define NO_IMPORT_ARRAY
#include "numpy_include.hpp"
#include "Python3Compat.hpp"

// Teuchos includes
#include "Teuchos_Array.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Tuple.hpp"

////////////////////////////////////////////////////////////////////////

// #define PYTRILINOS_DAP_VERBOSE

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

#if defined(PyTrilinos_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The PyTrilinos package is deprecated"
#endif
#endif

