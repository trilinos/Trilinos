/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/

#ifndef _fei_impl_utils_hpp_
#define _fei_impl_utils_hpp_


#include <fei_macros.hpp>
#include <fei_fwd.hpp>
#include <fei_mpi.h>

#include <Teuchos_ParameterList.hpp>

#include <string>
#include <vector>
#include <map>


/** The fei namespace contains public functions, classes and interfaces.
*/
namespace fei {

/** The impl_utils namespace contains implementation-utilities. Helpers
  for implementation code, not part of the public API.
*/
namespace impl_utils {

/** Given a vector of sources and a vector of targets, fill the offsets vector
   such that if offsets[i] >= 0, then sources[i] == targets[offsets[i]].
   For any i such that sources[i] is not found in targets, then offsets[i] == -1.
   The offsets vector will have length equal to the length of sources.

   For efficiency, it is assumed that both sources and targets are sorted.
*/
void find_offsets(const std::vector<int>& sources,
                  const std::vector<int>& targets,
                  std::vector<int>& offsets);

/** pack an fei::FillableMat object into a pair of std::vector objects.
*/
void pack_FillableMat(const fei::FillableMat& mat,
                      std::vector<int>& intdata,
                      std::vector<double>& doubledata);

void pack_FillableMat(const fei::FillableMat& mat,
                      std::vector<char>& buffer,
                      bool resize_buffer=true);

/** unpack a pair of std::vector objects into an fei::FillableMat object.
    The std::vector objects are assumed to have been produced by the
    function pack_FillableMat(...).
*/
void unpack_FillableMat(const std::vector<int>& intdata,
                        const std::vector<double>& doubledata,
                        fei::FillableMat& mat,
                        bool clear_mat_on_entry = true,
                        bool overwrite_entries = true);

void unpack_FillableMat(const std::vector<char>& intdata,
                        fei::FillableMat& mat,
                        bool clear_mat_on_entry = true,
                        bool overwrite_entries = true);

/** return-value is true if the unpacked matrix contains nothing but zeros...*/
bool unpack_CSRMat(const std::vector<char>& buffer, fei::CSRMat& mat);

void pack_indices_coefs(const std::vector<int>& indices,
                        const std::vector<double>& coefs,
                        std::vector<char>& buffer,
                        bool resize_buffer=true);

void unpack_indices_coefs(const std::vector<char>& buffer,
                          std::vector<int>& indices,
                          std::vector<double>& coefs);

void separate_BC_eqns(const fei::FillableMat& mat,
                    std::vector<int>& bcEqns,
                    std::vector<double>& bcVals);

void create_col_to_row_map(const fei::FillableMat& mat,
                           std::multimap<int,int>& crmap);

int remove_couplings(fei::FillableMat& mat);

void global_union(MPI_Comm comm,
                  const fei::FillableMat& localMatrix,
                  fei::FillableMat& globalUnionMatrix);

void global_union(MPI_Comm comm,
                  const fei::CSVec& localVec,
                  fei::CSVec& globalUnionVec);

void translate_to_reduced_eqns(const fei::Reducer& reducer, fei::CSRMat& mat);

void translate_to_reduced_eqns(const fei::Reducer& reducer, fei::CSVec& vec);

void add_to_graph(const fei::CSRMat& inmat, fei::Graph& graph);

void add_to_matrix(const fei::CSRMat& inmat, bool sum_into, fei::Matrix& matrix);

}//namespace impl_utils
}//namespace fei

#endif

