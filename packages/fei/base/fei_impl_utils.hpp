#ifndef _fei_impl_utils_hpp_
#define _fei_impl_utils_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2008 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

size_t num_bytes_FillableMat(const fei::FillableMat& mat);

void pack_FillableMat(const fei::FillableMat& mat, char* buffer);

/** unpack a pair of std::vector objects into an fei::FillableMat object.
    The std::vector objects are assumed to have been produced by the
    function pack_FillableMat(...).
*/
void unpack_FillableMat(const char* data_begin, const char* data_end,
                        fei::FillableMat& mat,
                        bool clear_mat_on_entry = true,
                        bool overwrite_entries = true);

/** return-value is true if the unpacked matrix contains nothing but zeros...*/
bool unpack_CSRMat(const char* buffer_begin, const char* buffer_end, fei::CSRMat& mat);

size_t num_bytes_indices_coefs(const std::vector<int>& indices,
                        const std::vector<double>& coefs);

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

