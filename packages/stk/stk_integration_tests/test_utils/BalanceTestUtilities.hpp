// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
 // Redistribution and use in source and binary forms, with or without
 // modification, are permitted provided that the following conditions are
 // met:
 // 
 //     * Redistributions of source code must retain the above copyright
 //       notice, this list of conditions and the following disclaimer.
 // 
 //     * Redistributions in binary form must reproduce the above
 //       copyright notice, this list of conditions and the following
 //       disclaimer in the documentation and/or other materials provided
 //       with the distribution.
 // 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
 // THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 // "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 // LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 // A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 // OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 // SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 // LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 // DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 // THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 // (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 // OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef BALANCETESTUTILITIES_HPP_
#define BALANCETESTUTILITIES_HPP_

#include <vector>
#include <string>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/SortAndUnique.hpp>
#include <test_utils/NemesisInfo.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class FieldBase; } }

namespace balance_utils
{
void putFieldDataOnMesh(stk::mesh::BulkData &bulkData, const std::vector<int>& coloring);
int getWidth(int numProcsDecomp);
void clearFiles(const std::string &baseFilename, int numProcs);
std::string getFilename(const std::string& filename, int numProcsDecomp, int subdomainId);
void putVectorDataIntoField(stk::mesh::BulkData &bulkData, const stk::mesh::EntityProcVec& vector_of_data, stk::mesh::FieldBase *field);
void putDecompFieldDataOnMesh(stk::mesh::BulkData &bulkData, int proc_rank);
void putEntityProcOnMeshField(stk::mesh::BulkData &bulkData, const stk::mesh::EntityProcVec& coloring);
}


#endif /* BALANCETESTUTILITIES_HPP_ */
