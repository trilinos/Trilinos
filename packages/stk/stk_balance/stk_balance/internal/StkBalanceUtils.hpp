// Copyright (c) 2013, Sandia Corporation.
 // Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 // the U.S. Government retains certain rights in this software.
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
 //     * Neither the name of Sandia Corporation nor the names of its
 //       contributors may be used to endorse or promote products derived
 //       from this software without specific prior written permission.
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
#ifndef STK_BALANCE_UTILITIES_HPP
#define STK_BALANCE_UTILITIES_HPP

#include <string>
#include <vector>

// stk search (order matters!)
#include <stk_search/IdentProc.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class MetaData; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { namespace mesh { class Entity; } }
namespace stk { namespace mesh { class FieldBase; } }

namespace stk { namespace balance { class BalanceSettings; } }

namespace stk { namespace balance { namespace internal {

typedef stk::search::Box<float> StkBox;
typedef stk::search::IdentProc<stk::mesh::EntityId, int> StkMeshIdent;
typedef std::pair<StkBox, StkMeshIdent> BoxWithStkId;
typedef std::vector< BoxWithStkId > BoxVectorWithStkId;
typedef std::vector<std::pair<StkMeshIdent, StkMeshIdent> > StkSearchResults;

std::string get_parallel_filename(int subdomainIndex, int numSubdomains, const std::string& baseFilename);

int getNumSharedNodesBetweenElements(const ::stk::mesh::BulkData& stkMeshBulkData,
                                     const ::stk::mesh::Entity element1,
                                     const ::stk::mesh::Entity element2);

StkSearchResults getSearchResultsForFacesParticles(stk::mesh::BulkData& stkMeshBulkData, const BalanceSettings &balanceSettings, const stk::mesh::Selector& searchSelector);

void addBoxForNodes(stk::mesh::BulkData& stkMeshBulkData,
                    unsigned numNodes,
                    const stk::mesh::Entity* nodes,
                    const stk::mesh::FieldBase* coord,
                    const double eps,
                    stk::mesh::EntityId elementId,
                    BoxVectorWithStkId& faceBoxes);


void fillFaceBoxesWithIds(stk::mesh::BulkData &stkMeshBulkData, const BalanceSettings & balanceSettings, const stk::mesh::FieldBase* coord, BoxVectorWithStkId &faceBoxes, const stk::mesh::Selector& searchSelector);

const stk::mesh::FieldBase * get_coordinate_field(const stk::mesh::MetaData& meta_data, const std::string& coordinateFieldName);

}}}

#endif
