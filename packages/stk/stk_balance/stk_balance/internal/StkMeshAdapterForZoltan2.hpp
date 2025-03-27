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
#ifndef _STKMESHADAPTERFORZOLTAN2_HPP_
#define _STKMESHADAPTERFORZOLTAN2_HPP_

#include <stddef.h>                     // for size_t, NULL
#include <sys/types.h>                  // for ssize_t
#include <Zoltan2_MeshAdapter.hpp>      // for MeshEntityType, etc
#include <ostream>                      // for basic_ostream, char_traits, etc
#include <vector>                       // for vector
#include "Zoltan2_Adapter.hpp"
#include "Zoltan2_InputTraits.hpp"      // for BasicUserTypes
#include "balanceTypes.hpp"
#include "Zoltan2ParallelGraph.hpp"

typedef Zoltan2::BasicUserTypes<double, BalanceLocalNumber, BalanceGlobalNumber> stkdata_t;

class StkMeshZoltanAdapter : public Zoltan2::MeshAdapter<stkdata_t>
{
public:
  size_t getGlobalNumOf(Zoltan2::MeshEntityType etype) const;

  typedef Zoltan2::MeshAdapter<stkdata_t> base_adapter_t;

  StkMeshZoltanAdapter(const Zoltan2ParallelGraph &graph);

  virtual ~StkMeshZoltanAdapter() { }

  virtual size_t getLocalNumOf(Zoltan2::MeshEntityType etype) const override;
  virtual void getIDsViewOf(Zoltan2::MeshEntityType etype, BalanceGlobalNumber const *&Ids) const override;

  virtual int getDimension() const override;

  virtual void getCoordinatesViewOf(Zoltan2::MeshEntityType etype, const scalar_t *&coords, int &stride, int coordDim) const override;

  virtual int getNumWeightsPerOf(Zoltan2::MeshEntityType etype) const override;

  virtual void getWeightsViewOf(Zoltan2::MeshEntityType etype, const scalar_t *&weights, int &stride, int idx = 0) const override;

  virtual bool avail2ndAdjs(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const override;

  virtual size_t getLocalNum2ndAdjs(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const override;

  virtual void get2ndAdjsView(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through, const BalanceLocalNumber *&offsets, const BalanceGlobalNumber *&adjacencyIds) const override;

  virtual int getNumWeightsPer2ndAdj(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through) const override;

  virtual void get2ndAdjWeightsView(Zoltan2::MeshEntityType sourcetarget, Zoltan2::MeshEntityType through, const scalar_t *&weights, int &stride, int idx) const override;

  void debuggingInfo(int proc_id,  std::ofstream& out) const;

private:
  const Zoltan2ParallelGraph &mGraph;

public: // defaultish
  virtual bool availAdjs(Zoltan2::MeshEntityType /*source*/, Zoltan2::MeshEntityType /*target*/) const override
  {
    return false;
  }

  virtual size_t getLocalNumAdjs(Zoltan2::MeshEntityType /*source*/, Zoltan2::MeshEntityType /*target*/) const override
  {
    return 0;
  }

  virtual void getAdjsView(Zoltan2::MeshEntityType /*source*/, Zoltan2::MeshEntityType /*target*/, const BalanceLocalNumber *&offsets, const BalanceGlobalNumber *& adjacencyIds) const override
  {
    offsets = NULL;
    adjacencyIds = NULL;
    Z2_THROW_NOT_IMPLEMENTED
  }

  virtual bool useDegreeAsWeightOf(Zoltan2::MeshEntityType /*etype*/, int /*idx*/) const override
  {
    return false;
  }
};

#endif
