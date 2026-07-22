// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef BasisTable_hpp
#define BasisTable_hpp

#include "Teuchos_RCP.hpp"
#include <Shards_BasicTopologies.hpp>
#include <Shards_CellTopologyData.h>
#include <Shards_CellTopology.hpp>
#include <Intrepid2_Basis.hpp>
#include <percept/function/MDArray.hpp>

  namespace percept {
    class BasisTable
    {
    public:
      using BasisType = Intrepid2::Basis<Kokkos::HostSpace, double, double >;
      using BasisTypeRCP =  Intrepid2::BasisPtr<Kokkos::HostSpace, double, double >;
      using BasisTableMap = std::map<unsigned, BasisTypeRCP >;

      BasisTypeRCP getBasis(shards::CellTopology& topo);

      static BasisTable* getInstance()
      {
        if(instance==nullptr){
          instance = new BasisTable();
        }
        return instance;
      }

    private:
      BasisTable() {}
      ~BasisTable() {m_basisTable.clear();}
      void setupBasisTable();
      static BasisTable * instance;
      BasisTableMap m_basisTable;
    };

  }

#endif
