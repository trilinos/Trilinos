// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_BasisDescriptor.hpp"

#include "Phalanx_KokkosDeviceTypes.hpp"

#include "Panzer_HashUtils.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_PointGenerator.hpp"   // includes Kokkos::DynRankView

namespace panzer
{

// Anonymous namespace that holds the coordinate generator,
// this hides any details about the generator from an external
// user and hopefully hides Kokkos from any file that doesn't need it.
namespace { 

/** A generator that builds basis values
  */
class BasisCoordsGenerator :  
    public PointGenerator {
public:
  BasisCoordsGenerator(const int basis_order, const std::string & basis_type)
    : _basis_type(basis_type), _basis_order(basis_order) {}
 
  virtual ~BasisCoordsGenerator() = default;

  virtual Kokkos::DynRankView<double> getPoints(const shards::CellTopology & topo) const override
  {
    Teuchos::RCP<Intrepid2::Basis<PHX::Device::execution_space,double,double> > 
        intrepid_basis = createIntrepid2Basis<PHX::Device::execution_space,double,double>(_basis_type,_basis_order,topo);
    
    Kokkos::DynRankView<double> view(_basis_type+"_ref_coords",intrepid_basis->getCardinality(),topo.getDimension());
  
    intrepid_basis->getDofCoords(view);
  
    return view;
  }

  virtual int numPoints(const shards::CellTopology & topo) const override
  {
    Teuchos::RCP<Intrepid2::Basis<PHX::Device::execution_space,double,double> > 
        intrepid_basis = createIntrepid2Basis<PHX::Device::execution_space,double,double>(_basis_type,_basis_order,topo);
    return intrepid_basis->getCardinality();
  }

  virtual bool hasPoints(const shards::CellTopology & topo) const override
  {
    return true;
  }

protected:
  std::string _basis_type;
  int _basis_order;

private:
  // hidden
  BasisCoordsGenerator();
  BasisCoordsGenerator(const BasisCoordsGenerator &);
};

} // end namespace <anonymous>


BasisDescriptor::BasisDescriptor():
  _basis_type("none"),
  _basis_order(-1)
{
  _key = std::hash<BasisDescriptor>{}(*this);
}

BasisDescriptor::BasisDescriptor(const int basis_order, const std::string & basis_type):
  _basis_type(basis_type),
  _basis_order(basis_order)
{
  _key = std::hash<BasisDescriptor>{}(*this);
}

PointDescriptor 
BasisDescriptor::
getPointDescriptor() const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  std::stringstream ss;
  ss << "cell_ref:" << _basis_type << "-" << _basis_order;

  RCP<PointGenerator> generator = rcp(new BasisCoordsGenerator(_basis_order,_basis_type));

  PointDescriptor pd(ss.str(),generator);

  return pd;
}

bool operator==(const panzer::BasisDescriptor& left,
                const panzer::BasisDescriptor& right)
{
  return ( (left.getType() == right.getType()) &&
           (left.getOrder() == right.getOrder()) );
}

} // end namespace panzer

std::size_t
std::hash<panzer::BasisDescriptor>::operator()(const panzer::BasisDescriptor& desc) const
{
  std::size_t seed = 0;

  panzer::hash_combine(seed,desc.getType());
  panzer::hash_combine(seed,desc.getOrder());

  return seed;
}

