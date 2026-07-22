// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_DOF_CURL_IMPL_HPP
#define PANZER_DOF_CURL_IMPL_HPP

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Phalanx_KokkosDeviceTypes.hpp"

namespace panzer {

namespace {

//**********************************************************************
template <typename ScalarT,typename Array,int spaceDim>
class EvaluateCurlWithSens_Vector {
  PHX::MDField<const ScalarT,Cell,Point> dof_value;
  PHX::MDField<ScalarT,Cell,Point,Dim> dof_curl;
  Array curl_basis;

  int numFields;
  int numPoints;

public:
  typedef typename PHX::Device execution_space;

  EvaluateCurlWithSens_Vector(PHX::MDField<const ScalarT,Cell,Point> in_dof_value,
                              PHX::MDField<ScalarT,Cell,Point,Dim> in_dof_curl,
                              Array in_curl_basis) 
    : dof_value(in_dof_value), dof_curl(in_dof_curl), curl_basis(in_curl_basis)
  {
    numFields = curl_basis.extent(1);
    numPoints = curl_basis.extent(2);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {
    for (int pt=0; pt<numPoints; pt++) {
      for (int d=0; d<spaceDim; d++) {
        // first initialize to the right thing (prevents over writing with 0)
        // then loop over one less basis function
        dof_curl(cell,pt,d) = dof_value(cell, 0) * curl_basis(cell, 0, pt, d);
        for (int bf=1; bf<numFields; bf++)
          dof_curl(cell,pt,d) += dof_value(cell, bf) * curl_basis(cell, bf, pt, d);
      }
    }
  }
};

template <typename ScalarT,typename ArrayT>
void evaluateCurl_withSens_vector(int numCells,
                           PHX::MDField<ScalarT,Cell,Point,Dim> & dof_curl, 
                           PHX::MDField<const ScalarT,Cell,Point> & dof_value,
                           const ArrayT & curl_basis)
{ 
  if(numCells>0) {
    // evaluate at quadrature points
    int numFields = curl_basis.extent(1);
    int numPoints = curl_basis.extent(2);
    int spaceDim  = curl_basis.extent(3);

    for (int cell=0; cell<numCells; cell++) {
      for (int pt=0; pt<numPoints; pt++) {
        for (int d=0; d<spaceDim; d++) {
          // first initialize to the right thing (prevents over writing with 0)
          // then loop over one less basis function
          dof_curl(cell,pt,d) = dof_value(cell, 0) * curl_basis(cell, 0, pt, d);
          for (int bf=1; bf<numFields; bf++)
            dof_curl(cell,pt,d) += dof_value(cell, bf) * curl_basis(cell, bf, pt, d);
        }
      }
    }
  }
}

//**********************************************************************
template <typename ScalarT,typename Array>
class EvaluateCurlWithSens_Scalar {
  PHX::MDField<const ScalarT,Cell,Point> dof_value;
  PHX::MDField<ScalarT,Cell,Point> dof_curl;
  Array curl_basis;

  int numFields;
  int numPoints;

public:
  typedef typename PHX::Device execution_space;

  EvaluateCurlWithSens_Scalar(PHX::MDField<const ScalarT,Cell,Point> in_dof_value,
                              PHX::MDField<ScalarT,Cell,Point> in_dof_curl,
                              Array in_curl_basis) 
    : dof_value(in_dof_value), dof_curl(in_dof_curl), curl_basis(in_curl_basis)
  {
    numFields = curl_basis.extent(1);
    numPoints = curl_basis.extent(2);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {
    for (int pt=0; pt<numPoints; pt++) {
      // first initialize to the right thing (prevents over writing with 0)
      // then loop over one less basis function
      dof_curl(cell,pt) = dof_value(cell, 0) * curl_basis(cell, 0, pt);
      for (int bf=1; bf<numFields; bf++)
        dof_curl(cell,pt) += dof_value(cell, bf) * curl_basis(cell, bf, pt);
    }
  }
};

template <typename ScalarT,typename ArrayT>
void evaluateCurl_withSens_scalar(int numCells,
                           PHX::MDField<ScalarT,Cell,Point> & dof_curl, 
                           PHX::MDField<const ScalarT,Cell,Point> & dof_value,
                           const ArrayT & curl_basis)
{ 
  if(numCells>0) {
    // evaluate at quadrature points
    int numFields = curl_basis.extent(1);
    int numPoints = curl_basis.extent(2);

    for (int cell=0; cell<numCells; cell++) {
      for (int pt=0; pt<numPoints; pt++) {
        // first initialize to the right thing (prevents over writing with 0)
        // then loop over one less basis function
        dof_curl(cell,pt) = dof_value(cell, 0) * curl_basis(cell, 0, pt);
        for (int bf=1; bf<numFields; bf++)
          dof_curl(cell,pt) += dof_value(cell, bf) * curl_basis(cell, bf, pt);
      }
    }
  }
}

//**********************************************************************
template <typename ScalarT,typename Array,int spaceDim>
class EvaluateCurlFastSens_Vector {
  PHX::MDField<const ScalarT,Cell,Point> dof_value;
  PHX::MDField<ScalarT,Cell,Point,Dim> dof_curl;
  PHX::View<const int*> offsets;
  Array curl_basis;

  int numFields;
  int numPoints;

public:
  typedef typename PHX::Device execution_space;

  EvaluateCurlFastSens_Vector(PHX::MDField<const ScalarT,Cell,Point> in_dof_value,
                              PHX::MDField<ScalarT,Cell,Point,Dim> in_dof_curl,
                              PHX::View<const int*> in_offsets,
                              Array in_curl_basis) 
    : dof_value(in_dof_value), dof_curl(in_dof_curl), offsets(in_offsets), curl_basis(in_curl_basis)
  {
    numFields = curl_basis.extent(1);
    numPoints = curl_basis.extent(2);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {
    for (int pt=0; pt<numPoints; pt++) {
      for (int d=0; d<spaceDim; d++) {
        // first initialize to the right thing (prevents over writing with 0)
        // then loop over one less basis function
        dof_curl(cell,pt,d) = dof_value(cell, 0).val() * curl_basis(cell, 0, pt, d);
        dof_curl(cell,pt,d).fastAccessDx(offsets(0)) = dof_value(cell, 0).fastAccessDx(offsets(0)) * curl_basis(cell, 0, pt, d);
        for (int bf=1; bf<numFields; bf++) {
          dof_curl(cell,pt,d).val() += dof_value(cell, bf).val() * curl_basis(cell, bf, pt, d);
          dof_curl(cell,pt,d).fastAccessDx(offsets(bf)) += dof_value(cell, bf).fastAccessDx(offsets(bf)) * curl_basis(cell, bf, pt, d);
        }
      }
    }
  }
};
template <typename ScalarT,typename ArrayT>
void evaluateCurl_fastSens_vector(int numCells,
                           PHX::MDField<ScalarT,Cell,Point,Dim> & dof_curl, 
                           PHX::MDField<const ScalarT,Cell,Point> & dof_value,
                           const std::vector<int> & offsets,
                           const ArrayT & curl_basis)
{ 
  if(numCells>0) {
    int numFields = curl_basis.extent(1);
    int numPoints = curl_basis.extent(2);
    int spaceDim  = curl_basis.extent(3);

    for (int cell=0; cell<numCells; cell++) {
      for (int pt=0; pt<numPoints; pt++) {
        for (int d=0; d<spaceDim; d++) {
          // first initialize to the right thing (prevents over writing with 0)
          // then loop over one less basis function
          dof_curl(cell,pt,d) = ScalarT(numFields, dof_value(cell, 0).val() * curl_basis(cell, 0, pt, d));
          dof_curl(cell,pt,d).fastAccessDx(offsets[0]) = dof_value(cell, 0).fastAccessDx(offsets[0]) * curl_basis(cell, 0, pt, d);
          for (int bf=1; bf<numFields; bf++) {
            dof_curl(cell,pt,d).val() += dof_value(cell, bf).val() * curl_basis(cell, bf, pt, d);
            dof_curl(cell,pt,d).fastAccessDx(offsets[bf]) += dof_value(cell, bf).fastAccessDx(offsets[bf]) * curl_basis(cell, bf, pt, d);
          }
        }
      }
    }
  }
}

//**********************************************************************
template <typename ScalarT,typename Array>
class EvaluateCurlFastSens_Scalar {
  PHX::MDField<const ScalarT,Cell,Point> dof_value;
  PHX::MDField<ScalarT,Cell,Point> dof_curl;
  PHX::View<const int*> offsets;
  Array curl_basis;

  int numFields;
  int numPoints;

public:
  typedef typename PHX::Device execution_space;

  EvaluateCurlFastSens_Scalar(PHX::MDField<const ScalarT,Cell,Point> in_dof_value,
                              PHX::MDField<ScalarT,Cell,Point> in_dof_curl,
                              PHX::View<const int*> in_offsets,
                              Array in_curl_basis) 
    : dof_value(in_dof_value), dof_curl(in_dof_curl), offsets(in_offsets), curl_basis(in_curl_basis)
  {
    numFields = curl_basis.extent(1);
    numPoints = curl_basis.extent(2);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned int cell) const
  {
    for (int pt=0; pt<numPoints; pt++) {
      // first initialize to the right thing (prevents over writing with 0)
      // then loop over one less basis function
      dof_curl(cell,pt) = dof_value(cell, 0).val() * curl_basis(cell, 0, pt);
      dof_curl(cell,pt).fastAccessDx(offsets(0)) = dof_value(cell, 0).fastAccessDx(offsets(0)) * curl_basis(cell, 0, pt);
      for (int bf=1; bf<numFields; bf++) {
        dof_curl(cell,pt).val() += dof_value(cell, bf).val() * curl_basis(cell, bf, pt);
        dof_curl(cell,pt).fastAccessDx(offsets(bf)) += dof_value(cell, bf).fastAccessDx(offsets(bf)) * curl_basis(cell, bf, pt);
      }
    }
  }
};
template <typename ScalarT,typename ArrayT>
void evaluateCurl_fastSens_scalar(int numCells,
                           PHX::MDField<ScalarT,Cell,Point> & dof_curl, 
                           PHX::MDField<const ScalarT,Cell,Point> & dof_value,
                           const std::vector<int> & offsets,
                           const ArrayT & curl_basis)
{ 
  if(numCells>0) {
    int numFields = curl_basis.extent(1);
    int numPoints = curl_basis.extent(2);

    for (int cell=0; cell<numCells; cell++) {
      for (int pt=0; pt<numPoints; pt++) {
        // first initialize to the right thing (prevents over writing with 0)
        // then loop over one less basis function
        dof_curl(cell,pt) = ScalarT(numFields, dof_value(cell, 0).val() * curl_basis(cell, 0, pt));
        dof_curl(cell,pt).fastAccessDx(offsets[0]) = dof_value(cell, 0).fastAccessDx(offsets[0]) * curl_basis(cell, 0, pt);
        for (int bf=1; bf<numFields; bf++) {
          dof_curl(cell,pt).val() += dof_value(cell, bf).val() * curl_basis(cell, bf, pt);
          dof_curl(cell,pt).fastAccessDx(offsets[bf]) += dof_value(cell, bf).fastAccessDx(offsets[bf]) * curl_basis(cell, bf, pt);
        }
      }
    }
  }
}

//**********************************************************************

}

//**********************************************************************
// MOST EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
DOFCurl<EvalT, TRAITS>::
DOFCurl(const Teuchos::ParameterList & p) :
  use_descriptors_(false),
  dof_value( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the curl operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsCurl(),std::logic_error,
                             "DOFCurl: Basis of type \"" << basis->name() << "\" does not support CURL");
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->requiresOrientations(),std::logic_error,
                             "DOFCurl: Basis of type \"" << basis->name() << "\" in DOF Curl should require orientations. So we are throwing.");

  // build dof_curl
  basis_dimension = basis->dimension();
  if(basis_dimension==2) {
     dof_curl_scalar = PHX::MDField<ScalarT,Cell,Point>(p.get<std::string>("Curl Name"), 
      	                              p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar );
     this->addEvaluatedField(dof_curl_scalar);
  }
  else if(basis_dimension==3) {
     dof_curl_vector = PHX::MDField<ScalarT,Cell,Point,Dim>(p.get<std::string>("Curl Name"), 
      	                              p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector );
     this->addEvaluatedField(dof_curl_vector);
  }
  else
  { TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFCurl only works for 2D and 3D basis functions"); } 

  // add to evaluation graph
  this->addDependentField(dof_value);
  
  std::string n = "DOFCurl: " + (basis_dimension==2 ? dof_curl_scalar.fieldTag().name() : dof_curl_vector.fieldTag().name())+ " ()";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
DOFCurl<EvalT, TRAITS>::
DOFCurl(const PHX::FieldTag & input,
        const PHX::FieldTag & output,
        const panzer::BasisDescriptor & bd,
        const panzer::IntegrationDescriptor & id,
        int basis_dim)
  : use_descriptors_(true)
  , bd_(bd) 
  , id_(id) 
  , dof_value(input)
{
  TEUCHOS_ASSERT(bd_.getType()=="HCurl");

  basis_dimension = basis_dim; // user specified

  // build dof_curl
  if(basis_dimension==2) {
     dof_curl_scalar = output;
     this->addEvaluatedField(dof_curl_scalar);
  }
  else if(basis_dimension==3) {
     dof_curl_vector = output;
     this->addEvaluatedField(dof_curl_vector);
  }
  else
  { TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFCurl only works for 2D and 3D basis functions"); } 

  // add to evaluation graph
  this->addDependentField(dof_value);
  
  std::string n = "DOFCurl: " + (basis_dimension==2 ? dof_curl_scalar.fieldTag().name() : dof_curl_vector.fieldTag().name())+ " ()";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
void DOFCurl<EvalT, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(dof_value,fm);
  if(basis_dimension==3)
    this->utils.setFieldData(dof_curl_vector,fm);
  else
    this->utils.setFieldData(dof_curl_scalar,fm);

  if(not use_descriptors_)
    basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
void DOFCurl<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  const panzer::BasisValues2<double> & basisValues = use_descriptors_ ?  this->wda(workset).getBasisValues(bd_,id_)
                                                                      : *this->wda(workset).bases[basis_index];

  if(basis_dimension==3) {
    EvaluateCurlWithSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,3> functor(dof_value,dof_curl_vector,basisValues.curl_basis_vector);
    Kokkos::parallel_for(workset.num_cells,functor);
  }
  else {
    EvaluateCurlWithSens_Scalar<ScalarT,typename BasisValues2<double>::Array_CellBasisIP> functor(dof_value,dof_curl_scalar,basisValues.curl_basis_scalar);
    Kokkos::parallel_for(workset.num_cells,functor);
  }
}

//**********************************************************************

//**********************************************************************
// JACOBIAN EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename TRAITS>                   
DOFCurl<typename TRAITS::Jacobian, TRAITS>::
DOFCurl(const Teuchos::ParameterList & p) :
  use_descriptors_(false),
  dof_value( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // do you specialize because you know where the basis functions are and can
  // skip a large number of AD calculations?
  if(p.isType<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector")) {
    offsets = *p.get<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector");

    // allocate and copy offsets vector to Kokkos array
    PHX::View<int*> offsets_array_nc("offsets",offsets.size());
    auto offsets_array_nc_h = Kokkos::create_mirror_view(offsets_array_nc);
    for(std::size_t i=0;i<offsets.size();i++)
      offsets_array_nc_h(i) = offsets[i];
    Kokkos::deep_copy(offsets_array_nc, offsets_array_nc_h);
    offsets_array = offsets_array_nc;

    accelerate_jacobian = true;  // short cut for identity matrix
  }
  else
    accelerate_jacobian = false; // don't short cut for identity matrix

  // Verify that this basis supports the curl operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsCurl(),std::logic_error,
                             "DOFCurl: Basis of type \"" << basis->name() << "\" does not support CURL");
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->requiresOrientations(),std::logic_error,
                             "DOFCurl: Basis of type \"" << basis->name() << "\" in DOF Curl should require orientations. So we are throwing.");

  // build dof_curl
  basis_dimension = basis->dimension();
  if(basis_dimension==2) {
     dof_curl_scalar = PHX::MDField<ScalarT,Cell,Point>(p.get<std::string>("Curl Name"), 
      	                              p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar );
     this->addEvaluatedField(dof_curl_scalar);
  }
  else if(basis_dimension==3) {
     dof_curl_vector = PHX::MDField<ScalarT,Cell,Point,Dim>(p.get<std::string>("Curl Name"), 
      	                              p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector );
     this->addEvaluatedField(dof_curl_vector);
  }
  else
  { TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFCurl only works for 2D and 3D basis functions"); } 

  // add to evaluation graph
  this->addDependentField(dof_value);
  
  std::string n = "DOFCurl: " + (basis_dimension==2 ? dof_curl_scalar.fieldTag().name() : dof_curl_vector.fieldTag().name())+ " (Jacobian)";
  this->setName(n);
}

//**********************************************************************
template<typename TRAITS>                   
DOFCurl<typename TRAITS::Jacobian, TRAITS>::
DOFCurl(const PHX::FieldTag & input,
        const PHX::FieldTag & output,
        const panzer::BasisDescriptor & bd,
        const panzer::IntegrationDescriptor & id,
        int basis_dim)
  : use_descriptors_(true)
  , bd_(bd) 
  , id_(id) 
  , dof_value(input)
{
  TEUCHOS_ASSERT(bd_.getType()=="HCurl");

  basis_dimension = basis_dim; // user specified

  accelerate_jacobian = false; // don't short cut for identity matrix

  // build dof_curl
  if(basis_dimension==2) {
     dof_curl_scalar = output;
     this->addEvaluatedField(dof_curl_scalar);
  }
  else if(basis_dimension==3) {
     dof_curl_vector = output;
     this->addEvaluatedField(dof_curl_vector);
  }
  else
  { TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"DOFCurl only works for 2D and 3D basis functions"); } 

  // add to evaluation graph
  this->addDependentField(dof_value);
  
  std::string n = "DOFCurl: " + (basis_dimension==2 ? dof_curl_scalar.fieldTag().name() : dof_curl_vector.fieldTag().name())+ " (Jacobian)";
  this->setName(n);
}

//**********************************************************************
template<typename TRAITS>                   
void DOFCurl<typename TRAITS::Jacobian, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(dof_value,fm);
  if(basis_dimension==3)
    this->utils.setFieldData(dof_curl_vector,fm);
  else
    this->utils.setFieldData(dof_curl_scalar,fm);

  if(not use_descriptors_) 
    basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);
}

template<typename TRAITS>                   
void DOFCurl<typename TRAITS::Jacobian,TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  const panzer::BasisValues2<double> & basisValues = use_descriptors_ ?  this->wda(workset).getBasisValues(bd_,id_)
                                                                      : *this->wda(workset).bases[basis_index];

  if(!accelerate_jacobian) {
    if(basis_dimension==3) {
      using Array=typename BasisValues2<double>::ConstArray_CellBasisIPDim;
      Array curl_basis_vector = use_descriptors_ ? basisValues.getCurlVectorBasis(false) : Array(basisValues.curl_basis_vector);
      EvaluateCurlWithSens_Vector<ScalarT,Array,3> functor(dof_value,dof_curl_vector,curl_basis_vector);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
    else {
      using Array=typename BasisValues2<double>::ConstArray_CellBasisIP;
      Array curl_basis_scalar = use_descriptors_ ? basisValues.getCurl2DVectorBasis(false) : Array(basisValues.curl_basis_scalar);
      EvaluateCurlWithSens_Scalar<ScalarT,Array> functor(dof_value,dof_curl_scalar,curl_basis_scalar);
      Kokkos::parallel_for(workset.num_cells,functor);
    }

    return;
  }
  else {

    if(basis_dimension==3) {
      using Array=typename BasisValues2<double>::ConstArray_CellBasisIPDim;
      Array curl_basis_vector = use_descriptors_ ? basisValues.getCurlVectorBasis(false) : Array(basisValues.curl_basis_vector);
      EvaluateCurlFastSens_Vector<ScalarT,Array,3> functor(dof_value,dof_curl_vector,offsets_array,curl_basis_vector);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
    else {
      using Array=typename BasisValues2<double>::ConstArray_CellBasisIP;
      Array curl_basis_scalar = use_descriptors_ ? basisValues.getCurl2DVectorBasis(false) : Array(basisValues.curl_basis_scalar);
      EvaluateCurlFastSens_Scalar<ScalarT,Array> functor(dof_value,dof_curl_scalar,offsets_array,curl_basis_scalar);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
  }
}

}

#endif
