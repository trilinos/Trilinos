// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PDEOPT_BATHCMANAGER_HPP
#define PDEOPT_BATHCMANAGER_HPP

#include "ROL_TeuchosBatchManager.hpp"
#include "pdevector.hpp"

template<class Real,
         class LO=Tpetra::Map<>::local_ordinal_type,
         class GO=Tpetra::Map<>::global_ordinal_type,
         class Node=Tpetra::Map<>::node_type> 
class PDE_OptVector_BatchManager : public ROL::TeuchosBatchManager<Real,GO> {
private:
  typedef Tpetra::MultiVector<Real,LO,GO,Node> FieldVector;
  typedef std::vector<Real>                    ParamVector;
  typedef PDE_OptVector<Real,LO,GO,Node>       OptVector;

public:
  PDE_OptVector_BatchManager(const ROL::Ptr<const Teuchos::Comm<int> > &comm)
    : ROL::TeuchosBatchManager<Real,GO>(comm) {}

  using ROL::TeuchosBatchManager<Real,GO>::sumAll;
  void sumAll(ROL::Vector<Real> &input, ROL::Vector<Real> &output) {
    // Sum all field components across processors
    ROL::Ptr<ROL::TpetraMultiVector<Real,LO,GO,Node> > input_field_ptr
      = dynamic_cast<OptVector&>(input).getField();
    ROL::Ptr<ROL::TpetraMultiVector<Real,LO,GO,Node> > output_field_ptr
      = dynamic_cast<OptVector&>(output).getField();

    if ( input_field_ptr != ROL::nullPtr ) {
      ROL::Ptr<FieldVector> input_field  = input_field_ptr->getVector();
      ROL::Ptr<FieldVector> output_field = output_field_ptr->getVector();
      size_t input_length  = input_field->getLocalLength();
      size_t output_length = output_field->getLocalLength();
      TEUCHOS_TEST_FOR_EXCEPTION(input_length != output_length, std::invalid_argument,
        ">>> (PDE_OptVector_BatchManager::sumAll): Field dimension mismatch!");

      size_t input_nvec  = input_field->getNumVectors();
      size_t output_nvec = output_field->getNumVectors();
      TEUCHOS_TEST_FOR_EXCEPTION(input_nvec != output_nvec, std::invalid_argument,
        ">>> (PDE_OptVector_BatchManager::sumAll): Field dimension mismatch!");

      for (size_t i = 0; i < input_nvec; ++i) {
        ROL::TeuchosBatchManager<Real,GO>::sumAll((input_field->getDataNonConst(i)).getRawPtr(),
                                                  (output_field->getDataNonConst(i)).getRawPtr(),
                                                  input_length);
      }
    }
    // Sum all parameter components across processors
    ROL::Ptr<ROL::StdVector<Real> > input_param_ptr
      = dynamic_cast<OptVector&>(input).getParameter();
    ROL::Ptr<ROL::StdVector<Real> > output_param_ptr
      = dynamic_cast<OptVector&>(output).getParameter();

    if ( input_param_ptr != ROL::nullPtr ) {
      ROL::Ptr<ParamVector> input_param  = input_param_ptr->getVector();
      ROL::Ptr<ParamVector> output_param = output_param_ptr->getVector();
      size_t input_size  = static_cast<size_t>(input_param->size());
      size_t output_size = static_cast<size_t>(output_param->size());
      TEUCHOS_TEST_FOR_EXCEPTION(input_size != output_size, std::invalid_argument,
        ">>> (PDE_OptVector_BatchManager::SumAll): Parameter dimension mismatch!");

      ROL::TeuchosBatchManager<Real,GO>::sumAll(&input_param->front(),
                                                &output_param->front(),
                                                input_size);
    }
  }
};

#endif
