// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_KOKKOSTUNINGINTERFACE_HPP
#define MUELU_KOKKOSTUNINGINTERFACE_HPP

#include <string>
#include <vector>
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "MueLu_BaseClass.hpp"

namespace MueLu {

/*! @class
  Manages the interface to the KokkosTuning interface

  Note only proc 0 (as defined by comm) will actually have kokkos-tuning instantiated.  The options determined
  by kokkostune will be broadcast to all processors
*/
class KokkosTuningInterface : public BaseClass {
 public:
  KokkosTuningInterface(Teuchos::RCP<const Teuchos::Comm<int> >& comm)
    : comm_(comm) {}

  //  KokkosTuningInterface(Teuchos::RCP<const Teuchos::Comm<int> >& comm, Teuchos::ParameterList& inParams)
  //    : comm_(comm)
  //    , params_(inParams){};

  virtual ~KokkosTuningInterface() { }

  Teuchos::RCP<const Teuchos::ParameterList> GetValidParameterList() const;

  // Sets the input parameters for the KokkosTuneInterface
  // NOTE: These are *not* the parameters which KokkosTuning is going to overwrite, rather the
  // parameters to tell KokkosTuning what to do.
  void SetParameterList(Teuchos::ParameterList& inParams) { params_ = inParams; }

  // Calls Kokkos Tuning to set MueLu Parameters
  void SetMueLuParameters(size_t kokkos_context_id, Teuchos::ParameterList& mueluParams, bool overwrite = true) const;

 private:
  // Utility functions
  void UnpackMueLuMapping();

  // Sets up Kokkos Tuning - This gets called from SetParameterList
  void Setup();

  // Cached data
  Teuchos::RCP<const Teuchos::Comm<int> > comm_;
  mutable Teuchos::ParameterList params_;  // The mutable is a hack to deal with issues in Teuchos


  // Tuning information
  std::vector<Kokkos::Tools::Experimental::VariableValue> in_variables;
  mutable std::vector<Kokkos::Tools::Experimental::VariableValue> out_variables;
  std::vector<std::string> out_names;
  std::vector<std::string> out_typenames;

  // FIXME: Dealing with Kokkos grabbing the pointer
  std::vector<std::vector<int64_t> > int64_ranges;
  std::vector<std::vector<double> > double_ranges;


};

}  // namespace MueLu

#endif  // MUELU_KOKKOSTUNEINTERFACE_HPP
