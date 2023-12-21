// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#include "Tempus_ExplicitTemplateInstantiation.hpp"

#ifdef HAVE_TEMPUS_EXPLICIT_INSTANTIATION
#include "CDR_Model.hpp"
#include "CDR_Model_impl.hpp"
#include "CDR_Model_Tpetra.hpp"
#include "CDR_Model_Tpetra_impl.hpp"


namespace Tempus_Test {
  // Get default Tpetra template types
  using SC = Tpetra::Vector<>::scalar_type;
  using LO = Tpetra::Vector<>::local_ordinal_type;
  using GO = Tpetra::Vector<>::global_ordinal_type;
  using Node = Tpetra::Vector<>::node_type;

  TEMPUS_INSTANTIATE_TEMPLATE_CLASS(CDR_Model)
  TEMPUS_INSTANTIATE_TEMPLATE_CLASS_TPETRA(CDR_Model_Tpetra, SC, LO, GO, Node)
}

#endif // HAVE_TEMPUS_EXPLICIT_INSTANTIATION
