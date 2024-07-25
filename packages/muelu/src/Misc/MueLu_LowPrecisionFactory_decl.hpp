// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_LOWPRECISIONFACTORY_DECL_HPP
#define MUELU_LOWPRECISIONFACTORY_DECL_HPP

#include <string>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_LowPrecisionFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"

namespace MueLu {

/*!
  @class LowPrecisionFactory class.
  @brief Factory for converting matrices to half precision operators
*/

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType>
class LowPrecisionFactory : public SingleLevelFactoryBase {
#undef MUELU_LOWPRECISIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  LowPrecisionFactory() {}

  //! Destructor.
  virtual ~LowPrecisionFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Converts a matrix to half precision operators and returns it in <tt>currentLevel</tt>.
    */
  void Build(Level& currentLevel) const;

  //@}

};  // class LowPrecisionFactory

#if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_FLOAT)
template <class LocalOrdinal, class GlobalOrdinal, class Node>
class LowPrecisionFactory<double, LocalOrdinal, GlobalOrdinal, Node> : public SingleLevelFactoryBase {
  typedef double Scalar;
#undef MUELU_LOWPRECISIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  LowPrecisionFactory() {}

  //! Destructor.
  virtual ~LowPrecisionFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Converts a matrix to half precision operators and returns it in <tt>currentLevel</tt>.
    */
  void Build(Level& currentLevel) const;

  //@}

};  // class LowPrecisionFactory
#endif

#if defined(HAVE_TPETRA_INST_COMPLEX_DOUBLE) && defined(HAVE_TPETRA_INST_COMPLEX_FLOAT)
template <class LocalOrdinal, class GlobalOrdinal, class Node>
class LowPrecisionFactory<std::complex<double>, LocalOrdinal, GlobalOrdinal, Node> : public SingleLevelFactoryBase {
  typedef std::complex<double> Scalar;
#undef MUELU_LOWPRECISIONFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  //! @name Constructors/Destructors.
  //@{

  LowPrecisionFactory() {}

  //! Destructor.
  virtual ~LowPrecisionFactory() {}

  RCP<const ParameterList> GetValidParameterList() const;

  //@}

  //! Input
  //@{

  void DeclareInput(Level& currentLevel) const;

  //@}

  //! @name Build methods.
  //@{

  /*!
    @brief Build method.

    Converts a matrix to half precision operators and returns it in <tt>currentLevel</tt>.
    */
  void Build(Level& currentLevel) const;

  //@}

};  // class LowPrecisionFactory
#endif

}  // namespace MueLu

#define MUELU_LOWPRECISIONFACTORY_SHORT
#endif  // MUELU_LOWPRECISIONFACTORY_DECL_HPP
