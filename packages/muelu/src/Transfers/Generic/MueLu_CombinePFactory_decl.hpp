// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_COMBINEPFACTORY_DECL_HPP
#define MUELU_COMBINEPFACTORY_DECL_HPP

#include "MueLu_ConfigDefs.hpp"
#include <Teuchos_SerialDenseVector.hpp>

#include <Xpetra_MultiVector.hpp>
#include <Xpetra_Matrix_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_CombinePFactory_fwd.hpp"

#include "MueLu_Level_fwd.hpp"
#include <type_traits>

namespace MueLuTests {
// Forward declaration of friend tester class used to UnitTest CombinePFactory
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class CombinePFactoryTester;
}  // namespace MueLuTests

namespace MueLu {

namespace Details {

/*!
 * \brief Compile-time trait indicating whether
 *        \c CombinePFactory::BuildPBlocked may use its Thyra-based blocked
 *        prolongator construction for a given
 *        \c Scalar/\c LocalOrdinal/\c GlobalOrdinal/\c Node tuple.
 *
 * This trait is used to guard the Thyra-dependent implementation of
 * \c BuildPBlocked against template instantiations for which the required
 * Thyra/Tpetra/Xpetra explicit template instantiations (ETI) are not available.
 * Without this guard, unsupported type combinations may compile successfully
 * but fail later at link time due to missing Thyra-side instantiations.
 *
 * The trait therefore answers a narrower question than "is this scalar type
 * generally supported by Tpetra?"  It is intended to represent whether the
 * specific interop path used by \c BuildPBlocked is available for the given
 * template arguments.
 *
 * The primary template defaults to \c std::false_type.  Known-good type
 * combinations should specialize this trait to \c std::true_type.
 *
 * Advanced users who provide the necessary Thyra ETI externally may also
 * specialize this trait for additional type combinations.  Such a
 * specialization must be visible at the point where
 * \c CombinePFactory<...>::BuildPBlocked is instantiated.
 *
 * \note Setting this trait to \c true does not create the required ETI; it
 *       only informs MueLu that the caller guarantees the necessary Thyra-side
 *       support exists.
 */
template <class Scalar, class LO, class GO, class Node>
struct has_build_p_blocked_thyra_eti : std::false_type {};
}  // namespace Details

/*!
  @class CombinePFactory
  @ingroup MueLuTransferClasses
  @brief Prolongator factory that replicates 'Psubblock' matrix to create new prolongator suitable for PDE systems


  Takes a set of previously generated prolongators each for a subsystem multiphysics PDE and effectively makes a block diagonal prolongator
  by combining the "combo: npdes" times so that it can be used with a PDE system. A normal use case
  would be to run an existing MueLu multigrid algorithm on a set of scalar PDEs to generate a hierarchy.  Then use
  something like

  MueLu::HierarchyUtils<SC,LO,GO,NO>::CopyBetweenHierarchies(*scalarHierarchy0,*systemHierarchy, "P",  "Psubblock0", "RCP<Matrix>");
                                                                             .                                   .
                                                                             .                                   .
                                                                             .                                   .
  MueLu::HierarchyUtils<SC,LO,GO,NO>::CopyBetweenHierarchies(*scalarHierarchy7,*systemHierarchy, "P",  "Psubblock6", "RCP<Matrix>");

  to copy them to a new hierarchy. This new hierarchy would then have <Parameter name="multigrid algorithm" type="string"  value="combo"/>
  in its parameter list to then invoke ComboPFactory. Currently, this is used in src/Operators/MueLu_Multiphysics_def.hpp with an example
  in test/multiphysics.

  The end result is a block diagonal matrix corresponding to
                   [Psubblock0                                                                  ]
                   [           Psubblock1                                                       ]
                   [                      Psubblock2                                            ]
          P   =    [                                 Psubblock3                                 ]
                   [                                            Psubblock4                      ]
                   [                                                       Psubblock5           ]
                   [                                                                  Psubblock6]

  ## Input/output of CombinePFactory ##

  ### User parameters of SemiCoarsenPFactory ###
  | Parameter   | type    | default   | master.xml | validated | requested | description                                                                      |
  | ------------|---------|-----------|:----------:|:---------:|:---------:|----------------------------------------------------------------------------------|
  |combo: npdes | int     |  1        |     *      |     *     |           | Specifies the number of Psubblock# matrices to expect                            |
  | Psubblock0   | Matrix |           |            |           |           | Matrix defining P(0,0) block                                                     |
  | Psubblock1   | Matrix |           |            |           |           | Matrix defining P(1,1) block                                                     |
  |      .       |    .   |           |            |           |           |            .                                                                     |
  |      .       |    .   |           |            |           |           |            .                                                                     |
  |      .       |    .   |           |            |           |           |            .                                                                     |
  | Psubblockk   | Matrix |           |            |           |           | Matrix defining P(k,k) block                                                     |



  The * in the @c master.xml column denotes that the parameter is defined in the @c master.xml file.<br>
  The * in the @c validated column means that the parameter is declared in the list of valid input parameters (see CombineCoarsenPFactory::GetValidParameters).<br>
  The * in the @c requested column states that the data is requested as input with all dependencies (see CombineCoarsenPFactory::DeclareInput).

  ### Variables provided by CombinePFactory ###

  After CombinePFactory::Build the following data is available (if requested)

  | Parameter         | generated by             | description
  |-------------------|--------------------------|------------------------------------------------------------------------------------------------------------------|
  | P                 | CombinePFactory          | Prolongator                                                                                                      |

*/
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class CombinePFactory : public PFactory {
#undef MUELU_COMBINEPFACTORY_SHORT
#include "MueLu_UseShortNames.hpp"

 public:
  friend class MueLuTests::CombinePFactoryTester<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  //! @name Constructors/Destructors.
  //@{

  //! Constructor
  CombinePFactory() {}

  //! Destructor.
  virtual ~CombinePFactory() {}
  //@}

  RCP<const ParameterList> GetValidParameterList() const;

  //! Input
  //@{

  void DeclareInput(Level& fineLevel, Level& coarseLevel) const;

  //@}

  //! @name Build methods.
  //@{

  void Build(Level& fineLevel, Level& coarseLevel) const;
  void BuildP(Level& fineLevel, Level& coarseLevel) const;

  //@}

 private:
  void BuildPBlocked(Level& fineLevel, Level& coarseLevel) const;
  void BuildPBlockedImpl(Level& fineLevel, Level& coarseLevel, std::false_type) const;
#ifdef HAVE_XPETRA_THYRA
  template <class S = Scalar,
            std::enable_if_t<
                MueLu::Details::has_build_p_blocked_thyra_eti<
                    S, LocalOrdinal, GlobalOrdinal, Node>::value,
                int> = 0>
  void BuildPBlockedImpl(Level& fineLevel, Level& coarseLevel, std::true_type) const;
#endif
  int numPDEs_;

};  // class CombinePFactory

}  // namespace MueLu

#define MUELU_COMBINEPFACTORY_SHORT
#endif  // MUELU_COMBINEPFACTORY_DECL_HPP
