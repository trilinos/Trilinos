// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   Panzer_Integrator_CurlBasisDotVector_impl_hpp
#define   Panzer_Integrator_CurlBasisDotVector_impl_hpp

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Intrepid2
#include "Intrepid2_FunctionSpaceTools.hpp"

// Panzer
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"

// Phalanx
#include "Phalanx_KokkosDeviceTypes.hpp"

namespace panzer
{
  /////////////////////////////////////////////////////////////////////////////
  //
  //  Main Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  Integrator_CurlBasisDotVector<EvalT, Traits>::
  Integrator_CurlBasisDotVector(
    const panzer::EvaluatorStyle&   evalStyle,
    const std::string&              resName,
    const std::string&              valName,
    const panzer::BasisIRLayout&    basis,
    const panzer::IntegrationRule&  ir,
    const double&                   multiplier, /* = 1 */
    const std::vector<std::string>& fmNames     /* =
      std::vector<std::string>() */)
    :
    evalStyle_(evalStyle),
    useDescriptors_(false),
    multiplier_(multiplier),
    basisName_(basis.name())
  {
    using PHX::View;
    using panzer::BASIS;
    using panzer::Cell;
    using panzer::Dim;
    using panzer::EvaluatorStyle;
    using panzer::IP;
    using panzer::PureBasis;
    using PHX::MDField;
    using std::invalid_argument;
    using std::logic_error;
    using std::string;
    using Teuchos::RCP;

    // Ensure the input makes sense.
    TEUCHOS_TEST_FOR_EXCEPTION(resName == "", invalid_argument, "Error:  "    \
      "Integrator_CurlBasisDotVector called with an empty residual name.")
    TEUCHOS_TEST_FOR_EXCEPTION(valName == "", invalid_argument, "Error:  "    \
      "Integrator_CurlBasisDotVector called with an empty value name.")
    RCP<const PureBasis> tmpBasis = basis.getBasis();
    TEUCHOS_TEST_FOR_EXCEPTION(not tmpBasis->isVectorBasis(), logic_error,
      "Error:  Integrator_CurlBasisDotVector:  Basis of type \""
      << tmpBasis->name() << "\" is not a vector basis.")
    TEUCHOS_TEST_FOR_EXCEPTION(not tmpBasis->requiresOrientations(),
      logic_error, "Error:  Integrator_CurlBasisDotVector:  Basis of type \""
      << tmpBasis->name() << "\" does not require orientations, though it "   \
      "should for its use in this Evaluator.")
    TEUCHOS_TEST_FOR_EXCEPTION(not tmpBasis->supportsCurl(), logic_error,
      "Error:  Integrator_CurlBasisDotVector:  Basis of type \""
      << tmpBasis->name() << "\" does not support curl.")
    TEUCHOS_TEST_FOR_EXCEPTION(not ((tmpBasis->dimension() == 2) or
      (tmpBasis->dimension() == 3)), logic_error,
      "Error:  Integrator_CurlBasisDotVector requires either a two- or "      \
      "three-dimensional basis.  The basis \"" << tmpBasis->name()
      << "\" is neither.")

    // Use a scalar field only if we're dealing with a two-dimensional case.
    spaceDim_ = tmpBasis->dimension();

    // Create the field for the vector-values quantity we're integrating.
    if (spaceDim_ == 2)
    {
      vector2D_ = MDField<const ScalarT, Cell, IP>(valName, ir.dl_scalar);
      this->addDependentField(vector2D_);
    }
    else // if (spaceDim_ == 3)
    {
      vector3D_ = MDField<const ScalarT, Cell, IP, Dim>(valName, ir.dl_vector);
      this->addDependentField(vector3D_);
    } // end if spaceDim_ is 2 or 3

    // Create the field that we're either contributing to or evaluating
    // (storing).
    field_ = MDField<ScalarT, Cell, BASIS>(resName, basis.functional);
    if (evalStyle == EvaluatorStyle::CONTRIBUTES)
      this->addContributedField(field_);
    else // if (evalStyle == EvaluatorStyle::EVALUATES)
      this->addEvaluatedField(field_);

    // Add the dependent field multipliers, if there are any.
    int i(0);
    fieldMults_.resize(fmNames.size());
    kokkosFieldMults_ = View<Kokkos::View<const ScalarT**, typename PHX::DevLayout<ScalarT>::type, Kokkos::MemoryUnmanaged>*>(
      "CurlBasisDotVector::KokkosFieldMultipliers", fmNames.size());
    for (const auto& name : fmNames)
    {
      fieldMults_[i++] = MDField<const ScalarT, Cell, IP>(name, ir.dl_scalar);
      this->addDependentField(fieldMults_[i - 1]);
    } // end loop over the field multipliers

    // Set the name of this object.
    string n("Integrator_CurlBasisDotVector (");
    if (evalStyle == EvaluatorStyle::CONTRIBUTES)
      n += "CONTRIBUTES";
    else // if (evalStyle == EvaluatorStyle::EVALUATES)
      n += "EVALUATES";
    n += "):  " + field_.fieldTag().name();
    this->setName(n);
  } // end of Main Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  ParameterList Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  Integrator_CurlBasisDotVector<EvalT, Traits>::
  Integrator_CurlBasisDotVector(
    const Teuchos::ParameterList& p)
    :
    Integrator_CurlBasisDotVector(
      panzer::EvaluatorStyle::EVALUATES,
      p.get<std::string>("Residual Name"),
      p.get<std::string>("Value Name"),
      (*p.get<Teuchos::RCP<panzer::BasisIRLayout>>("Basis")),
      (*p.get<Teuchos::RCP<panzer::IntegrationRule>>("IR")),
      p.get<double>("Multiplier"),
      p.isType<Teuchos::RCP<const std::vector<std::string>>>
        ("Field Multipliers") ?
        (*p.get<Teuchos::RCP<const std::vector<std::string>>>
        ("Field Multipliers")) : std::vector<std::string>())
  {
    using Teuchos::ParameterList;
    using Teuchos::RCP;

    // Ensure that the input ParameterList didn't contain any bogus entries.
    RCP<ParameterList> validParams = this->getValidParameters();
    p.validateParameters(*validParams);
  } // end of ParameterList Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  FieldTag Constructor
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename TRAITS>
  Integrator_CurlBasisDotVector<EvalT, TRAITS>::
  Integrator_CurlBasisDotVector(
    const panzer::EvaluatorStyle&        evalStyle,
    const PHX::FieldTag&                 resTag,
    const PHX::FieldTag&                 valTag,
    const panzer::BasisDescriptor&       bd,
    const panzer::IntegrationDescriptor& id,
    const int&                           spaceDim,   /* = 3 */
    const double&                        multiplier, /* = 1 */
    const std::vector<PHX::FieldTag>&    multipliers /* =
      std::vector<PHX::FieldTag>() */)
    :
    evalStyle_(evalStyle),
    useDescriptors_(true),
    bd_(bd),
    id_(id),
    multiplier_(multiplier),
    spaceDim_(spaceDim)
  {
    using PHX::View;
    using panzer::EvaluatorStyle;
    using std::logic_error;
    using std::string;

    // Ensure the input makes sense.
    TEUCHOS_TEST_FOR_EXCEPTION(bd_.getType() != "HCurl", logic_error,
      "Error:  Integrator_CurlBasisDotVector:  Basis of type \""
      << bd_.getType() << "\" does not support curl.")
    TEUCHOS_TEST_FOR_EXCEPTION(not ((spaceDim == 2) or (spaceDim == 3)),
      logic_error, "Error:  Integrator_CurlBasisDotVector works on either " \
      "two- or three-dimensional problems.  You've provided spaceDim = "
      << spaceDim << ".")

    // Create the field for the vector-valued quantity we're integrating.
    if (spaceDim_ == 2)
    {
      vector2D_ = valTag;
      this->addDependentField(vector2D_);
    }
    else // if (spaceDim_ == 3)
    {
      vector3D_ = valTag;
      this->addDependentField(vector3D_);
    } // end if spaceDim_ is 2 or 3

    // Create the field that we're either contributing to or evaluating
    // (storing).
    field_ = resTag;
    if (evalStyle == EvaluatorStyle::CONTRIBUTES)
      this->addContributedField(field_);
    else // if (evalStyle == EvaluatorStyle::EVALUATES)
      this->addEvaluatedField(field_);

    // Add the dependent field multipliers, if there are any.
    int i(0);
    fieldMults_.resize(multipliers.size());
    kokkosFieldMults_ = View<Kokkos::View<const ScalarT**, typename PHX::DevLayout<ScalarT>::type, Kokkos::MemoryUnmanaged>*>(
      "CurlBasisDotVector::KokkosFieldMultipliers", multipliers.size());
    for (const auto& fm : multipliers)
    {
      fieldMults_[i++] = fm;
      this->addDependentField(fieldMults_[i - 1]);
    } // end loop over the field multipliers

    // Set the name of this object.
    string n("Integrator_CurlBasisDotVector (");
    if (evalStyle == EvaluatorStyle::CONTRIBUTES)
      n += "CONTRIBUTES";
    else // if (evalStyle == EvaluatorStyle::EVALUATES)
      n += "EVALUATES";
    n += "):  " + field_.fieldTag().name();
    this->setName(n);
  } // end of Other Constructor

  /////////////////////////////////////////////////////////////////////////////
  //
  //  postRegistrationSetup()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  void
  Integrator_CurlBasisDotVector<EvalT, Traits>::
  postRegistrationSetup(
    typename Traits::SetupData sd,
    PHX::FieldManager<Traits>& fm)
  {
    using panzer::getBasisIndex;
    using PHX::MDField;
    using std::vector;

    // Get the PHX::Views of the field multipliers.
    auto kokkosFieldMults_h = Kokkos::create_mirror_view(kokkosFieldMults_);
    for (size_t i(0); i < fieldMults_.size(); ++i)
      kokkosFieldMults_h(i) = fieldMults_[i].get_static_view();
    Kokkos::deep_copy(kokkosFieldMults_, kokkosFieldMults_h);
    // Determine the index in the Workset bases for our particular basis
    // name.
    if (not useDescriptors_)
      basisIndex_ = getBasisIndex(basisName_, (*sd.worksets_)[0], this->wda);

    // Set up the field that will be used to build of the result of this
    // integration.
    MDFieldArrayFactory af("",
      fm.template getKokkosExtendedDataTypeDimensions<EvalT>(), true);
    if (spaceDim_ == 2)
      result2D_ = af.buildStaticArray<ScalarT, Cell, IP>(
        "Integrator_CurlBasisDotVector:  2-D Result", vector2D_.extent(0),
        vector2D_.extent(1));
    else // if (spaceDim_ == 3)
      result3D_ = af.buildStaticArray<ScalarT, Cell, IP, Dim>(
        "Integrator_CurlBasisDotVector:  3-D Result", vector3D_.extent(0),
        vector3D_.extent(1), 3);
  } // end of postRegistrationSetup()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  Anonymous namespace containing classes for performing the integration.
  //
  /////////////////////////////////////////////////////////////////////////////
  namespace
  {
    /**
     *  \brief Multiply the integrand by the scalar multiplier (\f$ M \f$) and
     *         any field multipliers (\f$ a(x) \f$, \f$ b(x) \f$, etc.) out in
     *         front of the integral.
     *
     *  This must happen before the integration itself is performed by
     *  `Integrate2D`.
     *
     *  \note This class is used when integrating a two-dimensional problem.
     *        Three-dimensional problems use the `PreMultiply3D` class.
     */
    template<typename ScalarT>
    class PreMultiply2D
    {
      public:

        /**
         *  \brief The Kokkos dispatch tag to pre-multiply the integrand by the
         *         scalar multiplier.
         */
        struct ScalarMultiplierTag
        {
        }; // end of struct ScalarMultiplierTag

        /**
         *  \brief The Kokkos dispatch tag to pre-multiply the integrand by a
         *         field multiplier.
         */
        struct FieldMultiplierTag
        {
        }; // end of struct FieldMultiplierTag

        /**
         *  \brief Multiply the integrand by the scalar multiplier (\f$ M \f$)
         *         out in front of the integral.
         *
         *  Loop over the quadrature points and scale the integrand by the
         *  scalar multiplier.
         *
         *  \note This must be called before the `FieldMultiplierTag`
         *        `operator()()` routine, as this initializes the resulting
         *        field.
         */
        KOKKOS_INLINE_FUNCTION
        void operator()(const ScalarMultiplierTag, const unsigned cell) const
        {
          int numQP(result.extent(1));
          for (int qp(0); qp < numQP; ++qp)
            result(cell, qp) = multiplier * vectorField(cell, qp);
        } // end of ScalarMultiplierTag operator()()

        /**
         *  \brief Multiply the integrand by one of the field multipliers
         *         (\f$ a(x) \f$, \f$ b(x) \f$, etc.) out in front of the
         *         integral.
         *
         *  Loop over the quadrature points and scale the integrand by the
         *  field multiplier.
         *
         *  \note This must be called after the `ScalarMultiplierTag`
         *        `operator()()` routine, as that one initializes the resulting
         *        field.
         */
        KOKKOS_INLINE_FUNCTION
        void operator()(const FieldMultiplierTag, const unsigned cell) const
        {
          int numQP(result.extent(1));
          for (int qp(0); qp < numQP; ++qp)
            result(cell, qp) *= fieldMult(cell, qp);
        } // end of FieldMultiplierTag operator()()

        /**
         *  \brief This tells Kokkos to only execute this functor in the
         *         `PHX::Device` execution space.
         */
        using execution_space = PHX::Device;

        /**
         *  \brief A field that will be used to build up the result of the
         *         integral we're performing.
         */
        PHX::MDField<ScalarT, panzer::Cell, panzer::IP> result;

        /**
         *  \brief The scalar multiplier out in front of the integral (\f$ M
         *         \f$).
         */
        double multiplier;

        /**
         *  \brief A field representing the vector-valued function we're
         *         integrating (\f$ \vec{s} \f$).
         */
        PHX::MDField<const ScalarT, panzer::Cell, panzer::IP> vectorField;

        /**
         *  \brief One of the field multipliers (\f$ a(x) \f$, \f$ b(x) \f$,
         *         etc.) out in front of the integral.
         */
        PHX::MDField<const ScalarT, panzer::Cell, panzer::IP> fieldMult;
    }; // end of class PreMultiply2D

    /**
     *  \brief Multiply the integrand by the scalar multiplier (\f$ M \f$) and
     *         any field multipliers (\f$ a(x) \f$, \f$ b(x) \f$, etc.) out in
     *         front of the integral.
     *
     *  This must happen before the integration itself is performed by
     *  `Integrate3D`.
     *
     *  \note This class is used when integrating a three-dimensional problem.
     *        Two-dimensional problems use the `PreMultiply2D` class.
     */
    template<typename ScalarT>
    class PreMultiply3D
    {
      public:

        /**
         *  \brief The Kokkos dispatch tag to pre-multiply the integrand by the
         *         scalar multiplier.
         */
        struct ScalarMultiplierTag
        {
        }; // end of struct ScalarMultiplierTag

        /**
         *  \brief The Kokkos dispatch tag to pre-multiply the integrand by a
         *         field multiplier.
         */
        struct FieldMultiplierTag
        {
        }; // end of struct FieldMultiplierTag

        /**
         *  \brief Multiply the integrand by the scalar multiplier (\f$ M \f$)
         *         out in front of the integral.
         *
         *  Loop over the quadrature points and dimensions of our vector field
         *  and scale the integrand by the scalar multiplier.
         *
         *  \note This must be called before the `FieldMultiplierTag`
         *        `operator()()` routine, as this initializes the resulting
         *        field.
         */
        KOKKOS_INLINE_FUNCTION
        void operator()(const ScalarMultiplierTag, const unsigned cell) const
        {
          int numQP(result.extent(1)), numDim(result.extent(2));
          for (int qp(0); qp < numQP; ++qp)
            for (int dim(0); dim < numDim; ++dim)
              result(cell, qp, dim) = multiplier * vectorField(cell, qp, dim);
        } // end of ScalarMultiplierTag operator()()

        /**
         *  \brief Multiply the integrand by one of the field multipliers
         *         (\f$ a(x) \f$, \f$ b(x) \f$, etc.) out in front of the
         *         integral.
         *
         *  Loop over the quadrature points and dimensions of our vector field
         *  and scale the integrand by the field multiplier.
         *
         *  \note This must be called after the `ScalarMultiplierTag`
         *        `operator()()` routine, as that one initializes the resulting
         *        field.
         */
        KOKKOS_INLINE_FUNCTION
        void operator()(const FieldMultiplierTag, const unsigned cell) const
        {
          int numQP(result.extent(1)), numDim(result.extent(2));
          for (int qp(0); qp < numQP; ++qp)
            for (int dim(0); dim < numDim; ++dim)
              result(cell, qp, dim) *= fieldMult(cell, qp);
        } // end of FieldMultiplierTag operator()()

        /**
         *  \brief This tells Kokkos to only execute this functor in the
         *         `PHX::Device` execution space.
         */
        using execution_space = PHX::Device;

        /**
         *  \brief A field that will be used to build up the result of the
         *         integral we're performing.
         */
        PHX::MDField<ScalarT, panzer::Cell, panzer::IP, panzer::Dim> result;

        /**
         *  \brief The scalar multiplier out in front of the integral (\f$ M
         *         \f$).
         */
        double multiplier;

        /**
         *  \brief A field representing the vector-valued function we're
         *         integrating (\f$ \vec{s} \f$).
         */
        PHX::MDField<const ScalarT, panzer::Cell, panzer::IP, panzer::Dim>
        vectorField;

        /**
         *  \brief One of the field multipliers (\f$ a(x) \f$, \f$ b(x) \f$,
         *         etc.) out in front of the integral.
         */
        PHX::MDField<const ScalarT, panzer::Cell, panzer::IP> fieldMult;
    }; // end of class PreMultiply3D

    /**
     *  \brief Perform the actual integration, looping over the bases,
     *         quadrature points, and dimensions of our vector field.
     *
     *  This must happen after the integrand has already been pre-multiplied
     *  by the scalar multiplier (\f$ M \f$) and any field multipliers
     *  (\f$ a(x) \f$, \f$ b(x) \f$, etc.) out in front of the integral.
     *
     *  \note This class is used when integrating a three-dimensional problem.
     *        Two-dimensional problems use the `Integrate2D` class.
     */
    template<typename ScalarT, int spaceDim>
    class Integrate3D
    {
      public:

        /**
         *  \brief Perform the actual integration.
         *
         *  Loop over the bases, quadrature points, and dimensions in our
         *  vector field, summing up the integrand times the appropriate basis
         *  information.
         *
         *  \note If this is a `EVALUATES` style `Evaluator`, the field will be
         *        initialized to zero.  If it's a `CONTRIBUTES` style one
         *        instead, the field will not be initialized, as that's assumed
         *        to have happened in some other `Evaluator`.
         */
        KOKKOS_INLINE_FUNCTION
        void operator()(const unsigned cell) const
        {
          using panzer::EvaluatorStyle;
          int numBases(weightedCurlBasis.extent(1)),
            numQP(weightedCurlBasis.extent(2));
          for (int basis(0); basis < numBases; ++basis)
          {
            if (evalStyle == EvaluatorStyle::EVALUATES)
              field(cell, basis) = 0.0;
            for (int qp(0); qp < numQP; ++qp)
              for (int dim(0); dim < spaceDim; ++dim)
                field(cell, basis) += result(cell, qp, dim) *
                  weightedCurlBasis(cell, basis, qp, dim);
          } // end loop over the basis functions
        } // end of operator()

        /**
         *  \brief This tells Kokkos to only execute this functor in the
         *         `PHX::Device` execution space.
         */
        using execution_space = PHX::Device;

        /**
         *  \brief A field representing the result of this integration.
         */
        PHX::MDField<ScalarT, panzer::Cell, panzer::IP, panzer::Dim> result;

        /**
         *  \brief A field to which we'll contribute, or in which we'll store,
         *         the result of computing this integral.
         */
        PHX::MDField<ScalarT, panzer::Cell, panzer::BASIS> field;

        /**
         *  \brief The vector basis information necessary for integration.
         */
        PHX::MDField<const double, panzer::Cell, panzer::BASIS, panzer::IP,
          panzer::Dim> weightedCurlBasis;

        /**
         *  \brief The `EvaluatorStyle` of the parent
         *         `Integrator_CurlBasisDotVector` object.
         */
        panzer::EvaluatorStyle evalStyle;
    }; // end of class Integrate3D

    /**
     *  \brief Perform the actual integration, looping over the bases,
     *         quadrature points, and dimensions of our vector field.
     *
     *  This must happen after the integrand has already been pre-multiplied
     *  by the scalar multiplier (\f$ M \f$) and any field multipliers
     *  (\f$ a(x) \f$, \f$ b(x) \f$, etc.) out in front of the integral.
     *
     *  \note This class is used when integrating a two-dimensional problem.
     *        Three-dimensional problems use the `Integrate3D` class.
     */
    template<typename ScalarT>
    class Integrate2D
    {
      public:

        /**
         *  \brief Perform the actual integration.
         *
         *  Loop over the bases, quadrature points, summing up the integrand
         *  times the appropriate basis information.
         *
         *  \note If this is a `EVALUATES` style `Evaluator`, the field will be
         *        initialized to zero.  If it's a `CONTRIBUTES` style one
         *        instead, the field will not be initialized, as that's assumed
         *        to have happened in some other `Evaluator`.
         */
        KOKKOS_INLINE_FUNCTION
        void operator()(const unsigned cell) const
        {
          using panzer::EvaluatorStyle;
          int numBases(weightedCurlBasis.extent(1)),
            numQP(weightedCurlBasis.extent(2));
          for (int basis(0); basis < numBases; ++basis)
          {
            if (evalStyle == EvaluatorStyle::EVALUATES)
              field(cell, basis) = 0.0;
            for (int qp(0); qp < numQP; ++qp)
              field(cell, basis) += result(cell, qp) *
                weightedCurlBasis(cell, basis, qp);
          } // end loop over the basis functions
        } // end of operator()

        /**
         *  \brief This tells Kokkos to only execute this functor in the
         *         `PHX::Device` execution space.
         */
        using execution_space = PHX::Device;

        /**
         *  \brief A field representing the result of this integration.
         */
        PHX::MDField<ScalarT, panzer::Cell, panzer::IP> result;

        /**
         *  \brief A field to which we'll contribute, or in which we'll store,
         *         the result of computing this integral.
         */
        PHX::MDField<ScalarT, panzer::Cell, panzer::BASIS> field;

        /**
         *  \brief The vector basis information necessary for integration.
         */
        PHX::MDField<const double, panzer::Cell, panzer::BASIS, panzer::IP>
        weightedCurlBasis;

        /**
         *  \brief The `EvaluatorStyle` of the parent
         *         `Integrator_CurlBasisDotVector` object.
         */
        panzer::EvaluatorStyle evalStyle;
    }; // end of class Integrate2D

  } // end of anonymous namespace

  /////////////////////////////////////////////////////////////////////////////
  //
  //  evaluateFields()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename Traits>
  void
  Integrator_CurlBasisDotVector<EvalT, Traits>::
  evaluateFields(
    typename Traits::EvalData workset)
  {
    using Kokkos::parallel_for;
    using Kokkos::RangePolicy;
    using panzer::BasisValues2;
    using PHX::Device;
    using PHX::MDField;
    using std::vector;

    // Grab the basis information.
    const BasisValues2<double>& bv = useDescriptors_ ?
      this->wda(workset).getBasisValues(bd_, id_) :
      *this->wda(workset).bases[basisIndex_];

    // If we're dealing with a two- or three-dimensional problem...
    if (spaceDim_ == 2)
    {
      // Create an object to multiply the integrand by the scalar multiplier
      // and any field multipliers out in front of the integral.
      using PreMultiply         = PreMultiply2D<ScalarT>;
      using ScalarMultiplierTag = typename PreMultiply::ScalarMultiplierTag;
      using FieldMultiplierTag  = typename PreMultiply::FieldMultiplierTag;
      PreMultiply preMultiply;
      preMultiply.result      = result2D_;
      preMultiply.multiplier  = multiplier_;
      preMultiply.vectorField = vector2D_;

      // Multiply the integrand by the scalar multiplier out in front of the
      // integral.
      parallel_for(RangePolicy<Device, ScalarMultiplierTag>(0,
        workset.num_cells), preMultiply);

      // Multiply the integrand by any field multipliers out in front of the
      // integral.
      for (const auto& field : fieldMults_)
      {
        preMultiply.fieldMult = field;
        parallel_for(RangePolicy<Device, FieldMultiplierTag>(0,
          workset.num_cells), preMultiply);
      } // end loop over the field multipliers

      // Create an object to do the actual integration and then do it.
      Integrate2D<ScalarT> integrate;
      integrate.result            = result2D_;
      integrate.field             = field_;
      using Array=typename BasisValues2<double>::ConstArray_CellBasisIP;
      integrate.weightedCurlBasis = useDescriptors_ ? bv.getCurl2DVectorBasis(true) : Array(bv.weighted_curl_basis_scalar);
      integrate.evalStyle         = evalStyle_;
      parallel_for(workset.num_cells, integrate);
    }
    else // if (spaceDim_ == 3)
    {
      // Create an object to multiply the integrand by the scalar multiplier
      // and any field multipliers out in front of the integral.
      using PreMultiply         = PreMultiply3D<ScalarT>;
      using ScalarMultiplierTag = typename PreMultiply::ScalarMultiplierTag;
      using FieldMultiplierTag  = typename PreMultiply::FieldMultiplierTag;
      PreMultiply preMultiply;
      preMultiply.result      = result3D_;
      preMultiply.multiplier  = multiplier_;
      preMultiply.vectorField = vector3D_;

      // Multiply the integrand by the scalar multiplier out in front of the
      // integral.
      parallel_for(RangePolicy<Device, ScalarMultiplierTag>(0,
        workset.num_cells), preMultiply);

      // Multiply the integrand by any field multipliers out in front of the
      // integral.
      for (const auto& field : fieldMults_)
      {
        preMultiply.fieldMult = field;
        parallel_for(RangePolicy<Device, FieldMultiplierTag>(0,
          workset.num_cells), preMultiply);
      } // end loop over the field multipliers

      // Create an object to do the actual integration and then do it.
      Integrate3D<ScalarT, 3> integrate;
      integrate.result            = result3D_;
      integrate.field             = field_;
      using Array=typename BasisValues2<double>::ConstArray_CellBasisIPDim;
      integrate.weightedCurlBasis = useDescriptors_ ? bv.getCurlVectorBasis(true) : Array(bv.weighted_curl_basis_vector);
      integrate.evalStyle         = evalStyle_;
      parallel_for(workset.num_cells, integrate);
    } // end if spaceDim_ is 2 or 3
  } // end of evaluateFields()

  /////////////////////////////////////////////////////////////////////////////
  //
  //  getValidParameters()
  //
  /////////////////////////////////////////////////////////////////////////////
  template<typename EvalT, typename TRAITS>
  Teuchos::RCP<Teuchos::ParameterList>
  Integrator_CurlBasisDotVector<EvalT, TRAITS>::
  getValidParameters() const
  {
    using panzer::BasisIRLayout;
    using panzer::IntegrationRule;
    using std::string;
    using std::vector;
    using Teuchos::ParameterList;
    using Teuchos::RCP;
    using Teuchos::rcp;

    RCP<ParameterList> p = rcp(new ParameterList);
    p->set<string>("Residual Name", "?");
    p->set<string>("Value Name", "?");
    RCP<BasisIRLayout> basis;
    p->set("Basis", basis);
    RCP<IntegrationRule> ir;
    p->set("IR", ir);
    p->set<double>("Multiplier", 1.0);
    RCP<const vector<string>> fms;
    p->set("Field Multipliers", fms);
    return p;
  } // end of getValidParameters()

} // end of namespace panzer

#endif // Panzer_Integrator_CurlBasisDotVector_impl_hpp
