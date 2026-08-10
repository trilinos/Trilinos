// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_TpetraLinearObjContainer_hpp__
#define __Panzer_TpetraLinearObjContainer_hpp__

#include "PanzerDiscFE_config.hpp"

#include <map>

// Tpetra includes
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"

#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_NodeType.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

/**
 * \brief Tpetra-backed implementation of LinearObjContainer.
 *
 * Holds the solution vector `x`, its time derivative `dxdt`, the
 * residual vector `f`, and the Jacobian matrix `A` as Tpetra objects,
 * and (via ThyraObjContainer) exposes Thyra-wrapped views of the same
 * data so panzer's linear-algebra-agnostic code can operate on them
 * without depending on Tpetra directly. Built and populated by a
 * matching Tpetra linear object factory.
 * \tparam ScalarT the field scalar type.
 * \tparam LocalOrdinalT local ordinal type.
 * \tparam GlobalOrdinalT global ordinal type.
 * \tparam NodeT Kokkos node type; defaults to panzer::TpetraNodeType.
 */
template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT=panzer::TpetraNodeType>
class TpetraLinearObjContainer : public LinearObjContainer
                               , public ThyraObjContainer<ScalarT> {
   TpetraLinearObjContainer();

public:
   typedef LinearObjContainer::Members Members;

   typedef Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> VectorType;
   typedef Tpetra::MultiVector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> MultiVectorType;
   typedef Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> CrsMatrixType;
   typedef Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> CrsGraphType;
   typedef Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> MapType;
   typedef Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> ImportType;
   typedef Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> ExportType;

   /** \brief Construct from the domain (solution) and range (residual) Tpetra maps, building the corresponding Thyra vector spaces.
     *
     * \param[in] domain map for the solution/time-derivative vectors (x, dxdt) and the Jacobian's domain space.
     * \param[in] range map for the residual vector (f) and the Jacobian's range space.
     */
   TpetraLinearObjContainer(const Teuchos::RCP<const Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > & domain,
                            const Teuchos::RCP<const Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > & range)
   {
      domainSpace = Thyra::createVectorSpace<ScalarT>(domain);
      rangeSpace = Thyra::createVectorSpace<ScalarT>(range);
   }

   /// \brief Zeroes out whichever of x, dxdt, f, and A are currently set.
   //
   // Uses the _mv() accessors so this works whether the vectors are plain Tpetra::Vectors
   // or FEMultiVectors; putScalar() is a MultiVector-level method available on both.
   virtual void initialize()
   {
      if(get_x_mv()!=Teuchos::null) get_x_mv()->putScalar(0.0);
      if(get_dxdt_mv()!=Teuchos::null) get_dxdt_mv()->putScalar(0.0);
      if(get_f_mv()!=Teuchos::null) get_f_mv()->putScalar(0.0);
      if(get_A()!=Teuchos::null) {
        Teuchos::RCP<CrsMatrixType> mat = get_A();
        mat->setAllToScalar(0.0);
      }
   }

   //! Wipe out stored data.
   void clear()
   {
      set_x(Teuchos::null);
      set_dxdt(Teuchos::null);
      set_f(Teuchos::null);
      set_A(Teuchos::null);
   }

   //! \name Direct Tpetra accessors for the solution, time derivative, residual, and Jacobian.
   //
   // x/dxdt/f are *stored* as MultiVectorType (not VectorType) so that this same container
   // class can hold either a plain Tpetra::Vector (the traditional path) or a
   // Tpetra::FEMultiVector (the opt-in FE-assembly path) -- Tpetra::FEMultiVector derives
   // from Tpetra::MultiVector directly, not from Tpetra::Vector, so a Vector-typed field
   // could not hold one.
   //
   // The getters come in two flavors:
   //   get_x()    -- returns RCP<VectorType>, preserving the long-standing API for all
   //                 existing callers. Narrows via rcp_dynamic_cast, which returns null for
   //                 an unset field and throws loudly if the field actually holds an
   //                 FEMultiVector (rather than silently misbehaving).
   //   get_x_mv() -- returns the stored RCP<MultiVectorType>, valid for BOTH a plain Vector
   //                 and an FEMultiVector. Assembly code that must work in either mode
   //                 (the gather/scatter evaluators, the LOF's ghost<->global helpers)
   //                 uses this form.
   // Setters take MultiVectorType, which accepts either kind via implicit upcast, so
   // existing set_x(rcp_to_a_Vector) call sites are unaffected.
   //@{
   inline void set_x(const Teuchos::RCP<MultiVectorType> & in) { x = in; }
   inline const Teuchos::RCP<VectorType> get_x() const { return Teuchos::rcp_dynamic_cast<VectorType>(x,true); }
   inline const Teuchos::RCP<MultiVectorType> get_x_mv() const { return x; }

   inline void set_dxdt(const Teuchos::RCP<MultiVectorType> & in) { dxdt = in; }
   inline const Teuchos::RCP<VectorType> get_dxdt() const { return Teuchos::rcp_dynamic_cast<VectorType>(dxdt,true); }
   inline const Teuchos::RCP<MultiVectorType> get_dxdt_mv() const { return dxdt; }

   inline void set_f(const Teuchos::RCP<MultiVectorType> & in) { f = in; }
   inline const Teuchos::RCP<VectorType> get_f() const { return Teuchos::rcp_dynamic_cast<VectorType>(f,true); }
   inline const Teuchos::RCP<MultiVectorType> get_f_mv() const { return f; }

   inline void set_A(const Teuchos::RCP<CrsMatrixType> & in) { A = in; }
   inline const Teuchos::RCP<CrsMatrixType> get_A() const { return A; }
   //@}

   /// \brief Sets every entry of the Jacobian matrix A to \p value.
   void initializeMatrix(ScalarT value)
   {
     A->setAllToScalar(value);
   }

   //! \name ThyraObjContainer overrides: Thyra-wrapped views of the same underlying Tpetra data.
   //@{
   virtual void set_x_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in)
   {
     if(in==Teuchos::null) { x = Teuchos::null; return; }

     Teuchos::RCP<const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > x_const 
         = TOE::getConstTpetraVector(in);
     x = Teuchos::rcp_const_cast<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(x_const); 
   } 
   // Thyra::createVector() requires a genuine Tpetra::Vector, so these go through the
   // narrowing get_x()/get_dxdt()/get_f() accessors. That is always satisfied in practice:
   // the Thyra bridge (used to hand vectors to/from NOX/Piro) only ever touches the *global*
   // (owned) container, which never holds an FE object -- only the ghosted container does,
   // and it never crosses this bridge.
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_x_th() const
   { return (x==Teuchos::null) ? Teuchos::null : Thyra::createVector(get_x(),domainSpace); }

   virtual void set_dxdt_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in)
   { 
     if(in==Teuchos::null) { dxdt = Teuchos::null; return; }

     Teuchos::RCP<const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > dxdt_const 
         = TOE::getConstTpetraVector(in);
     dxdt = Teuchos::rcp_const_cast<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(dxdt_const); 
   } 
   // see get_x_th() -- same narrowing-back-to-Vector reasoning applies here.
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_dxdt_th() const
   { return (dxdt==Teuchos::null) ? Teuchos::null : Thyra::createVector(get_dxdt(),domainSpace); }

   virtual void set_f_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in)
   { f = (in==Teuchos::null) ? Teuchos::null : TOE::getTpetraVector(in); } 
   // see get_x_th() -- same narrowing-back-to-Vector reasoning applies here.
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_f_th() const
   { return (f==Teuchos::null) ? Teuchos::null : Thyra::createVector(get_f(),rangeSpace); }

   virtual void set_A_th(const Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > & in) 
   { A = (in==Teuchos::null) ? Teuchos::null : Teuchos::rcp_dynamic_cast<CrsMatrixType>(TOE::getTpetraOperator(in),true); }
   virtual Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > get_A_th() const
   { return (A==Teuchos::null) ? Teuchos::null : Thyra::createLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(A,rangeSpace,domainSpace); }
   //@}

private:
   typedef Thyra::TpetraOperatorVectorExtraction<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> TOE;

   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > domainSpace;
   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > rangeSpace;

   Teuchos::RCP<Tpetra::MultiVector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > x, dxdt, f;
   Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > A;
};

}

#endif
