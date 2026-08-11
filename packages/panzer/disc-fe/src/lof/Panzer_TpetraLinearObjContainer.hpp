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
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_x_th() const
   { return toThyraVector(x,domainSpace); }

   virtual void set_dxdt_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in)
   { 
     if(in==Teuchos::null) { dxdt = Teuchos::null; return; }

     Teuchos::RCP<const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > dxdt_const 
         = TOE::getConstTpetraVector(in);
     dxdt = Teuchos::rcp_const_cast<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(dxdt_const); 
   } 
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_dxdt_th() const
   { return toThyraVector(dxdt,domainSpace); }

   virtual void set_f_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in)
   { f = (in==Teuchos::null) ? Teuchos::null : TOE::getTpetraVector(in); } 
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_f_th() const
   { return toThyraVector(f,rangeSpace); }

   virtual void set_A_th(const Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > & in) 
   { A = (in==Teuchos::null) ? Teuchos::null : Teuchos::rcp_dynamic_cast<CrsMatrixType>(TOE::getTpetraOperator(in),true); }
   virtual Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > get_A_th() const
   { return (A==Teuchos::null) ? Teuchos::null : Thyra::createLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(A,rangeSpace,domainSpace); }
   //@}

private:
   /** \brief Wrap one of x/dxdt/f as a Thyra vector, whatever concrete kind it is.
     *
     * Thyra::createVector() needs a genuine Tpetra::Vector, but these fields are stored as
     * MultiVector so they can also hold a Tpetra::FEMultiVector (which derives from
     * MultiVector, not Vector). Two cases:
     *
     *  - a plain Vector: passed through with the container's own space, exactly as before,
     *    so the traditional path is bit-for-bit unchanged.
     *  - anything else (an FEMultiVector): handed back as getVectorNonConst(0), which
     *    Tpetra documents as "a Vector that views a single column" -- a genuine alias, so
     *    writes through the Thyra wrapper (e.g. Thyra::assign(...,0.0), which is how
     *    panzer::ModelEvaluator zeroes the ghosted residual) reach the real data.
     *
     * The FE case deliberately derives its space from the vector's CURRENT map instead of
     * reusing the container's stored space. An FEMultiVector's active map switches between
     * its owned and owned+shared views across beginAssembly()/endAssembly(), and
     * Thyra::createVector does NOT check the supplied space against the vector's map -- it
     * uses the space verbatim whenever one is given (see getOrCreateTpetraVectorSpace in
     * Thyra_TpetraThyraWrappers_def.hpp). Passing the stored space blindly would therefore
     * yield a Thyra vector silently claiming a space its data does not match, rather than
     * an error.
     */
   Teuchos::RCP<Thyra::VectorBase<ScalarT> >
   toThyraVector(const Teuchos::RCP<MultiVectorType> & mv,
                 const Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > & space) const
   {
     if(mv==Teuchos::null) return Teuchos::null;

     Teuchos::RCP<VectorType> v = Teuchos::rcp_dynamic_cast<VectorType>(mv);
     if(v!=Teuchos::null)
       return Thyra::createVector(v,space);

     // explicitly typed: a bare Teuchos::null cannot deduce the space template parameter
     const Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > deriveFromMap = Teuchos::null;
     return Thyra::createVector(mv->getVectorNonConst(0),deriveFromMap);
   }

   typedef Thyra::TpetraOperatorVectorExtraction<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> TOE;

   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > domainSpace;
   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > rangeSpace;

   Teuchos::RCP<Tpetra::MultiVector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > x, dxdt, f;
   Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > A;
};

}

#endif
