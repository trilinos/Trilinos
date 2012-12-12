//@HEADER
// ***********************************************************************
// 
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#ifndef IFPACK2_ADDITIVESCHWARZ_DECL_HPP
#define IFPACK2_ADDITIVESCHWARZ_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Tpetra_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Time.hpp"

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Ifpack2_ReorderFilter.hpp"
#include "Ifpack2_SingletonFilter.hpp"
#include "Ifpack2_OverlappingRowMatrix.hpp"

#ifdef HAVE_IFPACK2_ZOLTAN2
#include "Zoltan2_OrderingSolution.hpp"
#endif


namespace Ifpack2 {

/** \class AdditiveSchwarz
\brief Additive Schwarz domain decomposition for Tpetra::RowMatrix objects.

\section Ifpack2_AdditiveSchwarz_Summary Summary

This class implements an Additive Schwarz (one-level overlapping
domain decomposition) preconditioner.  It operates on a given
Tpetra::RowMatrix.  This class implements Tpetra::Operator, like all
other Ifpack2 Preconditioner subclasses.  Thus, the apply() method
applies the preconditioner to a multivector.

\section Ifpack2_AdditiveSchwarz_Alg Algorithm

One-level overlapping domain decomposition preconditioners use local
solvers of Dirichlet type. This means that the inverse of the local
matrix (with minimal or wider overlap) is applied to the residual to
be preconditioned.

The preconditioner can be written as:
\f[
P_{AS}^{-1} = \sum_{i=1}^M P_i A_i^{-1} R_i,
\f]
where \f$M\f$ is the number of subdomains (that is, the number of
processors in the computation), \f$R_i\f$ is an operator that
restricts the global vector to the vector lying on subdomain \f$i\f$,
\f$P_i\f$ is the prolongator operator, and

  \f[
  A_i = R_i A P_i.
  \f]

  The construction of Schwarz preconditioners is mainly composed by
  two steps:
  - definition of the restriction and prolongation operator
  \f$R_i\f$ and \f$R_i^T\f$. If minimal overlap is chosen, their
  implementation is trivial, \f$R_i\f$ will return all the local
  components. For wider overlaps, instead, Tpetra::Import and
  Tpetra::Export will be used to import/export data. The user
  must provide both the matrix to be preconditioned (which is suppose
  to have minimal-overlap) and the matrix with wider overlap.
  - definition of a technique to apply the inverse of \f$A_i\f$.
  To solve on each subdomain, the user can adopt any class, derived
  from Ifpack2::Preconditioner. This can be easily accomplished, as
  Ifpack2::AdditiveSchwarz is templated with the solver for each subdomain.

  The local matrix \f$A_i\f$ can be filtered, to eliminate singletons, and
  reordered. At the present time, RCM and METIS can be used to reorder the
  local matrix.
  
  The complete list of supported parameters is reported in page \ref ifp_params.

  \date Last modified on 05-Sep-12.
*/
template<class MatrixType,class LocalInverseType>
class AdditiveSchwarz : virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> {
      
public:
  typedef typename MatrixType::scalar_type         Scalar;
  typedef typename MatrixType::local_ordinal_type  LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type           Node;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  typedef typename Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>  LocalMatrixType;

  //@{ \name Constructors/Destructors
  //! Ifpack2::AdditiveSchwarz constructor with given Tpetra::RowMatrix.
  /*! Creates an Ifpack2::AdditiveSchwarz preconditioner with overlap.
   * To use minimal-overlap, OverlappingMatrix is omitted
   * (as defaulted to 0).
   *
   * \param
   * Matrix - (In) Pointer to matrix to be preconditioned
   *
   * \param
   * OverlappingMatrix - (In) Pointer to the matrix extended with the
   *                     desired level of overlap.
   */
  AdditiveSchwarz(const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & Matrix_in, int OverlapLevel_in = 0);
  
  //! Destructor
  virtual ~AdditiveSchwarz();
  //@}

  /** \name Methods implementing Tpetra::Operator. */
  //@{
  
  //! Returns the Map associated with the domain of this operator, which must be compatible with X.getMap().
  virtual const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getDomainMap() const;
  
  //! Returns the Map associated with the range of this operator, which must be compatible with Y.getMap().
  virtual const Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > & getRangeMap() const;
  
  //! Returns a pointer to the input matrix.
  virtual Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getMatrix() const;



  //! Applies the effect of the preconditioner.
  virtual void apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X, 
		     Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y, 
		     Teuchos::ETransp mode = Teuchos::NO_TRANS,
		     Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
		     Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;
  
  //@}
  
  //! Sets all parameters for the preconditioner.
  virtual void setParameters(const Teuchos::ParameterList& List);
  
  //! Computes all (graph-related) data necessary to initialize the preconditioner.
  virtual void initialize();
  
  //! Returns true if the  preconditioner has been successfully initialized, false otherwise.
  virtual bool isInitialized() const;
  
  //! Computes all (coefficient) data necessary to apply the preconditioner.
  virtual void compute();
  
  //! Returns true if the  preconditioner has been successfully computed, false otherwise.
  virtual bool isComputed() const;
  
  //! Computes the condition number estimate and returns its value.
  virtual magnitudeType computeCondEst(CondestType CT = Ifpack2::Cheap,
				       LocalOrdinal MaxIters = 1550,
				       magnitudeType Tol = 1e-9,
				       const Teuchos::Ptr<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > &Matrix = Teuchos::null);
  
  //! Returns the computed condition number estimate, or -1.0 if not computed.
  virtual magnitudeType getCondEst() const;
  
  //! Returns the number of calls to initialize().
  virtual int getNumInitialize() const;
  
  //! Returns the number of calls to compute().
  virtual int getNumCompute() const;
  
  //! Returns the number of calls to apply().
  virtual int getNumApply() const;
  
  //! Returns the time spent in initialize().
  virtual double getInitializeTime() const;
  
  //! Returns the time spent in compute().
  virtual double getComputeTime() const;

  //! Returns the time spent in apply().
  virtual double getApplyTime() const;

  //! @name Overridden from Teuchos::Describable 
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

  //! Prints basic information on iostream. This function is used by operator<<.
  virtual std::ostream& print(std::ostream& os) const;

  //! Returns the level of overlap.
  virtual int getOverlapLevel() const;

protected:

  // @}

  // @{ Internal merhods.
  
  //! Copy constructor (should never be used)
  AdditiveSchwarz(const AdditiveSchwarz& RHS);

  //! Sets up the localized matrix and the singleton filter.
  void setup();
  
  // @}

  // @{ Internal data.  
  //! Pointers to the matrix to be preconditioned.
  const Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Matrix_;

  //! Pointers to the overlapping matrix.
  Teuchos::RCP<Ifpack2::OverlappingRowMatrix<LocalMatrixType> >OverlappingMatrix_;

  //! Localized version of Matrix_ or OverlappingMatrix_.
  // CMS: Probably not a great idea, but this will remove the Local/Subdomain conflict here
  Teuchos::RCP<LocalMatrixType > LocalizedMatrix_;
  //! Pointer to the reordered matrix.
  Teuchos::RCP<Ifpack2::ReorderFilter<LocalMatrixType> > ReorderedLocalizedMatrix_;

  //! If true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! If true, overlapping is used
  bool IsOverlapping_;
  //! Level of overlap among the processors.
  int OverlapLevel_;
  //! Stores a copy of the list given in SetParameters()
  Teuchos::ParameterList List_;
  //! Combine mode for off-process elements (only if overlap is used)
  Tpetra::CombineMode CombineMode_;
  //! Contains the estimated condition number.
  magnitudeType Condest_;
  //! If \c true, compute the condition number estimate each time Compute() is called.
  bool ComputeCondest_;
  //! If \c true, reorder the local matrix.
  bool UseReordering_;
  //! If true, subdomain filtering is used
  bool UseSubdomain_;

  //! Filter for singletons.
  bool FilterSingletons_;
  //! filtering object.
  Teuchos::RCP<Ifpack2::SingletonFilter<LocalMatrixType> > SingletonMatrix_;
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to apply().
  mutable int NumApply_;
  //! Contains the time for all successful calls to initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to apply().
  mutable double ApplyTime_;
  //! Contains the number of flops for Initialize().
  double InitializeFlops_;
  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable double ApplyFlops_;
  //! Object used for timing purposes.
  Teuchos::RCP<Teuchos::Time> Time_;
  //! Pointer to the local solver.
  Teuchos::RCP<LocalInverseType> Inverse_;

}; // class AdditiveSchwarz

}// end namespace

#endif // IFPACK2_ADDITIVESCHWARZ_DECL_HPP
