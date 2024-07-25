// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __IFPACK2_BASKERILU_HPP__
#define __IFPACK2_BASKERILU_HPP__
//Include needed Trilinos
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_KokkosCompat_DefaultNode.hpp>
#include <Kokkos_CrsMatrix.hpp>
#include <Ifpack2_Preconditioner.hpp>

//Include needed Basker
#include "../src/shylubasker_decl.hpp"

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class ExecSpace>
class Baskerilu : public Ifpack2::Preconditioner<Scalar,LocalOrdinal, GlobalOrdinal, Node>
{
public:
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal,GlobalOrdinal,Node> TCrsMatrix;
  typedef Kokkos::CrsMatrix<Scalar, LocalOrdinal,GlobalOrdinal,Node> KCrsMatrix;
  typedef Kokkos::View<LocalOrdinal*, ExecSpace> OrdinalArray;
  typedef Kokkos::View<Scalar*, ExecSpace> ScalarArray;
 
private:
  bool initFlag;
  bool computedFlag;
  int nInit;
  int nApply;
  int nComputed;
  double initTime;
  double computeTime;
  double applyTime;
  Teuchos::RCP<TCrsMatrix> mat;
  
public:
  Baskerilu()
  {
    //Made to fit the needs of your solve
  }
  
  //required by IfPack2
  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  getDomainMap() const
  {
    return mat->getDomainMap();
  }

  Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
  getRangeMap() const
  {
    return mat->getRangeMap();
  }

  void
  apply
  (
   const Tpetra::MultVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
   Teuchos::ETransp mode = Teuchos::NO_TRANS,
   Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
   Scalar beta   =Teuchos::ScalarTraits<Scalar>::two()
   )const
  {

    return;
  }//end of apply()

  void setParameters(const Teuchos::ParameterList& List)
  {
     return;
  }//end setParameters()
	
  void initialize()
  {
    
    return;
  }//end initialize()

  bool isInitialized() const
  {
    return false;
  }//end isInitialized()

  void compute()
  {
    return;
  }//end compute()

  bool isComputed() const
  {
    return false;
  }//end isComputed()

  Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > 
  getMatrix() const
  {

    
  }//end getMatrix()

  int getNumInitialize() const
  {

  }//end getNumInitialize()

  int getNumCompute() const
  {  
    
  }//end getNumCompute()

  int getNumApply() const
  {

  }//end getNumApply()

  double getInitializeTime() const
  {

  }//end getInitializeTime

  double getComputeTime() const
  {

  }//end getComputeTime

  double getApplyTime() const
  {
    
  }//end getApplyTime()

  void checkLocalILU()
  {
    
  }//end checkLocalILU()

  void checkLocalIC()
  {
    
  }//end checkIC()

  void printStatus()
  {

  }//end printStat()

};//end class Baskerilu

#endif //end __IFPACK2_BASKERILU_HPP__
