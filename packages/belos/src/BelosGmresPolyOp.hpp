//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef BELOS_GMRESPOLYOP_HPP
#define BELOS_GMRESPOLYOP_HPP

/*!     \file BelosGmresPolyOp.hpp
        \brief Defines the GMRES polynomial operator hybrid-GMRES iterative linear solver.
*/

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"

#include "BelosOperator.hpp"
#include "BelosMultiVec.hpp"
#include "BelosOperatorTraits.hpp"
#include "BelosMultiVecTraits.hpp"
#include "BelosLinearProblem.hpp"

#include "BelosGmresIteration.hpp"
#include "BelosBlockGmresIter.hpp"

#include "BelosDGKSOrthoManager.hpp"
#include "BelosICGSOrthoManager.hpp"
#include "BelosIMGSOrthoManager.hpp"

#include "BelosStatusTestMaxIters.hpp"
#include "BelosStatusTestGenResNorm.hpp"
#include "BelosStatusTestImpResNorm.hpp"
#include "BelosStatusTestCombo.hpp"
#include "BelosStatusTestOutputFactory.hpp"

#include "BelosOutputManager.hpp"

#include "Teuchos_BLAS.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_as.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Belos {
 
  /** \brief GmresPolyOpOrthoFailure is thrown when the orthogonalization manager is
   * unable to generate orthonormal columns from the initial basis vectors.
   *
   * This std::exception is thrown from the GmresPolyOp::generateArnoldiPoly() method.
   *
   */
  class GmresPolyOpOrthoFailure : public BelosError {public:
    GmresPolyOpOrthoFailure(const std::string& what_arg) : BelosError(what_arg)
      {}};

  // Create a shell class for the MV, inherited off MultiVec<> that will operate with the GmresPolyOp.
  template <class ScalarType, class MV>
  class GmresPolyMv : public MultiVec< ScalarType > 
  {
    public:

     GmresPolyMv ( const Teuchos::RCP<MV>& mv_in )
       : mv_(mv_in)
     {}
     GmresPolyMv ( const Teuchos::RCP<const MV>& mv_in )
     {
       mv_ = Teuchos::rcp_const_cast<MV>( mv_in );
     }
     Teuchos::RCP<MV> getMV() { return mv_; }
     Teuchos::RCP<const MV> getConstMV() const { return mv_; }

     GmresPolyMv * Clone ( const int numvecs ) const 
     {
       GmresPolyMv * newMV = new GmresPolyMv( MVT::Clone( *mv_, numvecs ) );
       return newMV;
     }
     GmresPolyMv * CloneCopy () const 
     {
       GmresPolyMv * newMV = new GmresPolyMv( MVT::CloneCopy( *mv_ ) );
       return newMV;
     }
     GmresPolyMv * CloneCopy ( const std::vector<int>& index ) const 
     {
       GmresPolyMv * newMV = new GmresPolyMv( MVT::CloneCopy( *mv_, index ) );
       return newMV;
     }
     GmresPolyMv * CloneViewNonConst ( const std::vector<int>& index ) 
     {
       GmresPolyMv * newMV = new GmresPolyMv( MVT::CloneViewNonConst( *mv_, index ) );
       return newMV;
     }
     const GmresPolyMv * CloneView ( const std::vector<int>& index ) const 
     {
       const GmresPolyMv * newMV = new GmresPolyMv( MVT::CloneView( *mv_, index ) );
       return newMV;
     }
     ptrdiff_t GetGlobalLength () const { return MVT::GetGlobalLength( *mv_ ); }
     int GetNumberVecs () const { return MVT::GetNumberVecs( *mv_ ); }
     void MvTimesMatAddMv (const ScalarType alpha,
                           const MultiVec<ScalarType>& A,
                           const Teuchos::SerialDenseMatrix<int,ScalarType>& B, const ScalarType beta) 
     {  
       const GmresPolyMv<ScalarType,MV>& A_in = dynamic_cast<const GmresPolyMv<ScalarType,MV>&>(A);
       MVT::MvTimesMatAddMv( alpha, *(A_in.getConstMV()), B, beta, *mv_ );
     }
     void MvAddMv ( const ScalarType alpha, const MultiVec<ScalarType>& A, const ScalarType beta, const MultiVec<ScalarType>& B ) 
     {
       const GmresPolyMv<ScalarType,MV>& A_in = dynamic_cast<const GmresPolyMv<ScalarType,MV>&>(A);
       const GmresPolyMv<ScalarType,MV>& B_in = dynamic_cast<const GmresPolyMv<ScalarType,MV>&>(B);
       MVT::MvAddMv( alpha, *(A_in.getConstMV()), beta, *(B_in.getConstMV()), *mv_ );
     }
     void MvScale ( const ScalarType alpha ) { MVT::MvScale( *mv_, alpha ); }
     void MvScale ( const std::vector<ScalarType>& alpha ) { MVT::MvScale( *mv_, alpha ); }
     void MvTransMv ( const ScalarType alpha, const MultiVec<ScalarType>& A, Teuchos::SerialDenseMatrix<int,ScalarType>& B) const 
     {
       const GmresPolyMv<ScalarType,MV>& A_in = dynamic_cast<const GmresPolyMv<ScalarType,MV>&>(A);
       MVT::MvTransMv( alpha, *(A_in.getConstMV()), *mv_, B );
     }
     void MvDot ( const MultiVec<ScalarType>& A, std::vector<ScalarType>& b ) const 
     {
       const GmresPolyMv<ScalarType,MV>& A_in = dynamic_cast<const GmresPolyMv<ScalarType,MV>&>(A);
       MVT::MvDot( *(A_in.getConstMV()), *mv_, b );
     }
     void MvNorm ( std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>& normvec, NormType type = TwoNorm ) const
     {
       MVT::MvNorm( *mv_, normvec, type ); 
     }
     void SetBlock ( const MultiVec<ScalarType>& A, const std::vector<int>& index ) 
     {
       const GmresPolyMv<ScalarType,MV>& A_in = dynamic_cast<const GmresPolyMv<ScalarType,MV>&>(A);
       MVT::SetBlock( *(A_in.getConstMV()), index, *mv_ );
     }
     void MvRandom () { MVT::MvRandom( *mv_ ); }
     void MvInit ( const ScalarType alpha ) { MVT::MvInit( *mv_, alpha ); }
     void MvPrint ( std::ostream& os ) const { MVT::MvPrint( *mv_, os ); }

    private:

      typedef MultiVecTraits<ScalarType,MV> MVT;

      Teuchos::RCP<MV> mv_;

   };
 
  /*!	\class GmresPolyOp

	\brief Belos's class for applying the GMRES polynomial operator that is used by the hybrid-GMRES linear solver.  

	This operator is used as the interface to the matrix polynomial (<tt>phi(A)</tt>), 
	solution (<tt>X</tt>), and right-hand side (<tt>B</tt>) of the linear system <tt>phi(A)X = B</tt>.
	Furthermore, it is also the interface to left/right preconditioning of the linear system.

	\author Heidi Thornquist and Jennifer Loe
  */
  template <class ScalarType, class MV, class OP>
  class GmresPolyOp : public Operator<ScalarType> {
  public:
    
    //! @name Constructor/Destructor
    //@{ 
    
    //! Basic contstructor 
    GmresPolyOp( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >& problem_in, 
                 const Teuchos::RCP<Teuchos::ParameterList>& params_in
               )
      : problem_(problem_in), 
        params_(params_in),
        LP_(problem_in->getLeftPrec()), 
        RP_(problem_in->getRightPrec()),
        outputStream_ (Teuchos::rcp(outputStream_default_,false)),
        polyTol_ (DefaultSolverParameters::polyTol),
        maxDegree_ (maxDegree_default_),
        verbosity_ (verbosity_default_),
        randomRHS_ (randomRHS_default_),
        label_ (label_default_),
        polyType_ (polyType_default_),
        orthoType_ (orthoType_default_),
        dim_(0)
    {
      setParameters( params_ );

      if (polyType_ == "Arnoldi")
        generateArnoldiPoly();
      else if (polyType_ == "Gmres")
        generateGmresPoly();
      else
        TEUCHOS_TEST_FOR_EXCEPTION(polyType_!="Arnoldi"&&polyType_!="Gmres",std::logic_error,
          "Belos::GmresPolyOp(): Invalid polynomial type.");
    }

    //! Create a simple polynomial of dimension 1.
    GmresPolyOp( const Teuchos::RCP<LinearProblem<ScalarType,MV,OP> >& problem_in )
      : problem_(problem_in)
    {
      // If dimension is zero, it will just apply the operator from problem_in in the Apply method.
      dim_ = 0;
    }

    //! Destructor.
    virtual ~GmresPolyOp() {};
    //@}
  
    //! @name Parameter processing method
    //@{
    
    /*! \brief Process the passed in parameters. */
    void setParameters( const Teuchos::RCP<Teuchos::ParameterList>& params_in );
    //@}
 
    //! @name Polynomial creation method
    //@{
   
    /*! \brief This routine takes the matrix, preconditioner, and vectors from the linear problem
        as well as the parameters to generate the Arnoldi polynomial.
    */
    void generateArnoldiPoly();

    /*! \brief This routine takes the matrix, preconditioner, and vectors from the linear problem
        as well as the parameters to generate the Gmres polynomial.
    */
    void generateGmresPoly();

    //@} 

    //! @name Operator application method
    //@{ 
    
    /*! \brief This routine takes the MV \c x and applies the polynomial 
        operator \c phi(OP) to it resulting in the MV \c y, which is returned.
        \note It is expected that any problem with applying this operator to \c x will be
	indicated by an std::exception being thrown.
    */
    void ApplyPoly ( const MV& x, MV& y ) const;
    void ApplyArnoldiPoly ( const MV& x, MV& y ) const;
    void ApplyGmresPoly ( const MV& x, MV& y ) const;

    /*! \brief This routine casts the MultiVec to GmresPolyMv to retrieve the MV.  Then the above
        apply method is called.
    */
    void Apply ( const MultiVec<ScalarType>& x, MultiVec<ScalarType>& y, ETrans trans=NOTRANS ) const 
    {
      const GmresPolyMv<ScalarType,MV>& x_in = dynamic_cast<const GmresPolyMv<ScalarType,MV>&>(x);
      GmresPolyMv<ScalarType,MV>& y_in = dynamic_cast<GmresPolyMv<ScalarType,MV>&>(y);
      ApplyPoly( *(x_in.getConstMV()), *(y_in.getMV()) ); 
    }

    int polyDegree() const { return dim_; }

    private:

    //! Default constructor
    GmresPolyOp() {}
   
    typedef int OT; //Ordinal type 
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT ;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;

    // Default polynomial parameters
    static constexpr int maxDegree_default_ = 25;
    static constexpr int verbosity_default_ = Belos::Errors;
    static constexpr bool randomRHS_default_ = true;
    static constexpr const char * label_default_ = "Belos";
    static constexpr const char * polyType_default_ = "Arnoldi";
    static constexpr const char * orthoType_default_ = "DGKS";
    static constexpr std::ostream * outputStream_default_ = &std::cout;

    // Variables for generating the polynomial
    Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<const OP> LP_, RP_;

    // Output manager.
    Teuchos::RCP<OutputManager<ScalarType> > printer_;
    Teuchos::RCP<std::ostream> outputStream_;
 
    // Orthogonalization manager.
    Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > ortho_;

    // Current polynomial parameters
    MagnitudeType polyTol_;
    int maxDegree_;
    int verbosity_;
    bool randomRHS_;
    std::string label_;
    std::string polyType_;
    std::string orthoType_;
    int dim_;

    // Variables for Arnoldi polynomial
    mutable Teuchos::RCP<MV> V_, wL_, wR_;
    Teuchos::SerialDenseMatrix<OT,ScalarType> H_, y_;
    Teuchos::SerialDenseVector<OT,ScalarType> r0_;

    // Variables for Gmres polynomial;
    bool autoDeg = false;
    Teuchos::SerialDenseMatrix< OT, ScalarType > pCoeff_;

  };
  
  template <class ScalarType, class MV, class OP>
  void GmresPolyOp<ScalarType, MV, OP>::setParameters( const Teuchos::RCP<Teuchos::ParameterList>& params_in )
  {
    // Check which Gmres polynomial to use
    if (params_in->isParameter("Polynomial Type")) {
      polyType_ = params_in->get("Polynomial Type", polyType_default_);
    }

    // Check for polynomial convergence tolerance
    if (params_in->isParameter("Polynomial Tolerance")) {
      if (params_in->isType<MagnitudeType> ("Polynomial Tolerance")) {
        polyTol_ = params_in->get ("Polynomial Tolerance",
                                 static_cast<MagnitudeType> (DefaultSolverParameters::polyTol));
      }     
      else {
        polyTol_ = params_in->get ("Polynomial Tolerance", DefaultSolverParameters::polyTol);
      }
    }

    // Check for maximum polynomial degree
    if (params_in->isParameter("Maximum Degree")) {
      maxDegree_ = params_in->get("Maximum Degree", maxDegree_default_);
    }

    // Check for maximum polynomial degree
    if (params_in->isParameter("Random RHS")) {
      randomRHS_ = params_in->get("Random RHS", randomRHS_default_);
    }
  
    // Check for a change in verbosity level
    if (params_in->isParameter("Verbosity")) {
      if (Teuchos::isParameterType<int>(*params_in,"Verbosity")) {
        verbosity_ = params_in->get("Verbosity", verbosity_default_);
      } 
      else {
        verbosity_ = (int)Teuchos::getParameter<Belos::MsgType>(*params_in,"Verbosity");
      }
    }

    if (params_in->isParameter("Orthogonalization")) {
      orthoType_ = params_in->get("Orthogonalization",orthoType_default_);
      TEUCHOS_TEST_FOR_EXCEPTION( orthoType_ != "DGKS" && orthoType_ != "ICGS" && orthoType_ != "IMGS",
                                  std::invalid_argument,
                                  "Belos::GmresPolyOp: \"Orthogonalization\" must be either \"DGKS\", \"ICGS\", or \"IMGS\".");
    }

    // Check for timer label
    if (params_in->isParameter("Timer Label")) {
      label_ = params_in->get("Timer Label", label_default_);
    }
 
    // Output stream
    if (params_in->isParameter("Output Stream")) {
      outputStream_ = Teuchos::getParameter<Teuchos::RCP<std::ostream> >(*params_in,"Output Stream");
    }
  }

  template <class ScalarType, class MV, class OP>
  void GmresPolyOp<ScalarType, MV, OP>::generateGmresPoly()
  {
    Teuchos::RCP< MV > V = MVT::Clone( *problem_->getRHS(), maxDegree_+2 );   

    //Make power basis:
    std::vector<int> index(1,0);
    Teuchos::RCP< MV > V0 = MVT::CloneViewNonConst(*V, index);
    if (randomRHS_)
      MVT::MvRandom( *V0 );
    else
      MVT::Assign( *problem_->getRHS(), *V0 );

    if ( !LP_.is_null() )
    {
      Teuchos::RCP< MV > Vtemp = MVT::CloneCopy(*V0);
      problem_->applyLeftPrec( *Vtemp, *V0);
    }

    for(int i=0; i< maxDegree_+1; i++)
    {
      index[0] = i;
      Teuchos::RCP< const MV > Vi = MVT::CloneView(*V, index);
      index[0] = i+1;
      Teuchos::RCP< MV > Vip1 = MVT::CloneViewNonConst(*V, index);
      problem_->apply( *Vi, *Vip1);
    }

    //Consider AV:
    Teuchos::Range1D range( 1, maxDegree_+1);
    Teuchos::RCP< const MV > AV = MVT::CloneView( *V, range);
   
    //Make lhs (AV)^T(AV)
    Teuchos::SerialDenseMatrix< OT, ScalarType > AVtransAV( maxDegree_+1, maxDegree_+1);
    MVT::MvTransMv( SCT::one(), *AV, *AV, AVtransAV);
    //This process adds pDeg*pDeg + pDeg inner products that aren't in the final count.

    Teuchos::LAPACK< OT, ScalarType > lapack;
    int infoInt;
    bool status = true; //Keep adjusting poly deg when true.
   
    dim_ = maxDegree_; 
    Teuchos::SerialDenseMatrix< OT, ScalarType > lhs;
    while( status && dim_ >= 1)
    {
      Teuchos::SerialDenseMatrix< OT, ScalarType > lhstemp(Teuchos::Copy, AVtransAV,  dim_+1, dim_+1);
      lapack.POTRF( 'U', dim_+1, lhstemp.values(), lhstemp.stride(), &infoInt);
      
      if(autoDeg == false)
      { 
        status = false;
        if(infoInt != 0)
        {
          std::cout << "BelosGmresPolyOp.hpp: LAPACK POTRF was not successful!!" << std::endl;
          std::cout << "Error code: " << infoInt << std::endl; 
        } 
      } 
      else
      {
        if(infoInt != 0)
        {//Had bad factor.  Reduce poly degree.
          dim_--;
        } 
        else
        {
          status = false;
        } 
      } 
      if(status == false)
      {
        lhs = lhstemp;
      } 
    } 
    if(dim_ == 0)
    {
      pCoeff_.shape( 1, 1);
      pCoeff_(0,0) = SCT::one(); 
      std::cout << "Poly Degree is zero.  No preconditioner created." << std::endl;
    } 
    else
    {
      pCoeff_.shape( dim_+1, 1);
      //Get correct submatrix of AV:
      Teuchos::Range1D rangeSub( 1, dim_+1);
      Teuchos::RCP< const MV > AVsub = MVT::CloneView( *V, rangeSub);
      
      //Compute rhs (AV)^T V0
      MVT::MvTransMv( SCT::one(), *AVsub, *V0, pCoeff_);
      lapack.POTRS( 'U', dim_+1, 1, lhs.values(), lhs.stride(), pCoeff_.values(), pCoeff_.stride(), &infoInt);
      if(infoInt != 0) 
      {
        std::cout << "BelosGmresPolyOp.hpp: LAPACK POTRS was not successful!!" << std::endl;
        std::cout << "Error code: " << infoInt << std::endl; 
      } 
    } 
  }

  template <class ScalarType, class MV, class OP>
  void GmresPolyOp<ScalarType, MV, OP>::generateArnoldiPoly()
  {
    std::string polyLabel = label_ + ": GmresPolyOp creation";

    // Create a copy of the linear problem that has a zero initial guess and random RHS.
    std::vector<int> idx(1,0);
    Teuchos::RCP<MV> newX  = MVT::Clone( *(problem_->getLHS()), 1 );
    Teuchos::RCP<MV> newB  = MVT::Clone( *(problem_->getRHS()), 1 );
    MVT::MvInit( *newX, SCT::zero() );
    if (randomRHS_) {
      MVT::MvRandom( *newB );
    }
    else {
      MVT::Assign( *(MVT::CloneView(*(problem_->getRHS()), idx)), *newB );
    }
    Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > newProblem =
      Teuchos::rcp( new LinearProblem<ScalarType,MV,OP>( problem_->getOperator(), newX, newB ) );
    newProblem->setInitResVec( newB );
    newProblem->setLeftPrec( problem_->getLeftPrec() );
    newProblem->setRightPrec( problem_->getRightPrec() );
    newProblem->setLabel(polyLabel);
    newProblem->setProblem();
    newProblem->setLSIndex( idx );

    // Create a parameter list for the GMRES iteration.
    Teuchos::ParameterList polyList;
 
    // Tell the block solver that the block size is one.
    polyList.set("Num Blocks",maxDegree_);
    polyList.set("Block Size",1);
    polyList.set("Keep Hessenberg", true);

    // Create orthogonalization manager if we need to.  
    if (ortho_.is_null()) {
      params_->set("Orthogonalization", orthoType_);
      if (orthoType_=="DGKS") {
        ortho_ = Teuchos::rcp( new DGKSOrthoManager<ScalarType,MV,OP>( polyLabel ) );
      }
      else if (orthoType_=="ICGS") {
        ortho_ = Teuchos::rcp( new ICGSOrthoManager<ScalarType,MV,OP>( polyLabel ) );
      }
      else if (orthoType_=="IMGS") {
        ortho_ = Teuchos::rcp( new IMGSOrthoManager<ScalarType,MV,OP>( polyLabel ) );
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(orthoType_!="ICGS"&&orthoType_!="DGKS"&&orthoType_!="IMGS",std::logic_error,
          "Belos::GmresPolyOp(): Invalid orthogonalization type.");
      }
    }

    // Create output manager.
    printer_ = Teuchos::rcp( new OutputManager<ScalarType>(verbosity_, outputStream_) );

    // Create a simple status test that either reaches the relative residual tolerance or maximum polynomial size.
    Teuchos::RCP<StatusTestMaxIters<ScalarType,MV,OP> > maxItrTst =
      Teuchos::rcp( new StatusTestMaxIters<ScalarType,MV,OP>( maxDegree_ ) );
  
    // Implicit residual test, using the native residual to determine if convergence was achieved.
    Teuchos::RCP<StatusTestGenResNorm<ScalarType,MV,OP> > convTst =
      Teuchos::rcp( new StatusTestGenResNorm<ScalarType,MV,OP>( polyTol_ ) );
    convTst->defineScaleForm( convertStringToScaleType("Norm of RHS"), Belos::TwoNorm );
  
    // Convergence test that stops the iteration when either are satisfied.
    Teuchos::RCP<StatusTestCombo<ScalarType,MV,OP> > polyTest =
      Teuchos::rcp( new StatusTestCombo<ScalarType,MV,OP>( StatusTestCombo<ScalarType,MV,OP>::OR, maxItrTst, convTst ) );
  
    // Create Gmres iteration object to perform one cycle of Gmres.
    Teuchos::RCP<BlockGmresIter<ScalarType,MV,OP> > gmres_iter;
    gmres_iter = Teuchos::rcp( new BlockGmresIter<ScalarType,MV,OP>(newProblem,printer_,polyTest,ortho_,polyList) );
 
    // Create the first block in the current Krylov basis (residual).
    Teuchos::RCP<MV> V_0 = MVT::CloneCopy( *newB );
    if ( !LP_.is_null() )
      problem_->applyLeftPrec( *newB, *V_0 );
  
    // Get a matrix to hold the orthonormalization coefficients.
    r0_.resize(1);

    // Orthonormalize the new V_0
    int rank = ortho_->normalize( *V_0, Teuchos::rcp( &r0_, false ) );
    TEUCHOS_TEST_FOR_EXCEPTION(rank != 1,GmresPolyOpOrthoFailure,
      "Belos::GmresPolyOp::generateArnoldiPoly(): Failed to compute initial block of orthonormal vectors for polynomial generation.");
  
    // Set the new state and initialize the solver.
    GmresIterationState<ScalarType,MV> newstate;
    newstate.V = V_0;
    newstate.z = Teuchos::rcp( &r0_,false );
    newstate.curDim = 0;
    gmres_iter->initializeGmres(newstate);

    // Perform Gmres iteration
    try {
      gmres_iter->iterate();
    }
    catch (GmresIterationOrthoFailure e) {
      // Try to recover the most recent least-squares solution
      gmres_iter->updateLSQR( gmres_iter->getCurSubspaceDim() );
    }
    catch (std::exception e) {
      using std::endl;
      printer_->stream(Errors) << "Error! Caught exception in BlockGmresIter::iterate() at iteration "
        << gmres_iter->getNumIters() << endl << e.what () << endl;
      throw;
    }

    // Get the solution for this polynomial, use in comparison below
    Teuchos::RCP<MV> currX = gmres_iter->getCurrentUpdate();
  
    // Record polynomial info, get current GMRES state
    GmresIterationState<ScalarType,MV> gmresState = gmres_iter->getState();
  
    // If the polynomial has no dimension, the tolerance is too low, return false
    dim_ = gmresState.curDim;
    if (dim_ == 0) {
      return;
    }
 
    //  Make a view and then copy the RHS of the least squares problem.
    //
    y_ = Teuchos::SerialDenseMatrix<OT,ScalarType>( Teuchos::Copy, *gmresState.z, dim_, 1 );
    H_ = *gmresState.H;

    //
    // Solve the least squares problem.
    //
    Teuchos::BLAS<OT,ScalarType> blas;
    blas.TRSM( Teuchos::LEFT_SIDE, Teuchos::UPPER_TRI, Teuchos::NO_TRANS,
               Teuchos::NON_UNIT_DIAG, dim_, 1, SCT::one(),
               gmresState.R->values(), gmresState.R->stride(),
               y_.values(), y_.stride() );

  }

  template <class ScalarType, class MV, class OP>
  void GmresPolyOp<ScalarType, MV, OP>::ApplyPoly( const MV& x, MV& y ) const 
  {
    if (dim_) {
      if (polyType_ == "Arnoldi")
        ApplyArnoldiPoly(x, y);
      else if (polyType_ == "Gmres")
        ApplyGmresPoly(x, y);
    }
    else {
      // Just apply the operator in problem_ to x and return y.
      problem_->applyOp( x, y );
    }  
  }

  template <class ScalarType, class MV, class OP>
  void GmresPolyOp<ScalarType, MV, OP>::ApplyGmresPoly( const MV& x, MV& y ) const 
  {
    Teuchos::RCP<MV> AX = MVT::CloneCopy(x);
    Teuchos::RCP<MV> AX2 = MVT::Clone( x, MVT::GetNumberVecs(x) );

    // Apply left preconditioner.
    if (!LP_.is_null()) {
      Teuchos::RCP<MV> Xtmp = MVT::Clone( x, MVT::GetNumberVecs(x) );
      problem_->applyLeftPrec( *AX, *Xtmp ); // Left precondition x into the first vector 
      AX = Xtmp;
    } 

    MVT::MvAddMv(pCoeff_(0,0), *AX, SCT::zero(), y, y); //y= coeff_i(A^ix)
    for( int i=1; i < dim_+1; i++)
    {
      Teuchos::RCP<MV> X, Y;
      if ( i%2 )
      {
        X = AX;
        Y = AX2;
      }
      else
      {
        X = AX2;
        Y = AX;
      }
      problem_->apply(*X, *Y);
      MVT::MvAddMv(pCoeff_(i,0), *Y, SCT::one(), y, y); //y= coeff_i(A^ix) +y
    }

    // Apply right preconditioner.
    if (!RP_.is_null()) {
      Teuchos::RCP<MV> Ytmp = MVT::CloneCopy(y);
      problem_->applyRightPrec( *Ytmp, y );
    }
  }

  template <class ScalarType, class MV, class OP>
  void GmresPolyOp<ScalarType, MV, OP>::ApplyArnoldiPoly( const MV& x, MV& y ) const 
  {
    // Initialize vector storage.
    if (V_.is_null()) {
      V_ = MVT::Clone( x, dim_ );
      if (!LP_.is_null()) {
        wL_ = MVT::Clone( y, 1 );
      }
      if (!RP_.is_null()) {
        wR_ = MVT::Clone( y, 1 );
      }
    }
    //
    // Apply polynomial to x.
    // 
    int n = MVT::GetNumberVecs( x );
    std::vector<int> idxi(1), idxi2, idxj(1);

    // Select vector x[j].
    for (int j=0; j<n; ++j) {

      idxi[0] = 0;
      idxj[0] = j;
      Teuchos::RCP<const MV> x_view = MVT::CloneView( x, idxj );
      Teuchos::RCP<MV> y_view = MVT::CloneViewNonConst( y, idxj );
      if (!LP_.is_null()) {
        Teuchos::RCP<MV> v_curr = MVT::CloneViewNonConst( *V_, idxi );
        problem_->applyLeftPrec( *x_view, *v_curr ); // Left precondition x into the first vector of V
      } else {
        MVT::SetBlock( *x_view, idxi, *V_ );  // Set x as the first vector of V
      }

      for (int i=0; i<dim_-1; ++i) {

        // Get views into the current and next vectors
        idxi2.resize(i+1);
        for (int ii=0; ii<i+1; ++ii) { idxi2[ii] = ii; }
        Teuchos::RCP<const MV> v_prev = MVT::CloneView( *V_, idxi2 );
        // the tricks below with wR_ and wL_ (potentially set to v_curr and v_next) unfortunately imply that 
        // v_curr and v_next must be non-const views.
        Teuchos::RCP<MV> v_curr = MVT::CloneViewNonConst( *V_, idxi );
        idxi[0] = i+1;
        Teuchos::RCP<MV> v_next = MVT::CloneViewNonConst( *V_, idxi );

        //---------------------------------------------
        // Apply operator to next vector
        //---------------------------------------------
        // 1) Apply right preconditioner, if we have one.
        if (!RP_.is_null()) {
          problem_->applyRightPrec( *v_curr, *wR_ );
        } else {
          wR_ = v_curr;
        }
        // 2) Check for left preconditioner, if none exists, point at the next vector.
        if (LP_.is_null()) {
          wL_ = v_next;
        }
        // 3) Apply operator A.
        problem_->applyOp( *wR_, *wL_ );
        // 4) Apply left preconditioner, if we have one.
        if (!LP_.is_null()) {
          problem_->applyLeftPrec( *wL_, *v_next );
        }

        // Compute A*v_curr - v_prev*H(1:i,i)
        Teuchos::SerialDenseMatrix<OT,ScalarType> h(Teuchos::View,H_,i+1,1,0,i);
        MVT::MvTimesMatAddMv( -SCT::one(), *v_prev, h, SCT::one(), *v_next );

        // Scale by H(i+1,i)
        MVT::MvScale( *v_next, SCT::one()/H_(i+1,i) );  
      }

      // Compute output y = V*y_./r0_
      if (!RP_.is_null()) {
        MVT::MvTimesMatAddMv( SCT::one()/r0_(0), *V_, y_, SCT::zero(), *wR_ );
        problem_->applyRightPrec( *wR_, *y_view );
      } 
      else {
        MVT::MvTimesMatAddMv( SCT::one()/r0_(0), *V_, y_, SCT::zero(), *y_view );
      }
    } // (int j=0; j<n; ++j)
  } // end Apply()
} // end Belos namespace

#endif

// end of file BelosGmresPolyOp.hpp
