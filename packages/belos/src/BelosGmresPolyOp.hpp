// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
#include "BelosOrthoManagerFactory.hpp"

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
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_ParameterList.hpp"

#ifdef BELOS_TEUCHOS_TIME_MONITOR
  #include "Teuchos_TimeMonitor.hpp"
#endif // BELOS_TEUCHOS_TIME_MONITOR

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
        RP_(problem_in->getRightPrec())
    {
      setParameters( params_ );

      polyUpdateLabel_ = label_ + ": Hybrid Gmres: Vector Update";
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      timerPolyUpdate_ = Teuchos::TimeMonitor::getNewCounter(polyUpdateLabel_);
#endif // BELOS_TEUCHOS_TIME_MONITOR

      if (polyType_ == "Arnoldi" || polyType_=="Roots")
        generateArnoldiPoly();
      else if (polyType_ == "Gmres")
        generateGmresPoly();
      else
        TEUCHOS_TEST_FOR_EXCEPTION(polyType_!="Arnoldi"&&polyType_!="Gmres"&&polyType_!="Roots",std::invalid_argument,
          "Belos::GmresPolyOp: \"Polynomial Type\" must be either \"Arnoldi\", \"Gmres\", or \"Roots\".");
    }

    //! Given no ParameterList, constructor creates no polynomial and only applies the given operator.
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
    void ApplyRootsPoly ( const MV& x, MV& y ) const;

    /*! \brief This routine casts the MultiVec to GmresPolyMv to retrieve the MV.  Then the above
        apply method is called.
    */
    void Apply ( const MultiVec<ScalarType>& x, MultiVec<ScalarType>& y, ETrans /* trans */=NOTRANS ) const
    {
      const GmresPolyMv<ScalarType,MV>& x_in = dynamic_cast<const GmresPolyMv<ScalarType,MV>&>(x);
      GmresPolyMv<ScalarType,MV>& y_in = dynamic_cast<GmresPolyMv<ScalarType,MV>&>(y);
      ApplyPoly( *(x_in.getConstMV()), *(y_in.getMV()) ); 
    }

    int polyDegree() const { return dim_; }

    private:

#ifdef BELOS_TEUCHOS_TIME_MONITOR
    Teuchos::RCP<Teuchos::Time> timerPolyUpdate_;
#endif // BELOS_TEUCHOS_TIME_MONITOR
    std::string polyUpdateLabel_;

    typedef int OT; //Ordinal type 
    typedef MultiVecTraits<ScalarType,MV> MVT;
    typedef Teuchos::ScalarTraits<ScalarType> SCT ;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef Teuchos::ScalarTraits<MagnitudeType> MCT ;

    // Default polynomial parameters
    static constexpr int maxDegree_default_ = 25;
    static constexpr int verbosity_default_ = Belos::Errors;
    static constexpr bool randomRHS_default_ = true;
    static constexpr const char * label_default_ = "Belos";
    static constexpr const char * polyType_default_ = "Roots";
    static constexpr const char * orthoType_default_ = "DGKS";
    static constexpr bool damp_default_ = false;
    static constexpr bool addRoots_default_ = true;

    // Variables for generating the polynomial
    Teuchos::RCP<LinearProblem<ScalarType,MV,OP> > problem_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<const OP> LP_, RP_;

    // Output manager.
    Teuchos::RCP<OutputManager<ScalarType> > printer_;
    Teuchos::RCP<std::ostream> outputStream_ = Teuchos::rcpFromRef(std::cout);
 
    // Orthogonalization manager.
    Teuchos::RCP<MatOrthoManager<ScalarType,MV,OP> > ortho_;

    // Current polynomial parameters
    MagnitudeType polyTol_ = DefaultSolverParameters::polyTol;
    int maxDegree_ = maxDegree_default_;
    int verbosity_ = verbosity_default_;
    bool randomRHS_ = randomRHS_default_;
    std::string label_ = label_default_;
    std::string polyType_ = polyType_default_;
    std::string orthoType_ = orthoType_default_;
    int dim_ = 0;
    bool damp_ = damp_default_;
    bool addRoots_ = addRoots_default_;
    
    // Variables for Arnoldi polynomial
    mutable Teuchos::RCP<MV> V_, wL_, wR_;
    Teuchos::SerialDenseMatrix<OT,ScalarType> H_, y_;
    Teuchos::SerialDenseVector<OT,ScalarType> r0_;

    // Variables for Gmres polynomial;
    bool autoDeg = false;
    Teuchos::SerialDenseMatrix< OT, ScalarType > pCoeff_;

    // Variables for Roots polynomial:
    Teuchos::SerialDenseMatrix< OT, MagnitudeType > theta_; 
    
    // Modified Leja sorting function. Takes a serial dense matrix of M harmonic Ritz values and an index
    // of values from 0 to M. Returns the sorted values and sorted index, similar to Matlab.
    void SortModLeja(Teuchos::SerialDenseMatrix< OT, MagnitudeType > &thetaN, std::vector<int> &index) const ;

    //Function determines whether added roots are needed and adds them if option is turned on.
    void ComputeAddedRoots();
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
    }

    // Check for timer label
    if (params_in->isParameter("Timer Label")) {
      label_ = params_in->get("Timer Label", label_default_);
    }
 
    // Output stream
    if (params_in->isParameter("Output Stream")) {
      outputStream_ = Teuchos::getParameter<Teuchos::RCP<std::ostream> >(*params_in,"Output Stream");
    }

    // Check for damped polynomial
    if (params_in->isParameter("Damped Poly")) {
      damp_ = params_in->get("Damped Poly", damp_default_);
    }

    // Check for root-adding
    if (params_in->isParameter("Add Roots")) {
      addRoots_ = params_in->get("Add Roots", addRoots_default_);
    }
  }

  template <class ScalarType, class MV, class OP>
  void GmresPolyOp<ScalarType, MV, OP>::generateGmresPoly()
  {
    Teuchos::RCP< MV > V = MVT::Clone( *problem_->getRHS(), maxDegree_+1 );   

    //Make power basis:
    std::vector<int> index(1,0);
    Teuchos::RCP< MV > V0 = MVT::CloneViewNonConst(*V, index);
    if (randomRHS_)
      MVT::MvRandom( *V0 );
    else
      MVT::Assign( *problem_->getRHS(), *V0 );

    if ( !LP_.is_null() ) {
      Teuchos::RCP< MV > Vtemp = MVT::CloneCopy(*V0);
      problem_->applyLeftPrec( *Vtemp, *V0);
    }
    if ( damp_ ) {
      Teuchos::RCP< MV > Vtemp = MVT::CloneCopy(*V0);
      problem_->apply( *Vtemp, *V0);
    }

    for(int i=0; i< maxDegree_; i++)
    {
      index[0] = i;
      Teuchos::RCP< const MV > Vi = MVT::CloneView(*V, index);
      index[0] = i+1;
      Teuchos::RCP< MV > Vip1 = MVT::CloneViewNonConst(*V, index);
      problem_->apply( *Vi, *Vip1);
    }

    //Consider AV:
    Teuchos::Range1D range( 1, maxDegree_);
    Teuchos::RCP< const MV > AV = MVT::CloneView( *V, range);
   
    //Make lhs (AV)^T(AV)
    Teuchos::SerialDenseMatrix< OT, ScalarType > AVtransAV( maxDegree_, maxDegree_);
    MVT::MvTransMv( SCT::one(), *AV, *AV, AVtransAV);
    //This process adds pDeg*pDeg + pDeg inner products that aren't in the final count.

    Teuchos::LAPACK< OT, ScalarType > lapack;
    int infoInt;
    bool status = true; //Keep adjusting poly deg when true.
   
    dim_ = maxDegree_; 
    Teuchos::SerialDenseMatrix< OT, ScalarType > lhs;
    while( status && dim_ >= 1)
    {
      Teuchos::SerialDenseMatrix< OT, ScalarType > lhstemp(Teuchos::Copy, AVtransAV,  dim_, dim_);
      lapack.POTRF( 'U', dim_, lhstemp.values(), lhstemp.stride(), &infoInt);
      
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
      pCoeff_.shape( dim_, 1);
      //Get correct submatrix of AV:
      Teuchos::Range1D rangeSub( 1, dim_);
      Teuchos::RCP< const MV > AVsub = MVT::CloneView( *V, rangeSub);
      
      //Compute rhs (AV)^T V0
      MVT::MvTransMv( SCT::one(), *AVsub, *V0, pCoeff_);
      lapack.POTRS( 'U', dim_, 1, lhs.values(), lhs.stride(), pCoeff_.values(), pCoeff_.stride(), &infoInt);
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

    // Create output manager.
    printer_ = Teuchos::rcp( new OutputManager<ScalarType>(verbosity_, outputStream_) );

    // Create orthogonalization manager if we need to.  
    if (ortho_.is_null()) {
      params_->set("Orthogonalization", orthoType_);
      Belos::OrthoManagerFactory<ScalarType, MV, OP> factory;
      Teuchos::RCP<Teuchos::ParameterList> paramsOrtho;   // can be null

      ortho_ = factory.makeMatOrthoManager (orthoType_, Teuchos::null, printer_, polyLabel, paramsOrtho);
    }

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
      newProblem->applyLeftPrec( *newB, *V_0 );
    if ( damp_ )
    {
      Teuchos::RCP< MV > Vtemp = MVT::CloneCopy(*V_0);
      newProblem->apply( *Vtemp, *V_0 );
    }
  
    // Get a matrix to hold the orthonormalization coefficients.
    r0_.resize(1);

    // Orthonormalize the new V_0
    int rank = ortho_->normalize( *V_0, Teuchos::rcpFromRef(r0_) );
    TEUCHOS_TEST_FOR_EXCEPTION(rank != 1,GmresPolyOpOrthoFailure,
      "Belos::GmresPolyOp::generateArnoldiPoly(): Failed to compute initial block of orthonormal vectors for polynomial generation.");
  
    // Set the new state and initialize the solver.
    GmresIterationState<ScalarType,MV> newstate;
    newstate.V = V_0;
    newstate.z = Teuchos::rcpFromRef( r0_);
    newstate.curDim = 0;
    gmres_iter->initializeGmres(newstate);

    // Perform Gmres iteration
    try {
      gmres_iter->iterate();
    }
    catch (GmresIterationOrthoFailure& e) {
      // Try to recover the most recent least-squares solution
      gmres_iter->updateLSQR( gmres_iter->getCurSubspaceDim() );
    }
    catch (std::exception& e) {
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
    if(polyType_ == "Arnoldi"){
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
    else{ //Generate Roots Poly
    //Find Harmonic Ritz Values to use as polynomial roots:
    
    //Copy of square H used to find poly roots:
    H_ = Teuchos::SerialDenseMatrix<OT,ScalarType>(Teuchos::Copy, *gmresState.H, dim_, dim_);
    //Zero out below subdiagonal of H:
    for(int i=0; i <= dim_-3; i++) {
      for(int k=i+2; k <= dim_-1; k++) {
        H_(k,i) = SCT::zero();
      }
    }
    //Extra copy of H because equilibrate changes the matrix: 
    Teuchos::SerialDenseMatrix<OT,ScalarType> Htemp (Teuchos::Copy, H_, dim_, dim_);

    //View the m+1,m element and last col of H:
    ScalarType Hlast = (*gmresState.H)(dim_,dim_-1);
    Teuchos::SerialDenseMatrix<OT,ScalarType> HlastCol (Teuchos::View, H_, dim_, 1, 0, dim_-1);

    //Set up linear system for H^{-*}e_m:
    Teuchos::SerialDenseMatrix< OT, ScalarType > F(dim_,1), E(dim_,1);
    E.putScalar(SCT::zero());
    E(dim_-1,0) = SCT::one();
      
    Teuchos::SerialDenseSolver< OT, ScalarType > HSolver;
    HSolver.setMatrix( Teuchos::rcpFromRef(Htemp));
    HSolver.solveWithTransposeFlag( Teuchos::CONJ_TRANS );
    HSolver.setVectors( Teuchos::rcpFromRef(F), Teuchos::rcpFromRef(E));
    HSolver.factorWithEquilibration( true );

    //Factor matrix and solve for F = H^{-*}e_m:
    int info = 0;
    info = HSolver.factor();
    if(info != 0){
      std::cout << "Hsolver factor: info = " << info << std::endl;
    }
    info = HSolver.solve();
    if(info != 0){
      std::cout << "Hsolver solve : info = " << info << std::endl;
    }

    //Scale F and adjust H for Harmonic Ritz value eigenproblem:
    F.scale(Hlast*Hlast); 
    HlastCol += F; 

    //Set up for eigenvalue problem to get Harmonic Ritz Values:
    Teuchos::LAPACK< OT, ScalarType > lapack;
    theta_.shape(dim_,2);//1st col for real part, 2nd col for imaginary 

    const int ldv = 1;
    ScalarType* vlr = 0;

    // Size of workspace and workspace for DGEEV
    int lwork = -1;
    std::vector<ScalarType> work(1);
    std::vector<MagnitudeType> rwork(2*dim_);

    //Find workspace size for DGEEV:
    lapack.GEEV('N','N',dim_,H_.values(),H_.stride(),theta_[0],theta_[1],vlr, ldv, vlr, ldv, &work[0], lwork, &rwork[0], &info);
    lwork = std::abs (static_cast<int> (Teuchos::ScalarTraits<ScalarType>::real (work[0])));
    work.resize( lwork );
    // Solve for Harmonic Ritz Values:
    lapack.GEEV('N','N',dim_,H_.values(),H_.stride(),theta_[0],theta_[1],vlr, ldv, vlr, ldv, &work[0], lwork, &rwork[0], &info);

    if(info != 0){
      std::cout << "GEEV solve : info = " << info << std::endl;
    }

    // Set index for sort function, verify roots are non-zero,
    // and sort Harmonic Ritz Values:
    const MagnitudeType tol = 10.0 * Teuchos::ScalarTraits<MagnitudeType>::eps();
    std::vector<int> index(dim_);
    for(int i=0; i<dim_; ++i){ 
      index[i] = i; 
      // Check if real + imag parts of roots < tol. 
      TEUCHOS_TEST_FOR_EXCEPTION(hypot(theta_(i,0),theta_(i,1)) < tol, std::runtime_error, "BelosGmresPolyOp Error: One of the computed polynomial roots is approximately zero.  This will cause a divide by zero error!  Your matrix may be close to singular.  Please select a lower polynomial degree or give a shifted matrix.");
    }
    SortModLeja(theta_,index);

    //Add roots if neded.
    ComputeAddedRoots();

   }
  }
  
  //Function determines whether added roots are needed and adds them if option is turned on.
  template <class ScalarType, class MV, class OP>
  void GmresPolyOp<ScalarType, MV, OP>::ComputeAddedRoots()
  {
    // Store theta (with cols for real and imag parts of Harmonic Ritz Vals) 
    // as one vector of complex numbers to perform arithmetic:
    std::vector<std::complex<MagnitudeType>> cmplxHRitz (dim_);
    for(unsigned int i=0; i<cmplxHRitz.size(); ++i){
      cmplxHRitz[i] = std::complex<MagnitudeType>( theta_(i,0), theta_(i,1) );
    }
   
    // Compute product of factors (pof) to determine added roots: 
    const MagnitudeType one(1.0);
    std::vector<MagnitudeType> pof (dim_,one);
    for(int j=0; j<dim_; ++j) {
      for(int i=0; i<dim_; ++i) {
        if(i!=j) {
          pof[j] = std::abs(pof[j]*(one-(cmplxHRitz[j]/cmplxHRitz[i])));
        }
      }
    }

    // Compute number of extra roots needed:
    std::vector<int> extra (dim_);
    int totalExtra = 0;
    for(int i=0; i<dim_; ++i){
      if (pof[i] > MCT::zero())
        extra[i] = ceil((log10(pof[i])-MagnitudeType(4.0))/MagnitudeType(14.0));
      else
        extra[i] = 0;
      if(extra[i] > 0){
        totalExtra += extra[i];
      }
    }
    if (totalExtra){
      printer_->stream(Warnings) << "Warning: Need to add " << totalExtra << " extra roots." << std::endl;}

    // If requested to add roots, append them to the theta matrix:
    if(addRoots_ && totalExtra>0)
    {
      theta_.reshape(dim_+totalExtra,2);
      // Make a matrix copy for perturbed roots:
      Teuchos::SerialDenseMatrix<OT,MagnitudeType> thetaPert (Teuchos::Copy, theta_, dim_+totalExtra, 2);

      //Add extra eigenvalues to matrix and perturb for sort:
      int count = dim_;
      for(int i=0; i<dim_; ++i){
        for(int j=0; j< extra[i]; ++j){
          theta_(count,0) = theta_(i,0);
          theta_(count,1) = theta_(i,1);
          thetaPert(count,0) = theta_(i,0)+(j+MCT::one())*MagnitudeType(5e-8);
          thetaPert(count,1) = theta_(i,1);
          ++count;
        }
      }

      // Update polynomial degree:
      dim_ += totalExtra;
      if (totalExtra){
        printer_->stream(Warnings) << "New poly degree is: " << dim_ << std::endl;}

      // Create a new index and sort perturbed roots:
      std::vector<int> index2(dim_);
      for(int i=0; i<dim_; ++i){ 
        index2[i] = i; 
      }
      SortModLeja(thetaPert,index2);
      //Apply sorting to non-perturbed roots:
      for(int i=0; i<dim_; ++i)
      { 
        thetaPert(i,0) = theta_(index2[i],0);
        thetaPert(i,1) = theta_(index2[i],1);
      }
      theta_ = thetaPert;

    } 
  }

  // Modified Leja sorting function. Takes a serial dense matrix of M harmonic Ritz values and an index
  // of values from 0 to M. Returns the sorted values and sorted index, similar to Matlab.
  template <class ScalarType, class MV, class OP>
  void GmresPolyOp<ScalarType, MV, OP>::SortModLeja(Teuchos::SerialDenseMatrix< OT, MagnitudeType > &thetaN, std::vector<int> &index) const 
  {
    //Sort theta values via Modified Leja Ordering:

    // Set up blank matrices to track sorting:
    int dimN = index.size();
    std::vector<int> newIndex(dimN);
    Teuchos::SerialDenseMatrix< OT, MagnitudeType > sorted (thetaN.numRows(), thetaN.numCols());
    Teuchos::SerialDenseVector< OT, MagnitudeType > absVal (thetaN.numRows());
    Teuchos::SerialDenseVector< OT, MagnitudeType > prod (thetaN.numRows());

    //Compute all absolute values and find maximum:
    for(int i = 0; i < dimN; i++){
      absVal(i) = hypot(thetaN(i,0), thetaN(i,1)); 
    }
    MagnitudeType * maxPointer = std::max_element(absVal.values(), (absVal.values()+dimN));
    int maxIndex = int (maxPointer- absVal.values());

    //Put largest abs value first in the list:
    sorted(0,0) = thetaN(maxIndex,0);
    sorted(0,1) = thetaN(maxIndex,1);
    newIndex[0] = index[maxIndex];

    int j;  
    // If largest value was complex (for real scalar type) put its conjugate in the next slot.
    if(sorted(0,1)!= SCT::zero() && !SCT::isComplex) 
    {
      sorted(1,0) = thetaN(maxIndex,0);
      sorted(1,1) = -thetaN(maxIndex,1);
      newIndex[1] = index[maxIndex+1];
      j = 2;
    }
    else
    {
      j = 1;
    }

    //Sort remaining values:
    MagnitudeType a, b;
    while( j < dimN )
    {
      //For each value, compute (a log of) a product of differences:
      for(int i = 0; i < dimN; i++) 
      {
        prod(i) = MCT::one();
        for(int k = 0; k < j; k++)
        {
          a = thetaN(i,0) - sorted(k,0);
          b = thetaN(i,1) - sorted(k,1);
          if (a*a + b*b > MCT::zero())
            prod(i) = prod(i) + log10(hypot(a,b));
          else {
            prod(i) = -std::numeric_limits<MagnitudeType>::infinity();
            break;
          }
        }
      }
      
      //Value with largest product goes in the next slot:
      maxPointer = std::max_element(prod.values(), (prod.values()+dimN));
      maxIndex = int (maxPointer- prod.values());
      sorted(j,0) = thetaN(maxIndex,0);
      sorted(j,1) = thetaN(maxIndex,1);
      newIndex[j] = index[maxIndex];

      //If it was complex (and scalar type real) put its conjugate in next slot:
      if(sorted(j,1)!= SCT::zero() && !SCT::isComplex) 
      {
        j++;
        sorted(j,0) = thetaN(maxIndex,0);
        sorted(j,1) = -thetaN(maxIndex,1);
        newIndex[j] = index[maxIndex+1];
      }
      j++;
    }

    //Return sorted values and sorted indices:
    thetaN = sorted;
    index = newIndex;
  } //End Modified Leja ordering

  template <class ScalarType, class MV, class OP>
  void GmresPolyOp<ScalarType, MV, OP>::ApplyPoly( const MV& x, MV& y ) const 
  {
    if (dim_) {
      if (polyType_ == "Arnoldi")
        ApplyArnoldiPoly(x, y);
      else if (polyType_ == "Gmres")
        ApplyGmresPoly(x, y);
      else if (polyType_ == "Roots")
        ApplyRootsPoly(x, y);
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

    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor updateTimer( *timerPolyUpdate_ );
#endif
      MVT::MvAddMv(pCoeff_(0,0), *AX, SCT::zero(), y, y); //y= coeff_i(A^ix)
    }
    for( int i=1; i < dim_; i++)
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
      {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor updateTimer( *timerPolyUpdate_ );
#endif
        MVT::MvAddMv(pCoeff_(i,0), *Y, SCT::one(), y, y); //y= coeff_i(A^ix) +y
      }
    }

    // Apply right preconditioner.
    if (!RP_.is_null()) {
      Teuchos::RCP<MV> Ytmp = MVT::CloneCopy(y);
      problem_->applyRightPrec( *Ytmp, y );
    }
  }

  template <class ScalarType, class MV, class OP>
  void GmresPolyOp<ScalarType, MV, OP>::ApplyRootsPoly( const MV& x, MV& y ) const 
  {
    MVT::MvInit( y, SCT::zero() ); //Zero out y to take the vector with poly applied.
    Teuchos::RCP<MV> prod = MVT::CloneCopy(x);
    Teuchos::RCP<MV> Xtmp = MVT::Clone( x, MVT::GetNumberVecs(x) );
    Teuchos::RCP<MV> Xtmp2 = MVT::Clone( x, MVT::GetNumberVecs(x) );

    // Apply left preconditioner.
    if (!LP_.is_null()) {
      problem_->applyLeftPrec( *prod, *Xtmp ); // Left precondition x into the first vector 
      prod = Xtmp;
    } 
    
    int i=0;   
    while(i < dim_-1)
    {
      if(theta_(i,1)== SCT::zero() || SCT::isComplex) //Real Harmonic Ritz value or complex scalars
      {
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerPolyUpdate_ );
#endif
          MVT::MvAddMv(SCT::one(), y, SCT::one()/theta_(i,0), *prod, y); //poly = poly + 1/theta_i * prod
        }
        problem_->apply(*prod, *Xtmp); // temp = A*prod
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerPolyUpdate_ );
#endif
          MVT::MvAddMv(SCT::one(), *prod, -SCT::one()/theta_(i,0), *Xtmp, *prod); //prod = prod - 1/theta_i * temp
        }
        i++;
      }
      else //Current theta is complex and has a conjugate; combine to preserve real arithmetic 
      {
        MagnitudeType mod = theta_(i,0)*theta_(i,0) + theta_(i,1)*theta_(i,1); //mod = a^2 + b^2
        problem_->apply(*prod, *Xtmp); // temp = A*prod
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerPolyUpdate_ );
#endif
          MVT::MvAddMv(2*theta_(i,0), *prod, -SCT::one(), *Xtmp, *Xtmp); //temp = 2a*prod-temp 
          MVT::MvAddMv(SCT::one(), y, SCT::one()/mod, *Xtmp, y); //poly = poly + 1/mod*temp
        }
        if( i < dim_-2 )
        {
          problem_->apply(*Xtmp, *Xtmp2); // temp2 = A*temp
          {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
            Teuchos::TimeMonitor updateTimer( *timerPolyUpdate_ );
#endif
            MVT::MvAddMv(SCT::one(), *prod, -SCT::one()/mod, *Xtmp2, *prod); //prod = prod - 1/mod * temp2
          }
        }
        i = i + 2; 
      }
    }
    if(theta_(dim_-1,1)== SCT::zero() || SCT::isComplex) 
    {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
      Teuchos::TimeMonitor updateTimer( *timerPolyUpdate_ );
#endif
      MVT::MvAddMv(SCT::one(), y, SCT::one()/theta_(dim_-1,0), *prod, y); //poly = poly + 1/theta_i * prod
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
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerPolyUpdate_ );
#endif
          MVT::MvTimesMatAddMv( -SCT::one(), *v_prev, h, SCT::one(), *v_next );
        }

        // Scale by H(i+1,i)
        MVT::MvScale( *v_next, SCT::one()/H_(i+1,i) );  
      }

      // Compute output y = V*y_./r0_
      if (!RP_.is_null()) {
        {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
          Teuchos::TimeMonitor updateTimer( *timerPolyUpdate_ );
#endif
          MVT::MvTimesMatAddMv( SCT::one()/r0_(0), *V_, y_, SCT::zero(), *wR_ );
        }
        problem_->applyRightPrec( *wR_, *y_view );
      } 
      else {
#ifdef BELOS_TEUCHOS_TIME_MONITOR
        Teuchos::TimeMonitor updateTimer( *timerPolyUpdate_ );
#endif
        MVT::MvTimesMatAddMv( SCT::one()/r0_(0), *V_, y_, SCT::zero(), *y_view );
      }
    } // (int j=0; j<n; ++j)
  } // end Apply()
} // end Belos namespace

#endif

// end of file BelosGmresPolyOp.hpp
