// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
 */

#ifndef IFPACK2_ILU_HPP
#define IFPACK2_ILU_HPP

#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Condest.hpp"
#include "Ifpack2_ScalingType.hpp"
#include "Ifpack2_IlukGraph.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Teuchos_Comm.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Teuchos_Time.hpp"

#include "Ifpack2_CondestType.hpp"
#include "Teuchos_ParameterList.hpp"

// Define this macro to see some timers for some of these functions
#define ENABLE_IFPACK2_ILU_TEUCHOS_TIMERS

#ifdef ENABLE_IFPACK2_ILU_TEUCHOS_TIMERS
#  include "Teuchos_TimeMonitor.hpp"
#endif

namespace Teuchos {
  class ParameterList;
}
namespace Ifpack2 {
	//! ILU: A class for constructing and using an incomplete lower/upper (ILU) factorization of a given Tpetra::RowMatrix.
	
	/*! The ILU class computes a "Relaxed" ILU factorization with level k fill of a given Tpetra::RowMatrix. 
	 <P> Please refer to \ref ifp_ilu for a general description of the ILU algorithm.
	 <P>The complete list of supported parameters is reported in page \ref ifp_params. 
	 \author Mike Heroux, SNL 1416.
	 \date Last modified on 11-May-2009
	 */    
	template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>	
	class ILU: public Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> {
	
public:
	// @{ Constructors and destructors.
	//! Constructor
	ILU(Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> * A);
	
	//! Destructor
	~ILU()
	{
		destroy();
	}
	
	// @}
	// @{ Construction methods
	
	//! Initialize the preconditioner, does not touch matrix values.
	int initialize();
	
	//! Returns \c true if the preconditioner has been successfully initialized.
	bool isInitialized() const
	{
		return(IsInitialized_);
	}
	
	//! Compute ILU factors L and U using the specified graph, diagonal perturbation thresholds and relaxation parameters.
	/*! This function computes the ILU(k) factors L and U using the current:
	 <ol>
	 <li> Ifpack2_IlukGraph specifying the structure of L and U.
	 <li> Value for the ILU(k) relaxation parameter.
	 <li> Value for the \e a \e priori diagonal threshold values.
	 </ol>
	 InitValues() must be called before the factorization can proceed.
	 */
	int compute();
	
	//! If factor is completed, this query returns true, otherwise it returns false.
	bool isComputed() const 
	{
		return(isComputed_);
	}
	
	//! Set parameters using a Teuchos::ParameterList object.
	/* This method is only available if the Teuchos package is enabled.
	 This method recognizes four parameter names: relax_value,
	 absolute_threshold, relative_threshold and overlap_mode. These names are
	 case insensitive, and in each case except overlap_mode, the ParameterEntry
	 must have type double. For overlap_mode, the ParameterEntry must have
	 type Tpetra::CombineMode.
	 */
	int setParameters(const Teuchos::ParameterList& parameterlist);
	
	//! If set true, transpose of this operator will be applied.
	/*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
	 affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	 does not support transpose use, this method should return a value of -1.
	 
	 \param
	 UseTranspose_in - (In) If true, multiply by the transpose of operator, otherwise just use operator.
	 
	 \return Always returns 0.
	 */
	int setUseTranspose(bool UseTranspose_in) {UseTranspose_ = UseTranspose_in; return(0);};
	// @}
	
	// @{ Mathematical functions.
	// Applies the matrix to X, returns the result in Y.
	int apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
			  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
	{
		return(Multiply(false,X,Y));
	}
	
	int multiply(bool Trans, const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
				 Tpetra::Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;
	
	//! Returns the result of a Tpetra::Operator inverse applied to an Tpetra::MultiVector X in Y.
	/*! In this implementation, we use several existing attributes to determine how virtual
	 method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(), 
	 the Tpetra::CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
	 if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
	 be accounted for when doing a triangular solve.
	 
	 \param 
	 X - (In) A Tpetra::MultiVector of dimension NumVectors to solve for.
	 \param Out
	 Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing result.
	 
	 \return Integer error code, set to 0 if successful.
	 */
	int applyInverse(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
					 Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;
	
	//! Computes the estimated condition number and returns the value.
	double getCondest(const Ifpack2_CondestType CT = Ifpack2_Cheap, 
					  const int MaxIters = 1550,
					  const double Tol = 1e-9,
					  Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* Matrix_in = 0);
	
	//! Returns the computed estimated condition number, or -1.0 if not computed.
	double getCondest() const
	{
		return(condest_);
	}
	
	// @}
	// @{ Query methods
	
	//! Returns the address of the L factor associated with this factored matrix.
	const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & getL() const {return(*L_);};
	
	//! Returns the address of the D factor associated with this factored matrix.
	const Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> & getD() const {return(*D_);};
	
	//! Returns the address of the L factor associated with this factored matrix.
	const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & getU() const {return(*U_);};
	
	//! Returns a character string describing the operator
	const char* getLabel() const {return(Label_);}
	
	//! Sets label for \c this object.
	void setLabel(const char* Label_in)
	{
		strcpy(Label_,Label_in);
		return(0);
	}
	
	//! Returns 0.0 because this class cannot compute Inf-norm.
	Teuchos::ScalarTraits<Scalar>::magnitudeType normInf() const {return(0.0);};
	
	//! Returns false because this class cannot compute an Inf-norm.
	bool hasNormInf() const {return(false);};
	
	//! Returns the current UseTranspose setting.
	bool useTranspose() const {return(UseTranspose_);};
	
	//! Returns the Tpetra::Map object associated with the domain of this operator.
	const Tpetra::Map<Scalar,LocalOrdinal,GlobalOrdinal,Node> & getOperatorDomainMap() const {return(U_->OperatorDomainMap());};
	
	//! Returns the Tpetra::Map object associated with the range of this operator.
	const Tpetra::Map<Scalar,LocalOrdinal,GlobalOrdinal,Node> & getOperatorRangeMap() const{return(L_->OperatorRangeMap());};
	
	//! Returns a reference to the matrix to be preconditioned.
	const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& getUserMatrix() const
	{ 
		return(*A_);
	}
	
	//! Returns the number of calls to Initialize().
	int getNumInitialize() const
	{
		return(numInitialize_);
	}
	
	//! Returns the number of calls to Compute().
	int getNumCompute() const
	{
		return(numCompute_);
	}
	
	//! Returns the number of calls to ApplyInverse().
	int getNumApplyInverse() const
	{
		return(numApplyInverse_);
	}
	
	//! Returns the time spent in Initialize().
	double getInitializeTime() const
	{
		return(initializeTime_);
	}
	
	//! Returns the time spent in Compute().
	double getComputeTime() const
	{
		return(computeTime_);
	}
	
	//! Returns the time spent in ApplyInverse().
	double getApplyInverseTime() const
	{
		return(applyInverseTime_);
	}
	
	//! Returns the number of flops in the initialization phase.
	double getInitializeFlops() const
	{
		return(0.0);
	}
	
	double getComputeFlops() const
	{
		return(computeFlops_);
	}
	
	double getApplyInverseFlops() const
	{
		return(applyInverseFlops_);
	}
	
private:
	
	// @}
	// @{ Private methods
	
	//! Copy constructor (should never be used)
	ILU(const ILU& RHS){}
	
	//! operator= (should never be used)
	ILU& operator=(const ILU& RHS)
	{
		return(*this);
	}
	
	//! Destroys all internal data
	void destroy();
	
	//! Returns the result of a ILU forward/back solve on a Tpetra::MultiVector X in Y.
	/*! 
	 \param In
	 Trans -If true, solve transpose problem.
	 \param 
	 X - (In) A Tpetra::MultiVector of dimension NumVectors to solve for.
	 \param Out
	 Y - (Out) A Tpetra::MultiVector of dimension NumVectorscontaining result.
	 
	 \return Integer error code, set to 0 if successful.
	 */
	void solve(bool Trans, const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
			   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const;
	
	void computeSetup();
	void initAllValues(const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> & A, int MaxNumEntries);
	
	//! Returns the level of fill.
	int getLevelOfFill() const {return levelOfFill_;}
	
	//! Get ILU(k) relaxation parameter
	Scalar  getRelaxValue() const {return relaxValue_;}
	
	//! Get absolute threshold value
	Teuchos::ScalarTraits<Scalar>::magnitudeType  getAbsoluteThreshold() const {return athresh_;}
	
	//! Get relative threshold value
	Teuchos::ScalarTraits<Scalar>::magnitudeType  getRelativeThreshold() const {return rthresh_;}
	
	//! Returns the number of global matrix rows.
	GlobalOrdinal getNumGlobalRows() const {return(Graph().NumGlobalRows());};
	
	//! Returns the number of global matrix columns.
	GlobalOrdinal getNumGlobalCols() const {return(Graph().NumGlobalCols());};
	
	//! Returns the number of nonzero entries in the global graph.
	GlobalOrdinal getNumGlobalNonzeros() const {return(L().NumGlobalNonzeros()+U().NumGlobalNonzeros());};
	
	//! Returns the number of diagonal entries found in the global input graph.
	GlobalOrdinal getNumGlobalBlockDiagonals() const {return(Graph().NumGlobalBlockDiagonals());};
	
	//! Returns the number of local matrix rows.
	LocalOrdinal getNumMyRows() const {return(Graph().NumMyRows());};
	
	//! Returns the number of local matrix columns.
	LocalOrdinal getNumMyCols() const {return(Graph().NumMyCols());};
	
	//! Returns the number of nonzero entries in the local graph.
	LocalOrdinal getNumMyNonzeros() const {return(L().NumMyNonzeros()+U().NumMyNonzeros());};
	
	//! Returns the number of diagonal entries found in the local input graph.
	LocalOrdinal getNumMyBlockDiagonals() const {return(Graph().NumMyBlockDiagonals());};
	
	//! Returns the number of nonzero diagonal values found in matrix.
	LocalOrdinal getNumMyDiagonals() const {return(NumMyDiagonals_);};
	
	//! Returns the index base for row and column indices for this graph.
	int IndexBase() const {return(Graph().IndexBase());};
	
	//! Returns the address of the Ifpack2_IlukGraph associated with this factored matrix.
	const Ifpack2_IlukGraph & Graph() const {return(*Graph_);};
	
	//! Returns a reference to the matrix.
	Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Matrix()
	{
		return(*A_);
	}
	
	// @}
	// @{ Internal data
	
	//! Pointer to the Tpetra::RowMatrix to factorize
	Teuchos::RCP<Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>> A_;
	Teuchos::RCP<Ifpack2_IlukGraph<LocalOrdinal,GlobalOrdinal,Node>> Graph_;
	Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>> CrsGraph_;
	Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>> IlukRowMap_;
	Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>> IlukDomainMap_;
	Teuchos::RCP<Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node>> IlukRangeMap_;
	const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> * U_DomainMap_;
	const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> * L_RangeMap_;
	//! Contains the L factors
	Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>> L_;
	//! Contains the U factors.
	Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>> U_;
	Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>> L_Graph_;
	Teuchos::RCP<Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node>> U_Graph_;
	//! Diagonal of factors
	Teuchos::RCP<Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>> D_;
	bool useTranspose_;
	
	int numMyDiagonals_;
	bool allocated_;
	bool valuesInitialized_;
	bool factored_;
	//! Relaxation value
	double relaxValue_;
	//! absolute threshold
	double athresh_;
	//! relative threshold
	double rthresh_;
	//! condition number estimate
	double condest_;
	//! Level of fill
	int levelOfFill_;
	//! If \c true, the preconditioner has been successfully initialized.
	bool isInitialized_;
	//! If \c true, the preconditioner has been successfully computed.
	bool isComputed_;
	//! Label of \c this object.
	char label_[160];
	//! Contains the number of successful calls to Initialize().
	int numInitialize_;
	//! Contains the number of successful call to Compute().
	int numCompute_;
	//! Contains the number of successful call to ApplyInverse().
	mutable int numApplyInverse_;
	//! Contains the time for all successful calls to Initialize().
	double initializeTime_;
	//! Contains the time for all successful calls to Compute().
	double computeTime_;
	//! Contains the time for all successful calls to ApplyInverse().
	mutable double applyInverseTime_;
	//! Contains the number of flops for Compute().
	double computeFlops_;
	//! Contain sthe number of flops for ApplyInverse().
	mutable double applyInverseFlops_;
	
};


////////////////////////////////////////
////////////////////////////////////////

//==============================================================================
ILU::ILU(Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>* Matrix_in) :
A_(rcp(Matrix_in,false)),
useTranspose_(false),
numMyDiagonals_(Teuchos::OrdinalTraits<LocalOrdinal>::zero()),
relaxValue_(Teuchos::ScalarTraits<Scalar>::zero()),
athresh_(Teuchos::ScalarTraits<Teuchos::ScalarTraits<Scalar>::magnitudeType>::zero()),
rthresh_(Teuchos::ScalarTraits<Teuchos::ScalarTraits<Scalar>::magnitudeType>::one()),
condest_(-Teuchos::ScalarTraits<Teuchos::ScalarTraits<Scalar>::magnitudeType>::one()),
levelOfFill_(0),
isInitialized_(false),
isComputed_(false),
numInitialize_(0),
numCompute_(0),
numApplyInverse_(0),
initializeTime_(0.0),
computeTime_(0.0),
applyInverseTime_(0.0),
computeFlops_(0.0),
applyInverseFlops_(0.0)
{
	Teuchos::ParameterList List;
	SetParameters(List);
}

//==============================================================================
void ILU::destroy()
{
	// reset pointers to already allocated stuff
	U_DomainMap_ = 0;
	L_RangeMap_ = 0;
}

//==========================================================================
int ILU::setParameters(const Teuchos::ParameterList& list)
{
	relaxValue_ = list.get("fact: relax value", RelaxValue_);
	athresh_ = list.get("fact: absolute threshold", Athresh_);
	rthresh_ = list.get("fact: relative threshold", Rthresh_);
	levelOfFill_ = list.get("fact: level-of-fill", LevelOfFill_);
	
	// set label
	sprintf(label_, "TIFPACK ILU (fill=%d, relax=%f, athr=%f, rthr=%f)",
			levelOfFill(), relaxValue(), absoluteThreshold(), 
			relativeThreshold());
	return(0);
}

//==========================================================================
void ILU::computeSetup() {
	
#ifdef ENABLE_IFPACK2_ILU_TEUCHOS_TIMERS
	TEUCHOS_FUNC_TIME_MONITOR("ILU::computeSetup");
#endif
	typedef Teuchos::ScalarTraits<Scalar> ST;
	typedef Teuchos::OrdinalTraits<LocalOrdinal> LOT;
	LocalOrdinal izero = LOT::zero();
	Scalar zero = ST::zero();
	
	L_ = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Copy, Graph().getL_Graph()));
	U_ = rcp(new Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Copy, Graph().getU_Graph()));
	D_ = rcp(new Tpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(Graph().getL_Graph().getRowMap()));
	TEUCHOS_TEST_FOR_EXCEPTION(((L_.get() == 0) || (U_.get() == 0) || (D_.get() == 0)), std::runtime_error,
					   "Ifpack2::ILU::computeSetup(): Local column value of "<< k << 
					   " for user matrix at local row " << i << " and index " << j << " less than zero on node " <<
					   getGraph().getL_Graph().getRowMap().getComm()->getRank());
	
	// Get Maximun Row length
	LocalOrdinal maxNumEntries = getMatrix().getMaxNumEntries();
	
	// Set L range map and U domain map
	U_DomainMap_ = &(Matrix().OperatorDomainMap());
	L_RangeMap_ = &(Matrix().OperatorRangeMap());
	
	// this is the old InitAllValues()
	LocalOrdinal i, j;
	LocalOrdinal numIn, numL, numU;
	bool diagFound;
	LocalOrdinal numNonzeroDiags = izero;
	
	vector<LocalOrdinal> InI(MaxNumEntries); // Allocate temp space
	vector<LocalOrdinal> LI(MaxNumEntries);
	vector<LocalOrdinal> UI(MaxNumEntries);
	vector<Scalar> InV(MaxNumEntries);
	vector<Scalar> LV(MaxNumEntries);
	vector<Scalar> UV(MaxNumEntries);
	
	bool ReplaceValues = (L_->getStaticGraph() || L_->getIndicesAreLocal()); // Check if values should be inserted or replaced
	
	if (ReplaceValues) {
		L_->putScalar(zero); // Zero out L and U matrices
		U_->putScalar(zero);
	}
	
	D_->putScalar(zero); // Set diagonal values to zero
	Scalar * DV;
	D_->extractView(&DV); // Get view of diagonal
	
	// First we copy the user's matrix into L and U, regardless of fill level
	
	for (i=0; i< getNumMyRows(); i++) {
		
		getMatrix().extractMyRowCopy(i, maxNumEntries, numIn, &InV[0], &InI[0]); // Get Values and Indices
		
		// Split into L and U (we don't assume that indices are ordered).
		
		numL = izero; 
		numU = izero; 
		diagFound = false;
		
		for (j=0; j< numIn; j++) {
			int k = InI[j];
			
			TEUCHOS_TEST_FOR_EXCEPTION((k<0 || k> getNumMyRows()), std::runtime_error,
							   "Ifpack2::ILU::computeSetup(): Local column value of "<< k 
							   << " for user matrix at local row " << i << " and index " << j << " out of column range " 
							   <<  izero << " and " << getNumMyRows() << " on node " <<
							   getGraph().getL_Graph().getRowMap().getComm()->getRank());
			if (k==i) {
				diagFound = true;
				if (InV[j]>=zero)
					DV[i] += Rthresh_ * InV[j] + Athresh_; // Store perturbed diagonal in Tpetra::Vector D_
				else
					DV[i] += Rthresh_ * InV[j] - Athresh_; // Store perturbed diagonal in Tpetra::Vector D_
				
			}
			
			
			else if (k < i) {
				LI[NumL] = k;
				LV[NumL] = InV[j];
				NumL++;
			}
			else { // if (k<NumMyRows())
				UI[NumU] = k;
				UV[NumU] = InV[j];
				NumU++;
			}
		}
		
		// Check in things for this row of L and U
		
		if (DiagFound) numNonzeroDiags++;
		else DV[i] = athresh_;
		
		if (numL) {
			if (replaceValues)
				L_->replaceMyValues(i, numL, &LV[0], &LI[0]);
			else
				L_->insertMyValues(i, numL, &LV[0], &LI[0]);
		}
		
		if (numU) {
			if (replaceValues)
				U_->replaceMyValues(i, numU, &UV[0], &UI[0]);
			else
				U_->insertMyValues(i, numU, &UV[0], &UI[0]);
		}
		
	}
	
	if (!replaceValues) {
		// The domain of L and the range of U are exactly their own row maps (there is no communication).
		// The domain of U and the range of L must be the same as those of the original matrix,
		// However if the original matrix is a VbrMatrix, these two latter maps are translation from
		// a block map to a point map.
		L_->fillComplete((L_->getRowMatrixColMap()), *L_RangeMap_);
		U_->fillComplete(*U_DomainMap_, U_->getRowMatrixRowMap());
	}
	
	// At this point L and U have the values of A in the structure of L and U, and diagonal vector D
	GlobalOrdinal totalNonzeroDiags = Teuchos::OrdinalTraits<GlobalOrdinal>::zero();
	GlobalOrdinal temp = (GlobalOrdinal)numNonzeroDiags;
	getGraph().getL_Graph().getRowMap().getComm().SumAll(&numNonzeroDiags, &totalNonzeroDiags, 1));
	numMyDiagonals_ = numNonzeroDiags;
	
	TEUCHOS_TEST_FOR_EXCEPTION(numNonzeroDiags != getNumMyRows(), std::runtime_error,
					   "Ifpack2::ILU::computeSetup(): Number of nonzero diagonal values "<< numNonzeroDiags << 
					   " for factorization is less than the number of rows " << getNumMyRows() << " in factor node " <<
					   getGraph().getL_Graph().getRowMap().getComm()->getRank());	 // Diagonals are not right, warn user
}

//==========================================================================
int ILU::Initialize() 
{
	
#ifdef ENABLE_IFPACK2_ILU_TEUCHOS_TIMERS
	TEUCHOS_FUNC_TIME_MONITOR("ILU::Initialize");
#endif
	
	Time_.ResetStartTime();
	IsInitialized_ = false;
	
	// reset this object
	Destroy();
	
	Tpetra::CrsMatrix* CrsMatrix;
	CrsMatrix = dynamic_cast<Tpetra::CrsMatrix*>(&*A_);
	if (CrsMatrix == 0) {
		// this means that we have to create
		// the graph from a given Tpetra::RowMatrix. Note
		// that at this point we are ignoring any possible
		// graph coming from VBR matrices.
		int size = A_->MaxNumEntries();
		CrsGraph_ = rcp(new Tpetra::CrsGraph(Copy,A_->RowMatrixRowMap(), size));
		if (CrsGraph_.get() == 0)
			IFPACK2_CHK_ERR(-5); // memory allocation error
		
		vector<int> Indices(size);
		vector<double> Values(size);
		
		// extract each row at-a-time, and insert it into
		// the graph, ignore all off-process entries
		for (int i = 0 ; i < A_->NumMyRows() ; ++i) {
			int NumEntries;
			int GlobalRow = A_->RowMatrixRowMap().GID(i);
			IFPACK2_CHK_ERR(A_->ExtractMyRowCopy(i, size, NumEntries, 
												 &Values[0], &Indices[0]));
			// convert to global indices
			for (int j = 0 ; j < NumEntries ; ++j) {
				Indices[j] = A_->RowMatrixColMap().GID(Indices[j]); 
			}
			IFPACK2_CHK_ERR(CrsGraph_->InsertGlobalIndices(GlobalRow,NumEntries,
														   &Indices[0]));
		}
		
		IFPACK2_CHK_ERR(CrsGraph_->FillComplete(A_->RowMatrixRowMap(),
												A_->RowMatrixRowMap()));
		
		// always overlap zero, wider overlap will be handled
		// by the AdditiveSchwarz preconditioner.
		Graph_ = rcp(new Ifpack2_IlukGraph(*CrsGraph_, LevelOfFill_, 0));
		
	}
	else {
		// see comment above for the overlap.
		Graph_ = rcp(new Ifpack2_IlukGraph(CrsMatrix->Graph(), LevelOfFill_, 0));
	}
	
	if (Graph_.get() == 0)
		IFPACK2_CHK_ERR(-5); // memory allocation error
	IFPACK2_CHK_ERR(Graph_->ConstructFilledGraph());
	
	IsInitialized_ = true;
	NumInitialize_++;
	InitializeTime_ += Time_.ElapsedTime();
	
	return(0);
}

//==========================================================================
int ILU::Compute() 
{
	
#ifdef ENABLE_IFPACK2_ILU_TEUCHOS_TIMERS
	TEUCHOS_FUNC_TIME_MONITOR("ILU::Compute");
#endif
	
	if (!IsInitialized()) 
		IFPACK2_CHK_ERR(Initialize());
	
	Time_.ResetStartTime();
	IsComputed_ = false;
	
	// convert Matrix() into L and U factors.
	IFPACK2_CHK_ERR(ComputeSetup());
	
	// MinMachNum should be officially defined, for now pick something a little 
	// bigger than IEEE underflow value
	
	double MinDiagonalValue = Tpetra::MinDouble;
	double MaxDiagonalValue = 1.0/MinDiagonalValue;
	
	int ierr = 0;
	int i, j, k;
	int *LI, *UI;
	double *LV, *UV;
	int NumIn, NumL, NumU;
	
	// Get Maximun Row length
	int MaxNumEntries = L_->MaxNumEntries() + U_->MaxNumEntries() + 1;
	
	vector<int> InI(MaxNumEntries+1);    // Allocate temp space, pad by one to 
	vector<double> InV(MaxNumEntries+1); // to avoid debugger complaints for pathological cases
	vector<int> colflag(NumMyCols());
	
	double *DV;
	ierr = D_->ExtractView(&DV); // Get view of diagonal
	
	int current_madds = 0; // We will count multiply-add as they happen
	
	// =========================== //
	// Now start the factorization //
	// =========================== //
	
	// Need some integer workspace and pointers
	int NumUU; 
	int * UUI;
	double * UUV;
	for (j = 0; j < NumMyCols(); ++j) colflag[j] = - 1;
	
	for (i = 0; i < NumMyRows(); ++i) {
		
		// Fill InV, InI with current row of L, D and U combined
		
		NumIn = MaxNumEntries;
		IFPACK2_CHK_ERR(L_->ExtractMyRowCopy(i, NumIn, NumL, &InV[0], &InI[0]));
		LV = &InV[0];
		LI = &InI[0];
		
		InV[NumL] = DV[i]; // Put in diagonal
		InI[NumL] = i;
		
		IFPACK2_CHK_ERR(U_->ExtractMyRowCopy(i, NumIn-NumL-1, NumU, &InV[NumL+1], &InI[NumL+1]));
		NumIn = NumL+NumU+1;
		UV = &InV[NumL+1];
		UI = &InI[NumL+1];
		
		// Set column flags
		for (j=0; j<NumIn; j++) colflag[InI[j]] = j;
		
		double diagmod = 0.0; // Off-diagonal accumulator
		
		for (int jj=0; jj<NumL; jj++) {
			j = InI[jj];
			double multiplier = InV[jj]; // current_mults++;
			
			InV[jj] *= DV[j];
			
			IFPACK2_CHK_ERR(U_->ExtractMyRowView(j, NumUU, UUV, UUI)); // View of row above
			
			if (RelaxValue_==0.0) {
				for (k=0; k<NumUU; k++) {
					int kk = colflag[UUI[k]];
					if (kk>-1) {
						InV[kk] -= multiplier*UUV[k];
						current_madds++;
					}
				}
			}
			else {
				for (k=0; k<NumUU; k++) {
					int kk = colflag[UUI[k]];
					if (kk>-1) InV[kk] -= multiplier*UUV[k];
					else diagmod -= multiplier*UUV[k];
					current_madds++;
				}
			}
		}
		if (NumL) {
			IFPACK2_CHK_ERR(L_->ReplaceMyValues(i, NumL, LV, LI));  // Replace current row of L
		}
		
		DV[i] = InV[NumL]; // Extract Diagonal value
		
		if (RelaxValue_!=0.0) {
			DV[i] += RelaxValue_*diagmod; // Add off diagonal modifications
			// current_madds++;
		}
		
		if (fabs(DV[i]) > MaxDiagonalValue) {
			if (DV[i] < 0) DV[i] = - MinDiagonalValue;
			else DV[i] = MinDiagonalValue;
		}
		else
			DV[i] = 1.0/DV[i]; // Invert diagonal value
		
		for (j=0; j<NumU; j++) UV[j] *= DV[i]; // Scale U by inverse of diagonal
		
		if (NumU) {
			IFPACK2_CHK_ERR(U_->ReplaceMyValues(i, NumU, UV, UI));  // Replace current row of L and U
		}
		
		// Reset column flags
		for (j=0; j<NumIn; j++) colflag[InI[j]] = -1;
	}
	
	// Validate that the L and U factors are actually lower and upper triangular
	
	if (!L_->LowerTriangular()) 
		IFPACK2_CHK_ERR(-4);
	if (!U_->UpperTriangular()) 
		IFPACK2_CHK_ERR(-4);
	
	// Add up flops
	
	double current_flops = 2 * current_madds;
	double total_flops = 0;
	
	IFPACK2_CHK_ERR(Graph().L_Graph().RowMap().Comm().SumAll(&current_flops, &total_flops, 1)); // Get total madds across all PEs
	
	// Now count the rest
	total_flops += (double) L_->NumGlobalNonzeros(); // Accounts for multiplier above
	total_flops += (double) D_->GlobalLength(); // Accounts for reciprocal of diagonal
	if (RelaxValue_!=0.0) total_flops += 2 * (double)D_->GlobalLength(); // Accounts for relax update of diag
	
	// add to this object's counter
	ComputeFlops_ += total_flops;
	
	IsComputed_ = true;
	NumCompute_++;
	ComputeTime_ += Time_.ElapsedTime();
	
	return(ierr);
	
}

//=============================================================================
// This function finds Y such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int ILU::Solve(bool Trans, const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
			   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const 
{
	
#ifdef ENABLE_IFPACK2_ILU_TEUCHOS_TIMERS
	TEUCHOS_FUNC_TIME_MONITOR("ILU::ApplyInverse - Solve");
#endif
	
	// in this function the overlap is always zero
	bool Upper = true;
	bool Lower = false;
	bool UnitDiagonal = true;
	
	if (!Trans) {
		
		IFPACK2_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, X, Y));
		// y = D*y (D_ has inverse of diagonal)
		IFPACK2_CHK_ERR(Y.Multiply(1.0, *D_, Y, 0.0)); 
		// Solve Uy = y
		IFPACK2_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, Y, Y)); 
	}
	else {
		// Solve Uy = y
		IFPACK2_CHK_ERR(U_->Solve(Upper, Trans, UnitDiagonal, X, Y)); 
		// y = D*y (D_ has inverse of diagonal)
		IFPACK2_CHK_ERR(Y.Multiply(1.0, *D_, Y, 0.0)); 
		IFPACK2_CHK_ERR(L_->Solve(Lower, Trans, UnitDiagonal, Y, Y));
	} 
	
	
	return(0);
}

//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int ILU::Multiply(bool Trans, const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
				  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const 
{
	
#ifdef ENABLE_IFPACK2_ILU_TEUCHOS_TIMERS
	TEUCHOS_FUNC_TIME_MONITOR("ILU::Multiply");
#endif
	
	if (!IsComputed())
		IFPACK2_CHK_ERR(-3);
	
	if (!Trans) {
		IFPACK2_CHK_ERR(U_->Multiply(Trans, X, Y)); 
		// Y1 = Y1 + X1 (account for implicit unit diagonal)
		IFPACK2_CHK_ERR(Y.Update(1.0, X, 1.0)); 
		// y = D*y (D_ has inverse of diagonal)
		IFPACK2_CHK_ERR(Y.ReciprocalMultiply(1.0, *D_, Y, 0.0)); 
		Tpetra::MultiVector Y1temp(Y); // Need a temp copy of Y1
		IFPACK2_CHK_ERR(L_->Multiply(Trans, Y1temp, Y));
		// (account for implicit unit diagonal)
		IFPACK2_CHK_ERR(Y.Update(1.0, Y1temp, 1.0)); 
	}
	else {
		
		IFPACK2_CHK_ERR(L_->Multiply(Trans, X, Y));
		// Y1 = Y1 + X1 (account for implicit unit diagonal)
		IFPACK2_CHK_ERR(Y.Update(1.0, X, 1.0)); 
		// y = D*y (D_ has inverse of diagonal)
		IFPACK2_CHK_ERR(Y.ReciprocalMultiply(1.0, *D_, Y, 0.0)); 
		Tpetra::MultiVector Y1temp(Y); // Need a temp copy of Y1
		IFPACK2_CHK_ERR(U_->Multiply(Trans, Y1temp, Y));
		// (account for implicit unit diagonal)
		IFPACK2_CHK_ERR(Y.Update(1.0, Y1temp, 1.0)); 
	} 
	
	return(0);
}

//=============================================================================
// This function finds X such that LDU Y = X or U(trans) D L(trans) Y = X for multiple RHS
int ILU::ApplyInverse(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X, 
					  Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y) const
{
	
#ifdef ENABLE_IFPACK2_ILU_TEUCHOS_TIMERS
	TEUCHOS_FUNC_TIME_MONITOR("ILU::ApplyInverse");
#endif
	
	if (!IsComputed())
		IFPACK2_CHK_ERR(-3);
	
	if (X.NumVectors() != Y.NumVectors())
		IFPACK2_CHK_ERR(-2);
	
	Time_.ResetStartTime();
	
	// AztecOO gives X and Y pointing to the same memory location,
	// need to create an auxiliary vector, Xcopy
	Teuchos::RCP< const Tpetra::MultiVector > Xcopy;
	if (X.Pointers()[0] == Y.Pointers()[0])
		Xcopy = Teuchos::rcp( new Tpetra::MultiVector(X) );
	else
		Xcopy = Teuchos::rcp( &X, false );
	
	IFPACK2_CHK_ERR(Solve(ILU::UseTranspose(), *Xcopy, Y));
	
	// approx is the number of nonzeros in L and U
	ApplyInverseFlops_ += X.NumVectors() * 4 * 
	(L_->NumGlobalNonzeros() + U_->NumGlobalNonzeros());
	
	++NumApplyInverse_;
	ApplyInverseTime_ += Time_.ElapsedTime();
	
	return(0);
	
}

//=============================================================================
double ILU::Condest(const Ifpack2_CondestType CT, 
					const int MaxIters, const double Tol,
					Tpetra::RowMatrix* Matrix_in)
{
	
#ifdef ENABLE_IFPACK2_ILU_TEUCHOS_TIMERS
	TEUCHOS_FUNC_TIME_MONITOR("ILU::Condest");
#endif
	
	if (!IsComputed()) // cannot compute right now
		return(-1.0);
	
	Condest_ = Ifpack2_Condest(*this, CT, MaxIters, Tol, Matrix_in);
	
	return(Condest_);
}

//=============================================================================
std::ostream&
ILU::Print(std::ostream& os) const
{
	if (!Comm().MyPID()) {
		os << endl;
		os << "================================================================================" << endl;
		os << "ILU: " << Label() << endl << endl;
		os << "Level-of-fill      = " << LevelOfFill() << endl;
		os << "Absolute threshold = " << AbsoluteThreshold() << endl;
		os << "Relative threshold = " << RelativeThreshold() << endl;
		os << "Relax value        = " << RelaxValue() << endl;
		os << "Condition number estimate = " << Condest() << endl;
		os << "Global number of rows            = " << A_->NumGlobalRows() << endl;
		if (IsComputed_) {
			os << "Number of rows of L, D, U       = " << L_->NumGlobalRows() << endl;
			os << "Number of nonzeros of L + U     = " << NumGlobalNonzeros() << endl;
			os << "nonzeros / rows                 = " 
			<< 1.0 * NumGlobalNonzeros() / U_->NumGlobalRows() << endl;
		}
		os << endl;
		os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
		os << "-----           -------   --------------       ------------     --------" << endl;
		os << "Initialize()    "   << std::setw(5) << NumInitialize() 
		<< "  " << std::setw(15) << InitializeTime() 
		<< "               0.0            0.0" << endl;
		os << "Compute()       "   << std::setw(5) << NumCompute() 
		<< "  " << std::setw(15) << ComputeTime()
		<< "  " << std::setw(15) << 1.0e-6 * ComputeFlops();
		if (ComputeTime() != 0.0) 
			os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops() / ComputeTime() << endl;
		else
			os << "  " << std::setw(15) << 0.0 << endl;
		os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse() 
		<< "  " << std::setw(15) << ApplyInverseTime()
		<< "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops();
		if (ApplyInverseTime() != 0.0) 
			os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops() / ApplyInverseTime() << endl;
		else
			os << "  " << std::setw(15) << 0.0 << endl;
		os << "================================================================================" << endl;
		os << endl;
	}
	
	return(os);
}
} // namespace Ifpack2
#endif /* IFPACK2_ILU_HPP */
