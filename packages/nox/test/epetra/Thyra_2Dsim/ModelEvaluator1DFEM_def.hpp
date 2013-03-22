#ifndef NOX_THYRA_MODEL_EVALUATOR_1DFEM_DEF_HPP
#define NOX_THYRA_MODEL_EVALUATOR_1DFEM_DEF_HPP

// Thyra support
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_PreconditionerBase.hpp"

// Epetra support
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

// Nonmember constuctors

template<class Scalar>
Teuchos::RCP<ModelEvaluator1DFEM<Scalar> >
modelEvaluator1DFEM(const Teuchos::RCP<const Epetra_Comm>& comm,
		    const int num_global_elements,
		    const Scalar z_min,
		    const Scalar z_max)
{
  return Teuchos::rcp(new ModelEvaluator1DFEM<Scalar>(comm,num_global_elements,z_min,z_max));
}

// Constructor

template<class Scalar>
ModelEvaluator1DFEM<Scalar>::
ModelEvaluator1DFEM(const Teuchos::RCP<const Epetra_Comm>& comm,
		    const int num_global_elements,
		    const Scalar z_min,
		    const Scalar z_max) :
  comm_(comm),
  num_global_elements_(num_global_elements),
  z_min_(z_min),
  z_max_(z_max),
  showGetInvalidArg_(false)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using ::Thyra::VectorBase;
  typedef ::Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  
  TEUCHOS_ASSERT(nonnull(comm_));

  const int num_nodes = num_global_elements_ + 1;

  // owned space
  x_owned_map_ = rcp(new Epetra_Map(num_nodes,0,*comm_));
  x_space_ = ::Thyra::create_VectorSpace(x_owned_map_);

  // ghosted space
  if (comm_->NumProc() == 1) {
    x_ghosted_map_ = x_owned_map_;
  } else {

    int OverlapNumMyElements;
    int OverlapMinMyGID;
    OverlapNumMyElements = x_owned_map_->NumMyElements() + 2;
    if ( (comm_->MyPID() == 0) || (comm_->MyPID() == (comm_->NumProc() - 1)) ) 
      OverlapNumMyElements --;
    
    if (comm_->MyPID() == 0) 
      OverlapMinMyGID = x_owned_map_->MinMyGID();
    else 
      OverlapMinMyGID = x_owned_map_->MinMyGID() - 1;
    
    int* OverlapMyGlobalElements = new int[OverlapNumMyElements];
    
    for (int i = 0; i < OverlapNumMyElements; i ++) 
      OverlapMyGlobalElements[i] = OverlapMinMyGID + i;
    
    x_ghosted_map_ = 
      Teuchos::rcp(new Epetra_Map(-1, 
				  OverlapNumMyElements, 
				  OverlapMyGlobalElements,
				  0,
				  *comm_));

    delete [] OverlapMyGlobalElements;

  }

  importer_ = Teuchos::rcp(new Epetra_Import(*x_ghosted_map_, *x_owned_map_));

  // residual space
  f_owned_map_ = x_owned_map_;
  f_space_ = x_space_;

  x0_ = ::Thyra::createMember(x_space_);
  V_S(x0_.ptr(), ST::zero());

//   set_p(Teuchos::tuple<Scalar>(p0, p1)());
//   set_x0(Teuchos::tuple<Scalar>(x0, x1)());

  // Initialize the graph for W CrsMatrix object
  W_graph_ = createGraph();

  // Create the nodal coordinates
  {
    node_coordinates_ = Teuchos::rcp(new Epetra_Vector(*x_owned_map_));
    Scalar length = z_max_ - z_min_;
    Scalar dx = length/((double) num_global_elements_ - 1);
    for (int i=0; i < x_owned_map_->NumMyElements(); i++) {
      (*node_coordinates_)[i] = z_min_ + dx*((double) x_owned_map_->MinMyGID() + i);
    }
  }
  
  

  MEB::InArgsSetup<Scalar> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(MEB::IN_ARG_x);
  prototypeInArgs_ = inArgs;
  
  MEB::OutArgsSetup<Scalar> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(MEB::OUT_ARG_f);
  outArgs.setSupports(MEB::OUT_ARG_W_op);
  outArgs.setSupports(MEB::OUT_ARG_W_prec);
//   outArgs.set_W_properties(DerivativeProperties(
// 			     DERIV_LINEARITY_NONCONST
// 			     ,DERIV_RANK_FULL
// 			     ,true // supportsAdjoint
// 			     ));
  prototypeOutArgs_ = outArgs;

  nominalValues_ = inArgs;
  nominalValues_.set_x(x0_);
}

// Initializers/Accessors

template<class Scalar>
Teuchos::RCP<Epetra_CrsGraph> 
ModelEvaluator1DFEM<Scalar>::createGraph()
{
  Teuchos::RCP<Epetra_CrsGraph> W_graph;

  // Create the shell for the 
  W_graph = Teuchos::rcp(new Epetra_CrsGraph(Copy, *x_owned_map_, 5));

  // Declare required variables
  int row;
  int column;
  int OverlapNumMyElements = x_ghosted_map_->NumMyElements();
  
  // Loop Over # of Finite Elements on Processor
  for (int ne=0; ne < OverlapNumMyElements-1; ne++) {
          
    // Loop over Nodes in Element
    for (int i=0; i< 2; i++) {
      row=x_ghosted_map_->GID(ne+i);
      
      // Loop over Trial Functions
      for(int j=0;j < 2; j++) {
	
	// If this row is owned by current processor, add the index
	if (x_owned_map_->MyGID(row)) {
	  column=x_ghosted_map_->GID(ne+j);
	  W_graph->InsertGlobalIndices(row, 1, &column);
	}
      } 	
    }
  }
  W_graph->FillComplete();
  return W_graph;
}

template<class Scalar>
void ModelEvaluator1DFEM<Scalar>::set_x0(const Teuchos::ArrayView<const Scalar> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(x_space_->dim(), x0_in.size());
#endif
  Thyra::DetachedVectorView<Scalar> x0(x0_);
  x0.sv().values()().assign(x0_in);
}


template<class Scalar>
void ModelEvaluator1DFEM<Scalar>::setShowGetInvalidArgs(bool showGetInvalidArg)
{
  showGetInvalidArg_ = showGetInvalidArg;
}

template<class Scalar>
void ModelEvaluator1DFEM<Scalar>::
set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar> >& W_factory)
{
  W_factory_ = W_factory;
}

// Public functions overridden from ModelEvaulator


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
ModelEvaluator1DFEM<Scalar>::get_x_space() const
{
  return x_space_;
}


template<class Scalar>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> >
ModelEvaluator1DFEM<Scalar>::get_f_space() const
{
  return f_space_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluator1DFEM<Scalar>::getNominalValues() const
{
  return nominalValues_;
}


template<class Scalar>
Teuchos::RCP<Thyra::LinearOpBase<Scalar> >
ModelEvaluator1DFEM<Scalar>::create_W_op() const
{
  Teuchos::RCP<Epetra_CrsMatrix> W_epetra =
    Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));

  return Thyra::nonconstEpetraLinearOp(W_epetra);
}

template<class Scalar>
Teuchos::RCP< ::Thyra::PreconditionerBase<Scalar> >
ModelEvaluator1DFEM<Scalar>::create_W_prec() const
{
  Teuchos::RCP<Epetra_CrsMatrix> W_epetra =
    Teuchos::rcp(new Epetra_CrsMatrix(::Copy,*W_graph_));

  const Teuchos::RCP<Thyra::LinearOpBase< Scalar > > W_op = 
    Thyra::nonconstEpetraLinearOp(W_epetra);

  Teuchos::RCP<Thyra::DefaultPreconditioner<Scalar> > prec = 
    Teuchos::rcp(new Thyra::DefaultPreconditioner<Scalar>);

  prec->initializeRight(W_op);
  
  return prec;
}

template<class Scalar>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
ModelEvaluator1DFEM<Scalar>::get_W_factory() const
{
  return W_factory_;
}


template<class Scalar>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ModelEvaluator1DFEM<Scalar>::createInArgs() const
{
  return prototypeInArgs_;
}


// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ModelEvaluator1DFEM<Scalar>::createOutArgsImpl() const
{
  return prototypeOutArgs_;
}


template<class Scalar>
void ModelEvaluator1DFEM<Scalar>::evalModelImpl(
  const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
  const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
  ) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  
  TEUCHOS_ASSERT(nonnull(inArgs.get_x()));
		 
  //const Thyra::ConstDetachedVectorView<Scalar> x(inArgs.get_x());

  const RCP<Thyra::VectorBase<Scalar> > f_out = outArgs.get_f();
  const RCP<Thyra::LinearOpBase<Scalar> > W_out = outArgs.get_W_op();
  const RCP<Thyra::PreconditionerBase<Scalar> > W_prec_out = outArgs.get_W_prec();


  if ( nonnull(f_out) || nonnull(W_out) || nonnull(W_prec_out) ) {

    // ****************
    // Get the underlying epetra objects
    // ****************

    RCP<Epetra_Vector> f;
    if (nonnull(f_out)) {
      f = Thyra::get_Epetra_Vector(*f_owned_map_,outArgs.get_f());
    }

    RCP<Epetra_CrsMatrix> J;
    if (nonnull(W_out)) {
      RCP<Epetra_Operator> W_epetra = Thyra::get_Epetra_Operator(*W_out);
      J = rcp_dynamic_cast<Epetra_CrsMatrix>(W_epetra);
      TEUCHOS_ASSERT(nonnull(J));
    }

    RCP<Epetra_CrsMatrix> M_inv;
    if (nonnull(W_prec_out)) {
      RCP<Epetra_Operator> M_epetra = Thyra::get_Epetra_Operator(*(W_prec_out->getNonconstRightPrecOp()));
      M_inv = rcp_dynamic_cast<Epetra_CrsMatrix>(M_epetra);
      TEUCHOS_ASSERT(nonnull(M_inv));
      J_diagonal_ = Teuchos::rcp(new Epetra_Vector(*x_owned_map_));
    }

    // ****************
    // Create ghosted objects
    // ****************

    if (is_null(u_ptr))
      u_ptr = Teuchos::rcp(new Epetra_Vector(*x_ghosted_map_));

    u_ptr->Import(*(Thyra::get_Epetra_Vector(*x_owned_map_,inArgs.get_x())), *importer_, Insert);
    
    if (is_null(x_ptr)) {
      x_ptr = Teuchos::rcp(new Epetra_Vector(*x_ghosted_map_));
      x_ptr->Import(*node_coordinates_, *importer_, Insert);
    }

    Epetra_Vector& u = *u_ptr;
    Epetra_Vector& x = *x_ptr;

    double factor = 1.0;

    int ierr = 0;
    int OverlapNumMyElements = x_ghosted_map_->NumMyElements();
    
    double xx[2];
    double uu[2];
    Basis basis;
    
    // Zero out the objects that will be filled
    if (nonnull(f)) 
      f->PutScalar(0.0);
    if (nonnull(J)) 
      J->PutScalar(0.0);
    if (nonnull(M_inv)) 
      M_inv->PutScalar(0.0);
    
    // Loop Over # of Finite Elements on Processor
    for (int ne=0; ne < OverlapNumMyElements-1; ne++) {
      
      // Loop Over Gauss Points
      for(int gp=0; gp < 2; gp++) {
	// Get the solution and coordinates at the nodes 
	xx[0]=x[ne];
	xx[1]=x[ne+1];
	uu[0]=u[ne];
	uu[1]=u[ne+1];
	// Calculate the basis function at the gauss point
	basis.computeBasis(gp, xx, uu);
	
	// Loop over Nodes in Element
	for (int i=0; i< 2; i++) {
	  int row=x_ghosted_map_->GID(ne+i);
	  //printf("Proc=%d GlobalRow=%d LocalRow=%d Owned=%d\n",
	  //     MyPID, row, ne+i,x_owned_map_.MyGID(row));
	  if (x_owned_map_->MyGID(row)) {
	    if (nonnull(f)) {
	      (*f)[x_owned_map_->LID(x_ghosted_map_->GID(ne+i))]+=
		+basis.wt*basis.dz
		*((1.0/(basis.dz*basis.dz))*basis.duu*
		  basis.dphide[i]+factor*basis.uu*basis.uu*basis.phi[i]);
	    }
	  }
	  // Loop over Trial Functions
	  if (nonnull(J)) {
	    for(int j=0;j < 2; j++) {
	      if (x_owned_map_->MyGID(row)) {
		int column=x_ghosted_map_->GID(ne+j);
		double jac=basis.wt*basis.dz*((1.0/(basis.dz*basis.dz))*
				       basis.dphide[j]*basis.dphide[i]
				       +2.0*factor*basis.uu*basis.phi[j]*
				       basis.phi[i]);  
		ierr=J->SumIntoGlobalValues(row, 1, &jac, &column);
	      }
	    }
	  }
	  if (nonnull(M_inv)) {
	    for(int j=0;j < 2; j++) {
	      if (x_owned_map_->MyGID(row)) {
		int column=x_ghosted_map_->GID(ne+j);
		if (row == column) {
		  double jac = basis.wt*basis.dz*((1.0/(basis.dz*basis.dz))*
						  basis.dphide[j]*basis.dphide[i]
						  +2.0*factor*basis.uu*basis.phi[j]*
						 basis.phi[i]);
		  ierr = M_inv->SumIntoGlobalValues(row, 1, &jac, &column);
		}
	      }
	    }
	  }
	}
      }
    } 
    
    // Insert Boundary Conditions and modify Jacobian and function (F)
    // U(0)=1
    if (comm_->MyPID() == 0) {
      if (nonnull(f)) 
	(*f)[0]= u[0] - 1.0;
      if (nonnull(J)) {
	int column=0;
	double jac=1.0;
	ierr = J->ReplaceGlobalValues(0, 1, &jac, &column);
	column=1;
	jac=0.0;
	ierr = J->ReplaceGlobalValues(0, 1, &jac, &column);
      }
      if (nonnull(M_inv)) {
	int column=0;
	double jac=1.0;
	ierr = M_inv->ReplaceGlobalValues(0, 1, &jac, &column);
	column=1;
	jac=0.0;
	ierr = M_inv->ReplaceGlobalValues(0, 1, &jac, &column);
      }
    }

    if (nonnull(J))
      J->FillComplete();
    
    if (nonnull(M_inv)) {

      // invert the Jacobian diagonal for the preconditioner    
      M_inv->ExtractDiagonalCopy(*J_diagonal_);
      
      for (int i=0; i < J_diagonal_->MyLength(); ++i)
	(*J_diagonal_)[i] = 1.0 / ((*J_diagonal_)[i]);

      M_inv->ReplaceDiagonalValues(*J_diagonal_);

      M_inv->FillComplete();
    }

    TEUCHOS_ASSERT(ierr > -1);

  }

}

//====================================================================
// Basis vector

// Constructor
Basis::Basis() {
  phi = new double[2];
  dphide = new double[2];
}

// Destructor
Basis::~Basis() {
  delete [] phi;
  delete [] dphide;
}

// Calculates a linear 1D basis
void Basis::computeBasis(int gp, double *z, double *u, double *uold) {
  int N = 2;
  if (gp==0) {eta=-1.0/sqrt(3.0); wt=1.0;}
  if (gp==1) {eta=1.0/sqrt(3.0); wt=1.0;}

  // Calculate basis function and derivatives at nodel pts
  phi[0]=(1.0-eta)/2.0;
  phi[1]=(1.0+eta)/2.0;
  dphide[0]=-0.5;
  dphide[1]=0.5;
  
  // Caculate basis function and derivative at GP.
  dz=0.5*(z[1]-z[0]);
  zz=0.0;
  uu=0.0;
  duu=0.0;
  uuold=0.0;
  duuold=0.0;
  for (int i=0; i < N; i++) {
    zz += z[i] * phi[i];
    uu += u[i] * phi[i];
    duu += u[i] * dphide[i];
    if (uold) {
      uuold += uold[i] * phi[i];
      duuold += uold[i] * dphide[i];
    }
  }

  return;
}

#endif
