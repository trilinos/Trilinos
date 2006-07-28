#include "GLpApp_AdvDiffReactOptModelCreator.hpp"

namespace GLpApp {

AdvDiffReactOptModelCreator::AdvDiffReactOptModelCreator()
  :len_x_(1.0)
  ,len_y_(1.0)
  ,local_nx_(3)
  ,local_ny_(4)
  ,geomFileBase_("")
  ,np_(-1)
  ,normalizeBasis_(false)
  ,beta_(1.0)
  ,reactionRate_(1.0)
  ,x0_(0.0)
  ,p0_(1.0)
{}

void AdvDiffReactOptModelCreator::setupCLP(
  Teuchos::CommandLineProcessor *clp
  )
{
  clp->setOption( "len-x", &len_x_, "Mesh dimension in the x direction (Overridden by --geom-file-base)." );
  clp->setOption( "len-y", &len_y_, "Mesh dimension in the y direction (Overridden by --geom-file-base)." );
  clp->setOption( "local-nx", &local_nx_, "Number of local discretization segments in the x direction (Overridden by --geom-file-base)." );
  clp->setOption( "local-ny", &local_ny_, "Number of local discretization segments in the y direction (Overridden by --geom-file-base)." );
  clp->setOption( "geom-file-base", &geomFileBase_, "Base name of geometry file to read the mesh from." );
  clp->setOption( "np", &np_, "The number of optimization parameters p (If < 0 then all of boundary is used)" );
  clp->setOption( "normalize-basis", "no-normalize-basis", &normalizeBasis_, "Normalize the basis for the parameters p or not." );
  clp->setOption( "beta", &beta_, "Regularization." );
  clp->setOption( "reaction-rate", &reactionRate_, "The rate of the reaction" );
  clp->setOption( "x0", &x0_, "Initial guess for the state." );
  clp->setOption( "p0", &p0_, "Initial guess or nonminal value for optimization parameters." );
}

Teuchos::RefCountPtr<AdvDiffReactOptModel>
AdvDiffReactOptModelCreator::createModel(
  const Teuchos::RefCountPtr<const Epetra_Comm>     &comm
  ,std::ostream                                     *out
  ) const
{
  return Teuchos::rcp(
    new GLpApp::AdvDiffReactOptModel(
      comm,beta_,len_x_,len_y_,local_nx_,local_ny_,geomFileBase_.c_str()
      ,np_,x0_,p0_,reactionRate_,normalizeBasis_
      )
    );
}

} // namespace GLpApp
