// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_ELASTICITY_H
#define ROL_PDEOPT_ELASTICITY_H

#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "../TOOLS/PDE_FEM.hpp"
#include "../TOOLS/materials.hpp"

#include "ROL_Types.hpp"

template<class Real>
class Elasticity : public PDE_FEM <Real> {
private:
  using GO = typename Tpetra::Map<>::global_ordinal_type;

protected:

  Real E_;
  Real poissonR_;
  bool planeStrain_;

// dbc cases
  int DBC_Case_;
// parametrized loads
  std::vector<Real> param_;
  std::vector<bool> stochParam_;
// geometry and loads information
  Real xmin_;
  Real ymin_;
  Real xmax_;
  Real ymax_;
// body force
  Real bodyforce_Magnitude_;
  Real bodyforce_Angle_;
// traction
  int  traction_Side_;
  Real traction_Magnitude_;
  Real traction_Angle_;
// point load
  Real pointload_loc_x_;
  Real pointload_loc_y_;
  Real pointload_Magnitude_;
  Real pointload_Angle_;
  std::vector<int > my_pointload_bc_;
//

  std::vector<ROL::Ptr<Material<Real> > > material_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > BMat_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > BMatWeighted_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > CBMat_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > NMat_;
  ROL::Ptr<Intrepid::FieldContainer<Real> > NMatWeighted_;
  int materialTensorDim_;

// for boundary integration of traction force
  ROL::Ptr<Intrepid::FieldContainer<Real> > NMatWeighted_Side;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cub_points_side;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cub_weights_side;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cub_points_side_refcell;
  ROL::Ptr<Intrepid::FieldContainer<Real> > cub_points_side_physical;
  ROL::Ptr<Intrepid::FieldContainer<Real> > jacobian_side_refcell;
  ROL::Ptr<Intrepid::FieldContainer<Real> > jacobian_det_side_refcell;
  ROL::Ptr<Intrepid::FieldContainer<Real> > weighted_measure_side_refcell;
  ROL::Ptr<Intrepid::FieldContainer<Real> > value_of_basis_at_cub_points_side_refcell;
  ROL::Ptr<Intrepid::FieldContainer<Real> > transformed_value_of_basis_at_cub_points_side_refcell;
  ROL::Ptr<Intrepid::FieldContainer<Real> > weighted_transformed_value_of_basis_at_cub_points_side_refcell;
  ROL::Ptr<Intrepid::FieldContainer<Real> > tractions;
  ROL::Ptr<Intrepid::FieldContainer<Real> > tractions_on_dofs;

  bool verbose_;

  Real DegreesToRadians(const Real deg) const {
    return deg * static_cast<Real>(M_PI) / static_cast<Real>(180);
  }

  void PrintLoadingInformation(void) {
    if ( verbose_ ) {
      *this->outStream_ << std::endl << std::string(80,'-') << std::endl;
      *this->outStream_ << std::string(20,' ') << "LOADING INFORMATION" << std::endl;
      *this->outStream_ << std::string(80,'-') << std::endl;
      *this->outStream_ << "  Volumetric Force Magnitude: " << param_[0] << std::endl;
      *this->outStream_ << "  Volumetric Force Angle:     " << param_[1] << std::endl;
      *this->outStream_ << "  Traction Side:              " << traction_Side_ << std::endl;
      *this->outStream_ << "  Traction Force Magnitude:   " << param_[2] << std::endl;
      *this->outStream_ << "  Traction Force Angle:       " << param_[3] << std::endl;
      *this->outStream_ << "  Point Force Location:       " << "(" << pointload_loc_x_
                       << ", " << pointload_loc_y_ << ")" << std::endl;
      *this->outStream_ << "  Point Force Magnitude:      " << param_[4] << std::endl;
      *this->outStream_ << "  Point Force Angle:          " << param_[5] << std::endl;
      *this->outStream_ << std::string(80,'-') << std::endl << std::endl;
    }
  }
//
public:

  Elasticity() {}
  virtual ~Elasticity() {}

  virtual void Initialize(const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                          const Teuchos::RCP<Teuchos::ParameterList> &parlist,
                          const ROL::Ptr<std::ostream> &outStream) {
    /****************************************************************************/
    /*** Initialize the base PDE_FEM class. *************************************/
    /****************************************************************************/
    PDE_FEM<Real>::Initialize(comm,parlist,outStream);

    /****************************************************************************/
    /*** Grab the elasticity information. ***************************************/
    /****************************************************************************/
    Teuchos::ParameterList &Elist = parlist->sublist("Elasticity");
    planeStrain_ = Elist.get<bool>("Plane Strain");
    E_ 	         = Elist.get<Real>("Young's Modulus");
    poissonR_    = Elist.get<Real>("Poisson Ratio");
    Teuchos::ParameterList &Glist = parlist->sublist("Geometry");
    xmin_        = Glist.get<Real>("X0");
    ymin_        = Glist.get<Real>("Y0");
    xmax_        = Glist.get<Real>("Width");
    ymax_        = Glist.get<Real>("Height");
    // DBC cases
    DBC_Case_    = Elist.get<int>("DBC Case");

    /****************************************************************************/
    /*** Initialize mesh / finite element fields / degree-of-freedom manager. ***/
    /****************************************************************************/
    Teuchos::ParameterList &Plist = parlist->sublist("PDE FEM");
    verbose_ = Plist.get("Verbose Output",false);
    int basisOrder = Plist.get<int>("Order of FE Discretization");
    ROL::Ptr<MeshManager<Real> > meshMgr = ROL::makePtr<MeshManager_Rectangle<Real>>(*parlist);
    int spaceDim = 2;
    std::vector<ROL::Ptr<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs(spaceDim,ROL::nullPtr);
    for (int k = 0; k < spaceDim; ++k) {
      if (basisOrder == 1) {
        basisPtrs[k]
          = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >>();
      }
      else if (basisOrder == 2) {
        basisPtrs[k]
          = ROL::makePtr<Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >>();
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
          ">>> (Elasticity::Initialize): Basis Order is out of bounds!");
      }
    }
    PDE_FEM<Real>::SetDiscretization(meshMgr,basisPtrs);
    PDE_FEM<Real>::printMeshData(*outStream);
  }

  // for rectangular domain
  virtual void process_loading_information(const Teuchos::RCP<Teuchos::ParameterList> &parlist) { 
    Teuchos::ParameterList &Elist = parlist->sublist("Elasticity");
    bodyforce_Magnitude_ = Elist.get<Real>("Bodyforce Magnitude");
    bodyforce_Angle_     = Elist.get<Real>("Bodyforce Angle");
    traction_Side_       = Elist.get<int>("Traction Side");
    traction_Magnitude_  = Elist.get<Real>("Traction Magnitude");
    traction_Angle_      = Elist.get<Real>("Traction Angle");
    pointload_loc_x_     = Elist.get<Real>("Pointload Location X");
    pointload_loc_y_     = Elist.get<Real>("Pointload Location Y");
    pointload_Magnitude_ = Elist.get<Real>("Pointload Magnitude");
    pointload_Angle_     = Elist.get<Real>("Pointload Angle");

    param_.clear(); param_.resize(6);
    param_[0] = bodyforce_Magnitude_;
    param_[1] = DegreesToRadians(bodyforce_Angle_);
    param_[2] = traction_Magnitude_;
    param_[3] = DegreesToRadians(traction_Angle_);
    param_[4] = pointload_Magnitude_;
    param_[5] = DegreesToRadians(pointload_Angle_);

    stochParam_.clear(); stochParam_.resize(6);
    stochParam_[0] = Elist.get("Stochastic Bodyforce Magnitude",false);
    stochParam_[1] = Elist.get("Stochastic Bodyforce Angle",false);
    stochParam_[2] = Elist.get("Stochastic Traction Magnitude",false);
    stochParam_[3] = Elist.get("Stochastic Traction Angle",false);
    stochParam_[4] = Elist.get("Stochastic Pointload Magnitude",false);
    stochParam_[5] = Elist.get("Stochastic Pointload Angle",false);

    PrintLoadingInformation();
    
    this->my_nbc_.push_back(traction_Side_);

    // calculate the distance of the point load to the boundary
    // and move it to the closesd boundary point only works for
    // rectangular domain first make sure the point load is in
    // the domain
    TEUCHOS_TEST_FOR_EXCEPTION( pointload_loc_x_ < xmin_, std::invalid_argument,
      ">>> (elasticity::process_loading_information): x location of point load is less than xmin!");
    TEUCHOS_TEST_FOR_EXCEPTION( pointload_loc_x_ > xmax_, std::invalid_argument,
      ">>> (elasticity::process_loading_information): x location of point load is greater than xmax!");
    TEUCHOS_TEST_FOR_EXCEPTION( pointload_loc_y_ < ymin_, std::invalid_argument,
      ">>> (elasticity::process_loading_information): y location of point load is less than ymin!");
    TEUCHOS_TEST_FOR_EXCEPTION( pointload_loc_y_ > ymax_, std::invalid_argument,
      ">>> (elasticity::process_loading_information): y location of point load is greater than ymax!");
    //    
    Real distx1 = std::abs(pointload_loc_x_ - xmin_);
    Real distx2 = std::abs(pointload_loc_x_ - xmax_);
    Real movx = std::min(distx1, distx2);
    Real disty1 = std::abs(pointload_loc_y_ - ymin_);
    Real disty2 = std::abs(pointload_loc_y_ - ymax_);
    Real movy = std::min(disty1, disty2);
    // slight perturbation 
    // pointload will be moved to the closest boundary node eventually
    // perturbation trick to avoid parrallel issues
    Real eps = 1e-8 * (std::min((xmax_-xmin_), (ymax_-ymin_)));
    if(movx <= movy && distx1 <= distx2) {
      // mov load to the left boundary
      pointload_loc_x_ = xmin_ + eps;
      my_pointload_bc_.push_back(3);
      // perturb y
      if(disty1 <= disty2) {
        pointload_loc_y_ = pointload_loc_y_ + eps;
      }
      else {
        pointload_loc_y_ = pointload_loc_y_ - eps;
      }
    }
    else if(movx <= movy && distx1 > distx2) {
      // mov load to the right boundary
      pointload_loc_x_ = xmax_ - eps;
      my_pointload_bc_.push_back(1);
      //perturb y	
      if(disty1 <= disty2) {
        pointload_loc_y_ = pointload_loc_y_ + eps;
      }
      else {
        pointload_loc_y_ = pointload_loc_y_ - eps;
      }
    }
    else if(movx > movy && disty1 <= disty2) {
      // mov load to the bottom boundary
      pointload_loc_y_ = ymin_ + eps;
      my_pointload_bc_.push_back(0);
      // perturb x
      if(distx1 <= distx2) {
        pointload_loc_x_ = pointload_loc_x_ + eps;
      }
      else {
        pointload_loc_x_ = pointload_loc_x_ - eps;
      }
    }
    else {
      // mov load to the top boundary
      pointload_loc_y_ = ymax_ - eps;
      my_pointload_bc_.push_back(2);
      // perturb x
      if(distx1 <= distx2) {
        pointload_loc_x_ = pointload_loc_x_ + eps;
      }
      else {
        pointload_loc_x_ = pointload_loc_x_ - eps;
      }
    }
   
    // print to check
    if ( verbose_ ) {
      *this->outStream_ << "Load processing finished." << std::endl;
      *this->outStream_ << "My nbc numbers: ";
      for(unsigned i=0; i<this->my_nbc_.size(); i++) {
        *this->outStream_ << this->my_nbc_[i];
      }
      *this->outStream_ << std::endl;
      *this->outStream_ << "My pbc numbers: ";
      for(unsigned i=0; i<my_pointload_bc_.size(); i++) {
        *this->outStream_ << my_pointload_bc_[i];
      }
      *this->outStream_ << std::endl;
      *this->outStream_ << "My pointload location: " << pointload_loc_x_
                        << ", " << pointload_loc_y_ << std::endl;
    }
  }

  virtual void CreateMaterial(void) {
    int numCells = PDE_FEM<Real>::GetNumCells();
    int spaceDim = PDE_FEM<Real>::GetSpaceDim();
    for(int i = 0; i < numCells; ++i) {
      ROL::Ptr<Material<Real> > CellMaterial = ROL::makePtr<Material<Real>>();
      CellMaterial-> InitializeMaterial(spaceDim, planeStrain_, E_, poissonR_);
      materialTensorDim_ = CellMaterial->GetMaterialTensorDim();
      material_.push_back(CellMaterial);
    }
  }

  virtual void ComputeLocalSystemMats(void) {
    int spaceDim = PDE_FEM<Real>::GetSpaceDim();
    int numCells = PDE_FEM<Real>::GetNumCells();
    int numCubPoints = PDE_FEM<Real>::GetNumCubPoints();
    int lfs = PDE_FEM<Real>::GetLocalFieldSize();
    int full_lfs = lfs * spaceDim;
    
    this->gradgradMats_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells, full_lfs, full_lfs);
    this->valvalMats_   = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells, full_lfs, full_lfs);
    
    // construct mats
    CreateMaterial();
    
    BMat_
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells, full_lfs, numCubPoints, materialTensorDim_);
    BMatWeighted_
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells, full_lfs, numCubPoints, materialTensorDim_);
    CBMat_
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells, full_lfs, numCubPoints, materialTensorDim_);
    NMat_
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells, full_lfs, numCubPoints, spaceDim);
    NMatWeighted_
      = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells, full_lfs, numCubPoints, spaceDim);
    Construct_N_B_mats();
    Construct_CBmats();

    // compute local grad.grad (stiffness) matrices
    Intrepid::FunctionSpaceTools::integrate<Real>(*this->gradgradMats_,
                                                  *CBMat_,
                                                  *BMatWeighted_,
                                                  Intrepid::COMP_CPP);
    // compute local val.val (mass) matrices
    Intrepid::FunctionSpaceTools::integrate<Real>(*this->valvalMats_,
                                                  *NMat_,
                                                  *NMatWeighted_,
                                                  Intrepid::COMP_CPP);
  } 
  
  
  // new
  // constructing Nmat on side
  void Construct_Nmat_on_Side(int numCub) {
    int spaceDim = PDE_FEM<Real>::GetSpaceDim();
    int lfs = PDE_FEM<Real>::GetLocalFieldSize();
    int full_lfs = lfs * spaceDim;
    NMatWeighted_Side = ROL::makePtr<Intrepid::FieldContainer<Real>>(1, full_lfs, numCub, spaceDim);
    NMatWeighted_Side -> initialize(0.0);
    
    if(spaceDim == 2) {
      for (int j=0; j < numCub; ++j) {
        for (int k=0; k<lfs; ++k) {
          (*NMatWeighted_Side)(0, k*spaceDim+0, j, 0)
            = (*weighted_transformed_value_of_basis_at_cub_points_side_refcell)(0, k, j);
          (*NMatWeighted_Side)(0, k*spaceDim+1, j, 1)
            = (*weighted_transformed_value_of_basis_at_cub_points_side_refcell)(0, k, j);
        }
      }
    }

    if(spaceDim == 3) {
      for (int j=0; j < numCub; ++j) {
        for (int k=0; k<lfs; ++k) {
          (*NMatWeighted_Side)(0, k*spaceDim+0, j, 0)
            = (*weighted_transformed_value_of_basis_at_cub_points_side_refcell)(0, k, j);
          (*NMatWeighted_Side)(0, k*spaceDim+1, j, 1)
            = (*weighted_transformed_value_of_basis_at_cub_points_side_refcell)(0, k, j);
          (*NMatWeighted_Side)(0, k*spaceDim+2, j, 2)
            = (*weighted_transformed_value_of_basis_at_cub_points_side_refcell)(0, k, j);
        }
      }
    } 
  }


  // adding point load to right hand side vector
  virtual void AddPointLoadToRHS(void) {
    int n_pbc = my_pointload_bc_.size(); 
    for (int i=0; i<n_pbc; ++i) {
      int pbc_num = my_pointload_bc_[i];
      for (int j=0; j<this->myBoundaryCellIds_[pbc_num].size(); ++j) {
	int myGlobalCellId = this->myBoundaryCellIds_[pbc_num][j];
	// apply possible point loads
	this->check_and_Apply_PointLoad_By_Coords(myGlobalCellId, pbc_num);
      }
    }
  }

  // adding traction boundary data into right hand side vector
  virtual void ModifyLocalForceVecWithSideTractions(void) {

    int cellDim = PDE_FEM<Real>::GetSpaceDim();
    int lfs = PDE_FEM<Real>::GetLocalFieldSize(); //number of dof each dimension
    int full_lfs = lfs*cellDim;
    int numNodesPerCell = this->numNodesPerCell_;

    shards::CellTopology sideType(shards::getCellTopologyData< shards::Line<> >());
    int cubDegree = 10;                                                             
    Intrepid::DefaultCubatureFactory<Real> cubFactory; // create cubature factory
    ROL::Ptr<Intrepid::Cubature<Real> > sideCub = cubFactory.create(sideType, cubDegree);
    int numCubPointsSide = sideCub->getNumPoints();
    if ( verbose_ ) {
      *this->outStream_<<"numCubPointsSide: "<<numCubPointsSide<<std::endl;
    }

    int sideDim = this->sideDim_;
    int n_nbc = this->my_nbc_.size(); 
    ROL::Ptr<Intrepid::FieldContainer<Real> > thiscellNodes;
    thiscellNodes = ROL::makePtr<Intrepid::FieldContainer<Real>>(1, numNodesPerCell, cellDim);
    Intrepid::FieldContainer<Real> &nodes = *(this->meshMgr_)->getNodes();

    if ( verbose_ ) {
      *this->outStream_<<"n_nbc: "<<n_nbc<<std::endl;
    }
    for (int i=0; i<n_nbc; i++) {
      int nbc_num = this->my_nbc_[i];
      //std::cout << "nbc_num: " << nbc_num << std::endl;
      for (int j=0; j<this->myBoundaryCellIds_[nbc_num].size(); ++j) {
        int myGlobalCellId = this->myBoundaryCellIds_[nbc_num][j];
        //*this->outStream_ << "myGlobalCellId: " << myGlobalCellId << std::endl;

        // apply traction        
        for (int m=0; m<numNodesPerCell; ++m) {
          for (int n=0; n<cellDim; ++n) {
            (*thiscellNodes)(0, m, n) = nodes(this->ctn_(myGlobalCellId, m), n);
          }
        }
        //std::cout << "first node coords: " << (*thiscellNodes)(0, 0, 0)
        //          << ", " << (*thiscellNodes)(0, 0, 1) << std::endl;
        cub_points_side
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCubPointsSide, sideDim);
        cub_weights_side
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCubPointsSide);
        cub_points_side_refcell
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCubPointsSide, cellDim);
        cub_points_side_physical
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(1, numCubPointsSide, cellDim);
        jacobian_side_refcell
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(1, numCubPointsSide, cellDim, cellDim);
        jacobian_det_side_refcell
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(1, numCubPointsSide);
        weighted_measure_side_refcell
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(1, numCubPointsSide);
        value_of_basis_at_cub_points_side_refcell
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(lfs, numCubPointsSide);
        transformed_value_of_basis_at_cub_points_side_refcell
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(1, lfs, numCubPointsSide);
        weighted_transformed_value_of_basis_at_cub_points_side_refcell
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(1, lfs, numCubPointsSide);
        tractions
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(1, numCubPointsSide, cellDim);
        tractions_on_dofs
          = ROL::makePtr<Intrepid::FieldContainer<Real>>(1, full_lfs);
           
        // compute traction b.c. contributions and adjust rhs
        sideCub->getCubature(*cub_points_side, *cub_weights_side);
         
        // compute geometric cell information
        Intrepid::CellTools<Real>::mapToReferenceSubcell(*cub_points_side_refcell,
                                                         *cub_points_side,
                                                          sideDim,
                                                          nbc_num,
                                                          this->cellType_);

        Intrepid::CellTools<Real>::setJacobian(*jacobian_side_refcell,
                                               *cub_points_side_refcell,
                                               *thiscellNodes,
                                                this->cellType_);

        Intrepid::CellTools<Real>::setJacobianDet(*jacobian_det_side_refcell,
                                                  *jacobian_side_refcell);
         
        // compute weighted edge measure
        Intrepid::FunctionSpaceTools::computeEdgeMeasure<Real>(*weighted_measure_side_refcell,
                                                               *jacobian_side_refcell,
                                                               *cub_weights_side,
                                                                nbc_num,
                                                                this->cellType_);
         
        // tabulate values of basis functions at side cubature points, in the reference parent cell domain
        (*this->basisPtrs_[0]).getValues(*value_of_basis_at_cub_points_side_refcell,
                                         *cub_points_side_refcell,
                                          Intrepid::OPERATOR_VALUE);
        // transform
        Intrepid::FunctionSpaceTools::HGRADtransformVALUE<Real>(
          *transformed_value_of_basis_at_cub_points_side_refcell,
          *value_of_basis_at_cub_points_side_refcell);
         
        // multiply with weighted measure
        Intrepid::FunctionSpaceTools::multiplyMeasure<Real>(
          *weighted_transformed_value_of_basis_at_cub_points_side_refcell,
          *weighted_measure_side_refcell,
          *transformed_value_of_basis_at_cub_points_side_refcell);

        Construct_Nmat_on_Side(numCubPointsSide);
         
        // compute Neumann data
        // map side cubature points in reference parent cell domain to physical space
        Intrepid::CellTools<Real>::mapToPhysicalFrame(*cub_points_side_physical,
                                                      *cub_points_side_refcell,
                                                      *thiscellNodes,
                                                       this->cellType_);
        //std::cout << "cub_points_side_physical:" << (*cub_points_side_physical)(0,0,0)
        //          << ", " << (*cub_points_side_physical)(0,0,1) << std::endl;
        // now compute data
        std::vector<Real> x(cellDim), F(cellDim);
        for (int m = 0; m < numCubPointsSide; ++m) {
          for (int n = 0; n < cellDim; ++n) {
            x[n] = (*cub_points_side_physical)(0,m,n);
          }
          funcRHS_NBC(F, x);
          for (int n = 0; n < cellDim; ++n) {
            (*tractions)(0,m,n) = F[n];
          }
        }
        Intrepid::FunctionSpaceTools::integrate<Real>(*tractions_on_dofs,
                                                      *tractions,
                                                      *NMatWeighted_Side,
                                                      Intrepid::COMP_CPP);
         
        // adjust RHS
        for (int m=0; m < full_lfs; ++m) {
          (*this->datavalVecF_)(this->find_local_index(myGlobalCellId), m)
            += (*tractions_on_dofs)(0, m);
        }
 
        //check tractions on dofs
        //*this->outStream_<<"tractions_on_dofs: ";
        //for(int m=0; m<full_lfs; ++m)
	//  *this->outStream_<<(*tractions_on_dofs)(0, m)<<", ";
        //*this->outStream_<<std::endl;
      }
    }
  }
   

  virtual void ComputeLocalForceVec(void) {
    int spaceDim = PDE_FEM<Real>::GetSpaceDim();
    int numCells = PDE_FEM<Real>::GetNumCells();
    int numCubPoints = PDE_FEM<Real>::GetNumCubPoints();
    int lfs = PDE_FEM<Real>::GetLocalFieldSize();
    int full_lfs = lfs * spaceDim;
    this->dataF_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells, numCubPoints, spaceDim);
    this->datavalVecF_ = ROL::makePtr<Intrepid::FieldContainer<Real>>(numCells, full_lfs);

    std::vector<Real> x(spaceDim), F(spaceDim);
    for (int i = 0; i < numCells; ++i) { // evaluate functions at these points
      for (int j = 0; j < numCubPoints; ++j) {
        for (int k = 0; k < spaceDim; ++k) {
          x[k] = (*this->cubPointsPhysical_)(i,j,k);
        }
        funcRHS_BodyForce(F, x);
        for (int k = 0; k < spaceDim; ++k) {
          (*this->dataF_)(i,j,k) = F[k];
        }
      }
    }
    Intrepid::FunctionSpaceTools::integrate<Real>(*this->datavalVecF_, // compute local data.val vectors for RHS F
                                                  *this->dataF_,
                                                  *NMatWeighted_,
                                                   Intrepid::COMP_CPP);

    // new
    if ( verbose_ ) {
      *this->outStream_ << "Modifying local force vectors using boundary tractions" << std::endl;
    }
    ModifyLocalForceVecWithSideTractions();
    AddPointLoadToRHS();
    if ( verbose_ ) {
      *this->outStream_ << "Modifying Done!" << std::endl;
    }
  }

  virtual void updateF(const std::vector<Real> &param) {
    int ind = 0;
    for (int i = 0; i < 6; ++i) {
      if ( stochParam_[i] ) {
        param_[i] = param[ind];
        ind++;
      }
    }

    ComputeLocalForceVec();
    PDE_FEM<Real>::AssembleRHSVector();
    PDE_FEM<Real>::VectorRemoveDBC();

    PrintLoadingInformation();
  }


  void Construct_N_B_mats(void) {
    //std::cout<<"Computing N and B mats."<<std::endl;
    BMat_->initialize(0.0);
    NMat_->initialize(0.0);
    BMatWeighted_->initialize(0.0);
    NMatWeighted_->initialize(0.0);
    
    if(this->spaceDim_==2) {
      for (int i=0; i<this->numCells_; ++i) { // evaluate functions at these points
        for (int j=0; j<this->numCubPoints_; ++j) {
      	  for (int k=0; k<this->lfs_; ++k) {
            (*NMat_)(i, k*this->spaceDim_+0, j, 0) = (*this->valPhysical_)(i, k, j);
            (*NMat_)(i, k*this->spaceDim_+1, j, 1) = (*this->valPhysical_)(i, k, j);
            (*NMatWeighted_)(i, k*this->spaceDim_+0, j, 0) = (*this->valPhysicalWeighted_)(i, k, j);
            (*NMatWeighted_)(i, k*this->spaceDim_+1, j, 1) = (*this->valPhysicalWeighted_)(i, k, j);
                
            (*BMat_)(i, k*this->spaceDim_+0, j, 0) = (*this->gradPhysical_)(i, k, j, 0);
            (*BMat_)(i, k*this->spaceDim_+1, j, 1) = (*this->gradPhysical_)(i, k, j, 1);
            (*BMat_)(i, k*this->spaceDim_+0, j, 2) = (*this->gradPhysical_)(i, k, j, 1);
            (*BMat_)(i, k*this->spaceDim_+1, j, 2) = (*this->gradPhysical_)(i, k, j, 0);
            
            (*BMatWeighted_)(i, k*this->spaceDim_+0, j, 0) = (*this->gradPhysicalWeighted_)(i, k, j, 0);
            (*BMatWeighted_)(i, k*this->spaceDim_+1, j, 1) = (*this->gradPhysicalWeighted_)(i, k, j, 1);
            (*BMatWeighted_)(i, k*this->spaceDim_+0, j, 2) = (*this->gradPhysicalWeighted_)(i, k, j, 1);
            (*BMatWeighted_)(i, k*this->spaceDim_+1, j, 2) = (*this->gradPhysicalWeighted_)(i, k, j, 0);
	  }
        }
      }
    }

    if(this->spaceDim_==3) {
      for (int i=0; i<this->numCells_; ++i) { // evaluate functions at these points
        for (int j=0; j<this->numCubPoints_; ++j) {
          for (int k=0; k<this->lfs_; ++k) {
            (*NMat_)(i, k*this->spaceDim_+0, j, 0) = (*this->valPhysical_)(i, k, j);
            (*NMat_)(i, k*this->spaceDim_+1, j, 1) = (*this->valPhysical_)(i, k, j);
            (*NMat_)(i, k*this->spaceDim_+2, j, 2) = (*this->valPhysical_)(i, k, j);
            (*NMatWeighted_)(i, k*this->spaceDim_+0, j, 0) = (*this->valPhysicalWeighted_)(i, k, j);
            (*NMatWeighted_)(i, k*this->spaceDim_+1, j, 1) = (*this->valPhysicalWeighted_)(i, k, j);
            (*NMatWeighted_)(i, k*this->spaceDim_+2, j, 2) = (*this->valPhysicalWeighted_)(i, k, j);
            
            (*BMat_)(i, k*this->spaceDim_+0, j, 0) = (*this->gradPhysical_)(i, k, j, 0);
            (*BMat_)(i, k*this->spaceDim_+1, j, 1) = (*this->gradPhysical_)(i, k, j, 1);
            (*BMat_)(i, k*this->spaceDim_+2, j, 2) = (*this->gradPhysical_)(i, k, j, 2);
            (*BMat_)(i, k*this->spaceDim_+1, j, 3) = (*this->gradPhysical_)(i, k, j, 2);
            (*BMat_)(i, k*this->spaceDim_+2, j, 3) = (*this->gradPhysical_)(i, k, j, 1);
            (*BMat_)(i, k*this->spaceDim_+0, j, 4) = (*this->gradPhysical_)(i, k, j, 2);
            (*BMat_)(i, k*this->spaceDim_+2, j, 4) = (*this->gradPhysical_)(i, k, j, 0);
            (*BMat_)(i, k*this->spaceDim_+0, j, 5) = (*this->gradPhysical_)(i, k, j, 1);
            (*BMat_)(i, k*this->spaceDim_+1, j, 5) = (*this->gradPhysical_)(i, k, j, 0);
                
            (*BMatWeighted_)(i, k*this->spaceDim_+0, j, 0) = (*this->gradPhysicalWeighted_)(i, k, j, 0);
            (*BMatWeighted_)(i, k*this->spaceDim_+1, j, 1) = (*this->gradPhysicalWeighted_)(i, k, j, 1);
            (*BMatWeighted_)(i, k*this->spaceDim_+2, j, 2) = (*this->gradPhysicalWeighted_)(i, k, j, 2);
            (*BMatWeighted_)(i, k*this->spaceDim_+1, j, 3) = (*this->gradPhysicalWeighted_)(i, k, j, 2);
            (*BMatWeighted_)(i, k*this->spaceDim_+2, j, 3) = (*this->gradPhysicalWeighted_)(i, k, j, 1);
            (*BMatWeighted_)(i, k*this->spaceDim_+0, j, 4) = (*this->gradPhysicalWeighted_)(i, k, j, 2);
            (*BMatWeighted_)(i, k*this->spaceDim_+2, j, 4) = (*this->gradPhysicalWeighted_)(i, k, j, 0);
            (*BMatWeighted_)(i, k*this->spaceDim_+0, j, 5) = (*this->gradPhysicalWeighted_)(i, k, j, 1);
            (*BMatWeighted_)(i, k*this->spaceDim_+1, j, 5) = (*this->gradPhysicalWeighted_)(i, k, j, 0);
          }
        }
      }
    } 
  }

  // test matrices
  virtual void test_mats(void) {
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_Jaco_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_Grad_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_N_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_B_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_K_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_M_Mat;
    ROL::Ptr<Intrepid::FieldContainer<Real> > test_F_Vec;
    
    test_Jaco_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->spaceDim_, this->spaceDim_);
    test_Grad_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->lfs_, this->spaceDim_);
    test_N_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numLocalDofs_, this->spaceDim_);
    test_B_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numLocalDofs_, materialTensorDim_);
    test_K_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numLocalDofs_, this->numLocalDofs_);
    test_M_Mat = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numLocalDofs_, this->numLocalDofs_);
    test_F_Vec = ROL::makePtr<Intrepid::FieldContainer<Real>>(this->numLocalDofs_, 1);
    
    for(int i=0; i<this->spaceDim_; i++) {
       for(int j=0; j<this->spaceDim_; j++) {
         (*test_Jaco_Mat)(i, j) = (*this->cellJac_)(0, 0, i, j);
       }
    }
    for(int i=0; i<this->numLocalDofs_; i++) {
      for(int j=0; j<this->spaceDim_; j++) {
        if(i<this->lfs_) {
          (*test_Grad_Mat)(i, j) = (*this->gradReference_)(i, 0, j);
        }
        (*test_N_Mat)(i, j) = (*NMat_)(0, i, 0, j);	
      }
      for(int j=0; j<materialTensorDim_; j++) {
        (*test_B_Mat)(i, j) = (*BMat_)(0, i, 0, j);	
      }
      for(int j=0; j<this->numLocalDofs_; j++) {
        (*test_K_Mat)(i, j) = (*this->gradgradMats_)(0, i, j);
        (*test_M_Mat)(i, j) =   (*this->valvalMats_)(0, i, j);
      }
      (*test_F_Vec)(i, 0) = (*this->datavalVecF_)(0, i);
    }
    std::cout<<*test_Jaco_Mat<<std::endl;
    std::cout<<*test_Grad_Mat<<std::endl;
    std::cout<<*test_N_Mat<<std::endl;
    std::cout<<*test_B_Mat<<std::endl;
    std::cout<<*test_K_Mat<<std::endl;
    std::cout<<*test_M_Mat<<std::endl;
    std::cout<<*test_F_Vec<<std::endl;
  }


  virtual void Construct_CBmats(void) {
    //std::cout<<"Computing CB mats."<<std::endl;
    CBMat_->initialize(0.0);
    for (int i=0; i<this->numCells_; ++i) {
      ROL::Ptr<Intrepid::FieldContainer<Real> > materialMat = material_[i]->GetMaterialTensor();
      for (int j=0; j<this->numCubPoints_; ++j) {
        for (int m=0; m<this->lfs_*this->spaceDim_; m++) {
          for (int n=0; n<materialTensorDim_; n++) {
            for (int k=0; k<materialTensorDim_; k++) {
              (*CBMat_)(i, m, j, n) += (*BMat_)(i, m, j, k) * (*materialMat)(k, n);
            }
          }
        }
      }
    }
  }

// load functions, can be parametrized
// modification for stochastic loads should be made here
  virtual void funcRHS_BodyForce(std::vector<Real> &F,
                           const std::vector<Real> &x) const {
    F.clear(); F.resize(this->spaceDim_);
    F[0] = param_[0]*std::cos(param_[1]);
    F[1] = param_[0]*std::sin(param_[1]);
  }

  virtual void funcRHS_NBC(std::vector<Real> &F,
                     const std::vector<Real> &x) const {
    F.clear(); F.resize(this->spaceDim_);
    F[0] = param_[2]*std::cos(param_[3]);
    F[1] = param_[2]*std::sin(param_[3]);
  }

  virtual void funcRHS_PtLoad(std::vector<Real> &F) const {
    F.clear(); F.resize(this->spaceDim_);
    F[0] = param_[4]*std::cos(param_[5]);
    F[1] = param_[4]*std::sin(param_[5]);
  }
//
//


  virtual void ApplyPointLoad(const int pointload_bc,
                              const int globalCellNum,
                              const std::vector<int> &localNodeNum,
                              const std::vector<Real> &coord1,
                              const std::vector<Real> &coord2) { 
    ROL::Ptr<Intrepid::FieldContainer<GO> > nodeDofs = this->dofMgr_->getNodeDofs();
    bool isLoadPosContainedInCurrentSegment = false;
    int whichNodeIsCloser = -1;
    // if update F, provides parametrized computation of F[0] and F[1]
    std::vector<Real> F;
    funcRHS_PtLoad(F);
    //
    Real x11 = coord1[0], x12 = coord1[1];
    Real x21 = coord2[0], x22 = coord2[1];
    Real fx = pointload_loc_x_;
    Real fy = pointload_loc_y_;
    int pbc_num = pointload_bc;

    if(pbc_num == 0 || pbc_num == 2) {
      if ( ((x11-fx)*(x21-fx)<0) && ((x12-fy)*(x22-fy)>0) ) {
        isLoadPosContainedInCurrentSegment = true;
        if(std::abs(x11-fx) <= std::abs(x21-fx)) {
          whichNodeIsCloser = 1;
        }
        else {
          whichNodeIsCloser = 2;
        }
      }
    }
    if(pbc_num == 1 || pbc_num == 3) {
      if ( ((x11-fx)*(x21-fx)>0) && ((x12-fy)*(x22-fy)<0) ) {
        isLoadPosContainedInCurrentSegment = true;
        if(std::abs(x12-fy) <= std::abs(x22-fy)) {
          whichNodeIsCloser = 1;
        }
        else {
          whichNodeIsCloser = 2;
        }
      }
    }
    if(isLoadPosContainedInCurrentSegment) {
      int ctn = this->ctn_(globalCellNum, localNodeNum[whichNodeIsCloser-1]);
      this->myPointLoadDofs_.resize(2);
      this->myPointLoadVals_.resize(2);
      this->myPointLoadDofs_[0] = (*nodeDofs)(ctn, 0);
      this->myPointLoadDofs_[1] = (*nodeDofs)(ctn, 1);
      this->myPointLoadVals_[0] = F[0];
      this->myPointLoadVals_[1] = F[1];
      // print to check
      if ( verbose_ ) {
        *this->outStream_ << "Point load applied, at cell: " << globalCellNum << std::endl;
        *this->outStream_ << "Point load values: " << F[0] << ", " << F[1] << std::endl;
        if(whichNodeIsCloser == 1) {
          *this->outStream_ << "Point load position: " << x11 << ", " << x12 << std::endl;
        }
        else {
          *this->outStream_ << "Point load position: " << x21 << ", " << x22 << std::endl;
        }
      }
    }
  }
// ApplyPointLoad


// DBC cases, add more to extend
  virtual int check_DBC_Coords_Range( const std::vector<Real> &x ) const {
    // return value :
    // -1: not a DBC node
    //  0: x direction fixed
    //  1: y direction fixed
    //  5: both x, y direction are fixed
    //
    Real x1 = x[0], x2 = x[1], eps(1e-6);
    switch(DBC_Case_) {
      case 0: { // Fix bottom two corners, left corner x, y, right corner only y
        if ( (x2 < ymin_ + eps) && (x1 < xmin_ + eps) ) return 5;
        if ( (x2 < ymin_ + eps) && (x1 > xmax_ - eps) ) return 1;
        break;
      }
      case 1: { // Fix bottom two corners, both x, y
        if ( (x2 < ymin_ + eps) &&
             ( (x1 < xmin_ + eps) || (x1 > xmax_ - eps) ) ) return 5;
        break;
      }
      case 2: { // Fix left boundary, both x, y
        if ( x1 < xmin_ + eps ) return 5;
        break;
      }
      case 3: { // Fix bottom boundary, both x, y
        if ( x2 < ymin_ + eps ) return 5;
        break;
      }
      case 4: { // Fix ALL boundary, both x, y
        return 5;
        break;
      }
    }
    return -1;
  }


}; // class Elasticity

#endif
