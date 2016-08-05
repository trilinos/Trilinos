// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  pde.hpp
    \brief Implements the local PDE interface for the Stefan_Boltzmann control problem.
*/

#ifndef PDEOPT_PDE_STEFAN_BOLTZMANN_HPP
#define PDEOPT_PDE_STEFAN_BOLTZMANN_HPP

#include "../TOOLS/pde.hpp"
#include "../TOOLS/fe.hpp"

#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "Intrepid_CellTools.hpp"

#include "Teuchos_RCP.hpp"

template <class Real>
class PDE_Stefan_Boltzmann : public PDE<Real> {
private:

  Real alpha_;
  
  int lfs_;
  int spaceDim_;
  int numCubPerCell_;
  int numCubPerSide_;
  int numNodesPerCell_;
  int numEdgesPerCell_;
   
  int numCells_vol_;
  int n_bc_segments_;
  int n_bc_segments_w_fe_;
  std::vector<int > n_bc_sub_segments_;
  std::vector<std::vector<int> > sideIds_;
  
  std::vector<int > boundary_type_;
  std::vector<std::vector<int > > numCells_bc_;
  std::vector<std::vector<std::vector<int> > > bc_cell_local_id_;
  
  std::vector<std::vector<int > > local_dofs_on_sides_;
  
  //Teuchos::RCP<shards::CellTopology> vol_cellTopo_;   // base (parent) cell topology
  //Teuchos::RCP<shards::CellTopology>  bc_sideTopo_;   // side (subcell) topology
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > basisPtrs_;
  
  Teuchos::RCP<FE<Real> > fe_vol_; //volume
  std::vector< int > fe_bc_idx_; // bc segments that need fe objects (not including dirichlet bc)
  std::vector<std::vector<Teuchos::RCP<FE<Real> > > > fe_bc_;
 
  // get the physical coords of dof points
  Teuchos::RCP<Intrepid::FieldContainer<Real> > dof_points_physical_;
  Teuchos::RCP<Intrepid::FieldContainer<Real> > vol_cub_points_physical_;
  std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > side_cub_points_physical_;
   
  //side_cub_points_physical
  //        = Teuchos::rcp(new Intrepid::FieldContainer<Real>(nCellBC, numCubPointsSide, cellDim));


public:

  PDE_Stefan_Boltzmann(Teuchos::ParameterList &parlist)
  {
    int basisOrder = parlist.sublist("Problem").get("Order of FE discretization", 1);
    alpha_ = parlist.sublist("Problem").get("Robin boundary alpha", 0.1);
    
    Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > basisPtr;
    if (basisOrder == 1) {
      basisPtr = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C1_FEM<Real, Intrepid::FieldContainer<Real> >);
    }
    else if (basisOrder == 2) {
      basisPtr = Teuchos::rcp(new Intrepid::Basis_HGRAD_QUAD_C2_FEM<Real, Intrepid::FieldContainer<Real> >);
    }
    basisPtrs_.resize(1, Teuchos::null);
    basisPtrs_[0] = basisPtr;
    lfs_ = basisPtr -> getCardinality();
  }
  
  
  std::vector<Teuchos::RCP<Intrepid::Basis<Real, Intrepid::FieldContainer<Real> > > > getFields()
  {
    return basisPtrs_;
  }
  
	
  void setCellNodes(const Teuchos::RCP<Intrepid::FieldContainer<Real> > &cellNodes,
                    const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &bdryCellNodes, 
                    const std::vector<std::vector<std::vector<int> > > &bdryCellLocIds ) {
    bc_cell_local_id_ = bdryCellLocIds;
    std::cout << "setCellNodes before assign_boundary!" << std::endl;
    assign_boundary_type(bdryCellNodes);
    std::cout << "setCellNodes before construct FE!" << std::endl;
    construct_FE_objects(cellNodes, bdryCellNodes);
  }


// user assign boundary types to each boundary segments
   void assign_boundary_type( const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &bdryCells)
  {
    n_bc_segments_ = bdryCells.size();
    boundary_type_.resize(n_bc_segments_);
    sideIds_.resize(n_bc_segments_);
    sideIds_[0].resize(2);
    sideIds_[0][0] = 0;
    sideIds_[0][1] = 3;
    sideIds_[1].resize(1);
    sideIds_[1][0] = 1;
    sideIds_[2].resize(1);
    sideIds_[2][0] = 2;
    //    N
    // D     R
    //    D
    boundary_type_[0] = 0; // Dirichlet
    boundary_type_[1] = 2; // Robin 
    boundary_type_[2] = 1; // Neumann
    fe_bc_idx_.resize(n_bc_segments_);
    int k=0;
    for(int i=0; i<n_bc_segments_;i++)
      {
	fe_bc_idx_[i] = -1; // 
	if(boundary_type_[i] > 0)
	  {
	    fe_bc_idx_[i] = k;
	    k++;
	  }
      }
    n_bc_segments_w_fe_ = k;
  }
  


   void construct_FE_objects ( const Teuchos::RCP<Intrepid::FieldContainer<Real> > & volume_cellNodes, const std::vector<std::vector<Teuchos::RCP<Intrepid::FieldContainer<Real> > > > &bdryCells)
  {
    //vol_cellTopo_ = Teuchos::rcp(new shards::CellTopology(basisPtrs_[0]->getBaseCellTopology()));
    shards::CellTopology vol_cellTopo_ = basisPtrs_[0]->getBaseCellTopology();
    spaceDim_ = vol_cellTopo_.getDimension();
    numNodesPerCell_ = vol_cellTopo_.getNodeCount(); 
    numEdgesPerCell_ = vol_cellTopo_.getEdgeCount(); 
    
    // volume cubature data.
    Intrepid::DefaultCubatureFactory<Real> cubFactory;
    int vol_cubDegree = 4;                               
    Teuchos::RCP<Intrepid::Cubature<Real> > vol_Cub = cubFactory.create(vol_cellTopo_, vol_cubDegree);
    numCubPerCell_ = vol_Cub->getNumPoints();     
    numCells_vol_ = volume_cellNodes->dimension(0);
    fe_vol_ = Teuchos::rcp(new FE<Real > (volume_cellNodes, basisPtrs_[0], vol_Cub));
    //
    get_local_dofs_for_bc_sides();
    dof_points_physical_ = Teuchos::rcp( new Intrepid::FieldContainer<Real >(numCells_vol_, lfs_, spaceDim_));
    fe_vol_->computeDofCoords(dof_points_physical_,volume_cellNodes);
    std::cout << "PDE:construct_FE_objects: fe_vol_ success" << std::endl;
    vol_cub_points_physical_ = fe_vol_->cubPts();
    
    // boundary cubature data
    //bc_sideTopo_ = Teuchos::rcp(new shards::CellTopology(vol_cellTopo_.getCellTopologyData(spaceDim_-1, 0)));
    shards::CellTopology bc_sideTopo_ = vol_cellTopo_.getCellTopologyData(spaceDim_-1, 0);
    int bc_cubDegree = 8; 
    const Teuchos::RCP<Intrepid::Cubature<Real> > bc_Cub = cubFactory.create(bc_sideTopo_, bc_cubDegree);
    numCubPerSide_ = bc_Cub->getNumPoints();
    
    numCells_bc_.resize(n_bc_segments_);
    fe_bc_.resize(n_bc_segments_w_fe_);

    n_bc_sub_segments_.resize(n_bc_segments_);
    
    side_cub_points_physical_.resize(n_bc_segments_w_fe_);
    std::cout << "PDE:construct_FE_objects: n_bc_segments_"<<n_bc_segments_<<" n_bc_seg_w_fe_ "<<n_bc_segments_w_fe_<<" " << std::endl;
    int k = 0;
    for(int i=0; i<n_bc_segments_; i++)
      {
	int bcType = boundary_type_[i];
	n_bc_sub_segments_[i] = bdryCells[i].size();
	std::cout << "n_bc_sub_segments_[" << i << "] = " << n_bc_sub_segments_[i] << std::endl;
	numCells_bc_[i].resize(n_bc_sub_segments_[i]);
	if(bcType > 0)
	  {
	    fe_bc_[k].resize(n_bc_sub_segments_[i]);
	    side_cub_points_physical_[k].resize(n_bc_sub_segments_[i]);
	  }
	
	//actually, fe objects need not be built for dbc
	for (int j=0; j<n_bc_sub_segments_[i]; j++)
	  {
	    const int sideId = sideIds_[i][j];
	    numCells_bc_[i][j] = bdryCells[i][j]->dimension(0);
	    if (bcType > 0)
	      {
		std::cout << "bcType = " << bcType << std::endl;
		fe_bc_[k][j] = Teuchos::rcp(new FE<Real > (bdryCells[i][j], basisPtrs_[0], bc_Cub, sideId));
		std::cout << "PDE:construct_FE_objects: fe_bc[" << k << "][" << j << "_ success" << std::endl;
		side_cub_points_physical_[k][j] = fe_bc_[k][j]->cubPts();
	      }
	  }
	if (bcType > 0)
	  k++;
      } 
  }
  
  
  void get_local_dofs_for_bc_sides()
  {
    local_dofs_on_sides_ = fe_vol_->getBoundaryDofs();
  }


  std::vector<Real > get_local_dof_coord(int local_cell_idx, int local_dof_idx)
  {
    std::vector<Real > coords(2);
    coords[0] = (*dof_points_physical_)(local_cell_idx, local_dof_idx, 0);
    coords[1] = (*dof_points_physical_)(local_cell_idx, local_dof_idx, 1);
    return coords;
  }
  

  Real neumann_flux_func(const std::vector<Real>  &x) 
  {
    return 1.0;
  }  
  
  Real dirichlet_data_func(const std::vector<Real>  &x) 
  {
    return 0.0;
  }  
  
  // specify neumann boundary condition by boundary object idx in the objects vector
  void evaluate_g (Teuchos::RCP<Intrepid::FieldContainer<Real> > &g_cub_points, int bc_seg_i, int bc_sub_seg_j)
  {
    int bc_type = boundary_type_[bc_seg_i];
    int numcell = numCells_bc_[bc_seg_i][bc_sub_seg_j];
    g_cub_points->initialize(0.0);
    
    // if neumann condition, return g
    if(bc_type == 1) 	
      {
	int zero = 0;
	int one = 1;
	for (int k=0; k<numcell; k++) {
	  for (int l=0; l<numCubPerSide_; l++) {
	    std::vector<Real > x(2);
	    x[0] = (*(side_cub_points_physical_[fe_bc_idx_[bc_seg_i]][bc_sub_seg_j]))(k, l, zero);	
	    x[1] = (*(side_cub_points_physical_[fe_bc_idx_[bc_seg_i]][bc_sub_seg_j]))(k, l, one);	
	    Real val_g = neumann_flux_func(x);
	    (*g_cub_points)(k, l) = -1.0 * val_g;
	  }
	}
      } 
    
  }
  
  // compute U^4 given U on the boundary cubature points
  void evaluate_pow4_U(Teuchos::RCP<Intrepid::FieldContainer<Real> > &U4_cub_points, Intrepid::FieldContainer<Real> &U_at_rbc)
  {
    int c = U_at_rbc.dimension(0);
    int p = U_at_rbc.dimension(1);
    for (int i=0; i < c; i++) {
      for (int j=0; j < p; j++) {
	(*U4_cub_points)(i, j) = alpha_ * std::pow( U_at_rbc(i, j), 4);
      }
    }
  }
  
	
  // compute grad U^4 given U on the boundary cubature points
  void evaluate_grad_pow4_U(Teuchos::RCP< Intrepid::FieldContainer<Real > > &grad_U4_cub, Intrepid::FieldContainer<Real> &U_at_rbc)
  { 
    int c = U_at_rbc.dimension(0);
    int p = U_at_rbc.dimension(1);
    for (int i=0; i < c; i++) {
      for (int j=0; j < p; j++) {
	(*grad_U4_cub)(i, j) = 4 * alpha_ * std::pow( U_at_rbc(i, j), 3);
      }
    }
  }
  
  
  // K at cubature points, input T has to be on the cubacture
  // T^2 + 1
  void evaluate_K (Teuchos::RCP< Intrepid::FieldContainer<Real > > &K_cub, Intrepid::FieldContainer<Real > &U)
  {    
    for (int i=0; i < numCells_vol_; i++) {
      for (int j=0; j < numCubPerCell_; j++) {
	(*K_cub)(i, j) = U(i, j) * U(i, j) + 1;
      }
    }
  } 
  
  // grad K at cubature points
  // 2 * T
  void evaluate_gradK (Teuchos::RCP<Intrepid::FieldContainer<Real> > &gradK_cub, Intrepid::FieldContainer<Real> &U)
  {
    for (int i=0; i < numCells_vol_; i++) {
      for (int j=0; j < numCubPerCell_; j++) {
	(*gradK_cub)(i, j) = 2 * U(i, j);
      }
    }
  } 


// not sure correct 
  void multiply_gradK_N_gradU (Teuchos::RCP< Intrepid::FieldContainer<Real> > &gradK_N_gradU, Intrepid::FieldContainer<Real> &gradK_cub, Intrepid::FieldContainer<Real> &N_cub, Intrepid::FieldContainer<Real> &gradU_cub)
  {
    for (int i=0; i < numCells_vol_; i++) {
      for (int j=0; j < lfs_; j++) {
	for (int k=0; k<numCubPerCell_; k++) {
	  for (int l=0; l<spaceDim_; l++) {
	    (*gradK_N_gradU)(i, j, k, l) = gradK_cub(i, k) * N_cub(i, j, k) * gradU_cub(i, k, l);
	  }
	}
      }
    }
  }
  
/////////////////////////////////


  void residual(Teuchos::RCP<Intrepid::FieldContainer<Real> > & res,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff) {
	
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_vol_, numCubPerCell_, spaceDim_));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_vol_, numCubPerCell_));
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valZ_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_vol_, numCubPerCell_));

    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    fe_vol_->evaluateValue(valZ_eval, z_coeff);
    Teuchos::RCP< Intrepid::FieldContainer<Real > > K_cub = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCells_vol_, numCubPerCell_));
    evaluate_K(K_cub,*valU_eval);
    
    //Intrepid::FieldContainer<Real > K_gradNdetJ = multiply_K_gradNdetJ (K_cub_points, 
    //								      *(fe_vol_->gradNdetJ()));
    res = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_vol_, lfs_));
    Teuchos::RCP<Intrepid::FieldContainer<Real > > K_gradU = Teuchos::rcp( new Intrepid::FieldContainer<Real >(numCells_vol_, numCubPerCell_, spaceDim_));
    Intrepid::FunctionSpaceTools::scalarMultiplyDataData<Real >(*K_gradU, *K_cub, *gradU_eval);	
    Intrepid::FunctionSpaceTools::integrate<Real >(*res, *K_gradU, *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);

    Intrepid::RealSpaceTools<Real >::scale(*valZ_eval, -1.0);
    Intrepid::FunctionSpaceTools::integrate<Real >(*res, *valZ_eval, *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, true);
	
    add_BC_terms_to_residual(res, u_coeff);
  }


  Teuchos::RCP<Intrepid::FieldContainer<Real > > get_boundary_u_coeff(const Intrepid::FieldContainer<Real> & cell_u_coeff, int bc_seg_i, int bc_sub_seg_j)
  {
    std::vector<int> bc_cell_local_id = bc_cell_local_id_[bc_seg_i][bc_sub_seg_j];
    int numcell = numCells_bc_[bc_seg_i][bc_sub_seg_j];
    
    Teuchos::RCP<Intrepid::FieldContainer<Real > > boundary_u_coeff = 
      Teuchos::rcp(new Intrepid::FieldContainer<Real > (numcell, lfs_));
    for (int i=0; i<numcell; i++) {
      for (int j=0; j<lfs_; j++){
	(*boundary_u_coeff)(i, j) = cell_u_coeff(bc_cell_local_id[i], j);
      }
    }
    return boundary_u_coeff;
  }
  
	
  void add_BC_terms_to_residual(Teuchos::RCP<Intrepid::FieldContainer<Real > > & res, const Teuchos::RCP< const Intrepid::FieldContainer<Real > > & u_coeff)
  {
    for (int i=1; i<n_bc_segments_; i++)
      {
	int bc_type = boundary_type_[i];
	
	if(bc_type==1) // Neumann
	  {
	    for (int j=0; j<n_bc_sub_segments_[i]; j++)
	      {
		int numCells = numCells_bc_[i][j];
		Teuchos::RCP< Intrepid::FieldContainer<Real > > g_nbc
                  = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCells_bc_[i][j], numCubPerSide_));
		evaluate_g(g_nbc, i, j);
                Teuchos::RCP<Intrepid::FieldContainer<Real> > neumannRes
                  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells,lfs_));
		Intrepid::FunctionSpaceTools::integrate<Real> (*neumannRes,
                                                               *g_nbc,
                                                               *(fe_bc_[fe_bc_idx_[i]][j]->NdetJ()),
                                                               Intrepid::COMP_CPP, true);
		// Add Robin residual to volume residual
		for (int k = 0; k < numCells; ++k) {
		  int cidx = bc_cell_local_id_[i][j][k];
                  for (int l = 0; l < lfs_; ++l) { 
                    (*res)(cidx,l) += (*neumannRes)(k,l);
                  }
		}
	      }
	  }
	else if(bc_type==2) // Robin
	  {
	    for (int j=0; j<n_bc_sub_segments_[i]; j++)
	      {
		int numCells = numCells_bc_[i][j];
		Teuchos::RCP<Intrepid::FieldContainer<Real > > bc_u_coeff
                  = get_boundary_u_coeff(*u_coeff, i, j);
		Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bc
                  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_bc_[i][j], numCubPerSide_));
		fe_bc_[fe_bc_idx_[i]][j]->evaluateValue(valU_eval_bc, bc_u_coeff);
		Teuchos::RCP< Intrepid::FieldContainer<Real> > valU_pow4
                  = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCells_bc_[i][j], numCubPerSide_));
		evaluate_pow4_U(valU_pow4,*valU_eval_bc);
                Teuchos::RCP<Intrepid::FieldContainer<Real> > robinRes
                  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells,lfs_));
		Intrepid::FunctionSpaceTools::integrate<Real>(*robinRes,
                                                              *valU_pow4,
                                                              *(fe_bc_[fe_bc_idx_[i]][j]->NdetJ()),
                                                              Intrepid::COMP_CPP, false);
		// Add Robin residual to volume residual
		for (int k = 0; k < numCells; ++k) {
		  int cidx = bc_cell_local_id_[i][j][k];
                  for (int l = 0; l < lfs_; ++l) { 
                    (*res)(cidx,l) += (*robinRes)(k,l);
                  }
		}
	      }
	  }
      }
      // Apply Dirichlet conditions
	    for (int j=0; j<n_bc_sub_segments_[0]; j++)
	      {
		int numCells = numCells_bc_[0][j];
		std::vector<int> fidx = local_dofs_on_sides_[sideIds_[0][j]];
		int numBdryDofs = fidx.size();
		for (int k = 0; k < numCells; ++k) {
		  int cidx = bc_cell_local_id_[0][j][k];
		  for (int l = 0; l < numBdryDofs; ++l) {
		    std::vector<Real > x = get_local_dof_coord(cidx, fidx[l]);
		    Real dbc_data = dirichlet_data_func(x);
          (*res)(cidx,fidx[l]) = (*u_coeff)(cidx,fidx[l]) - dbc_data;
        }
      }
    }
  }
  
  void Jacobian_1(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff) {
    
    Teuchos::RCP<Intrepid::FieldContainer<Real> > gradU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_vol_, numCubPerCell_, spaceDim_));
    fe_vol_->evaluateGradient(gradU_eval, u_coeff);
    Teuchos::RCP<Intrepid::FieldContainer<Real> > valU_eval =
      Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_vol_, numCubPerCell_));
    fe_vol_->evaluateValue(valU_eval, u_coeff);
    
    Teuchos::RCP< Intrepid::FieldContainer<Real > > K_cub
      = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCells_vol_, numCubPerCell_));
    evaluate_K( K_cub, *valU_eval);
    Teuchos::RCP< Intrepid::FieldContainer<Real > > gradK_cub
      = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCells_vol_, numCubPerCell_));
    evaluate_gradK( gradK_cub, *valU_eval );
    
    Teuchos::RCP< Intrepid::FieldContainer<Real > > gradK_N_gradU
      = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCells_vol_, lfs_, numCubPerCell_, spaceDim_));
    multiply_gradK_N_gradU(gradK_N_gradU, *gradK_cub, *(fe_vol_->N()), *(gradU_eval));
    
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_vol_, lfs_, lfs_));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac, *gradK_N_gradU, *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, false);
    
    //Intrepid::FieldContainer<Real > K_gradNdetJ = multiply_K_gradNdetJ(K_cub, *(fe_vol_->gradN()));
    Teuchos::RCP< Intrepid::FieldContainer<Real > > K_gradN
      = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_vol_, lfs_, numCubPerCell_, spaceDim_));
    
    Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*K_gradN, *K_cub, *(fe_vol_->gradN()));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac, *K_gradN, *(fe_vol_->gradNdetJ()), Intrepid::COMP_CPP, true);
    
    for (int i=1; i<n_bc_segments_; i++)
      {
	int bc_type = boundary_type_[i];
	if(bc_type==2) //robin
	  {
	    for (int j=0; j<n_bc_sub_segments_[i]; j++)
	      {
		int numCells = numCells_bc_[i][j];
		Teuchos::RCP<Intrepid::FieldContainer<Real > > bc_u_coeff
                  = get_boundary_u_coeff(*u_coeff, i, j);
		Teuchos::RCP<Intrepid::FieldContainer<Real > > valU_eval_bc
                  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_bc_[i][j], numCubPerSide_));
		fe_bc_[fe_bc_idx_[i]][j]->evaluateValue(valU_eval_bc, bc_u_coeff);
		
		int numcell = valU_eval_bc->dimension(0);
		Teuchos::RCP< Intrepid::FieldContainer<Real > > valU_4pow3
                  = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numcell, numCubPerSide_));
		evaluate_grad_pow4_U(valU_4pow3, *valU_eval_bc);
		//Intrepid::FieldContainer<Real> valU_4pow3_N = 
		//			multiply_valU_w_N_rbc(valU_4pow3, *(fe_bc_[i][j]->N()));
		
		Teuchos::RCP< Intrepid::FieldContainer<Real> > valU_4pow3_N
                  = Teuchos::rcp( new Intrepid::FieldContainer<Real>(numCells_bc_[i][j], lfs_, numCubPerSide_));
		Intrepid::FunctionSpaceTools::scalarMultiplyDataField<Real>(*valU_4pow3_N, *valU_4pow3, *(fe_bc_[fe_bc_idx_[i]][j]->N()));
                Teuchos::RCP<Intrepid::FieldContainer<Real> > robinJac
                  = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells,lfs_,lfs_));
		Intrepid::FunctionSpaceTools::integrate<Real>(*robinJac,
                                                              *valU_4pow3_N,
                                                              *(fe_bc_[fe_bc_idx_[i]][j]->NdetJ()),
                                                              Intrepid::COMP_CPP, true);
          // Add Robin Jacobian to volume Jacobian
          for (int k = 0; k < numCells; ++k) {
            int cidx = bc_cell_local_id_[i][j][k];
            for (int l = 0; l < lfs_; ++l) { 
              for (int m = 0; m < lfs_; ++m) {
                (*jac)(cidx,l,m) += (*robinJac)(k,l,m);
              }
            }
          }
        }
      }
    }
    // Apply Dirichlet conditions
    for (int j=0; j<n_bc_sub_segments_[0]; ++j) {
      int numCells = numCells_bc_[0][j];
      std::vector<int> fidx = local_dofs_on_sides_[sideIds_[0][j]];
      int numBdryDofs = fidx.size();
      for (int k = 0; k < numCells; ++k) {
        int cidx = bc_cell_local_id_[0][j][k];
        for (int l = 0; l < numBdryDofs; ++l) {
          // Modify the local jacobian matrix
          for(int n=0; n<lfs_; ++n) {
            (*jac)(cidx, fidx[l], n) = 0.0;
          }
          (*jac)(cidx, fidx[l], fidx[l]) = 1.0;
        }
      }
    }
  }
  
  
/*
  Intrepid::FieldContainer<Real > multiply_valU_w_N_rbc ( int boundary_object_idx,
							  Intrepid::FieldContainer<Real > rbc_val_U_4pow3, 
							  Intrepid::FieldContainer<Real > rbc_N )
 	{
		int lfs = rbc_N.dimension(1);
		int numCell = numCells_bc_[boundary_object_idx];

		Intrepid::FieldContainer<Real > rbc_valU_N
			= new Intrepid::FieldContainer<Real > (numCell, lfs, numCubPerSide_);
		for (int i=0; i<numCell; i++) {
			for (int j=0; j<lfs; j++) {
				for (int k=0; k<numCubPerSide_; k++) {
					rbc_valU_N(i, j, k) = rbc_val_U_4pow3(i, k) * rbc_N(i, j, k);	
				}
			}
		}
		return rbc_valU_N;
	}
*/

  void Jacobian_2(Teuchos::RCP<Intrepid::FieldContainer<Real> > & jac,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & u_coeff,
                  const Teuchos::RCP<const Intrepid::FieldContainer<Real> > & z_coeff) {
    Teuchos::RCP< Intrepid::FieldContainer<Real > > negative_N = Teuchos::rcp( new Intrepid::FieldContainer<Real >(numCells_vol_, lfs_, numCubPerCell_));
    Intrepid::RealSpaceTools<Real >::scale(*negative_N, *(fe_vol_->N()), -1.0);
    jac = Teuchos::rcp(new Intrepid::FieldContainer<Real>(numCells_vol_, lfs_, lfs_));
    Intrepid::FunctionSpaceTools::integrate<Real>(*jac, *negative_N, *(fe_vol_->NdetJ()), Intrepid::COMP_CPP, false);
	    for (int j=0; j<n_bc_sub_segments_[0]; j++)
	      {
		// Apply Dirichlet conditions
		int numCells = numCells_bc_[0][j];
		std::vector<int> fidx = local_dofs_on_sides_[sideIds_[0][j]];
		int numBdryDofs = fidx.size();
		for (int k = 0; k < numCells; ++k) {
		  int cidx = bc_cell_local_id_[0][j][k];
		  for (int l = 0; l < numBdryDofs; ++l) {
		    // Modify the local jacobian matrix
		    for(int n=0; n<lfs_; ++n) {
		      (*jac)(cidx, fidx[l], n) = 0.0;
		    }
		  }
		}
	      }
  }


  void Hessian_11(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & u_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & z_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::NotImplemented(">>> Hessian_11 not implemented.");
  }

  void Hessian_12(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & u_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & z_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::NotImplemented(">>> Hessian_12 not implemented.");
  }

  void Hessian_21(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & u_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & z_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::NotImplemented(">>> Hessian_21 not implemented.");
  }

  void Hessian_22(Teuchos::RCP<Intrepid::FieldContainer<Real> > & hess,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & u_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & z_coeff,
                  Teuchos::RCP<Intrepid::FieldContainer<Real> > & l_coeff) {
    throw Exception::NotImplemented(">>> Hessian_22 not implemented.");
  }
  

  const Teuchos::RCP<FE<Real> > getFE_VOL(void) const {
    return fe_vol_;
  }

}; // PDE_stefan_boltzmann

#endif
