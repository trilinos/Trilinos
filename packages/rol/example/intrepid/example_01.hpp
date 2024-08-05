// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "discretization.hpp"
#include "coefficient.hpp"

// ROL Includes
#include "ROL_StdVector.hpp"
#include "ROL_Objective_SimOpt.hpp"
#include "ROL_Constraint_SimOpt.hpp"

// Teuchos Includes for Linear Algebra
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseSolver.hpp"



// Quadratic tracking objective class

using namespace ROL;
template<class Real>
class TrackingObjective : public Objective_SimOpt<Real> {

  template <typename T> using ROL::Ptr  = ROL::Ptr<T>;
  template <typename T> using FC   = Intrepid::FieldContainer<T>;            

  typedef Vector<Real>                  V;
  typedef StdVector<Real>               SV;
  typedef std::vector<Real>             vec;

  private:
    ROL::Ptr<Discretization<Real>> disc_;

    int numCells_;              // Number of cells (elements)
    int numCubPts_;             // Number of cubature points per cells
    int numFields_;             // Number of basis functions per cell
    int spaceDim_;              // Number of spatial dimensions (currently 1)
    int nDoF_;                  // Total number of degrees of freedom

    Real gamma_;

    ROL::Ptr<FC<Real>> massMatrices_;     

    ROL::Ptr<V> target_;


    void applyMass(V &Mv, const V &v) {

      

      ROL::Ptr<vec> Mvp = dynamic_cast<SV&>(Mv).getVector(); 
      ROL::Ptr<const vec> vp = dynamic_cast<const SV&>(v).getVector();

      for(int cell=0;cell<numCells_;++cell) {
        for(int rfield=0;rfield<numFields_;++rfield) {
          int i = cell*(numFields_-1) + rfield;
          for(int cfield=0;cfield<numFields_;++cfield) {
            int j = cell*(numFields_-1) + cfield;
            (*Mvp)[i] += (*massMatrices_)(cell,rfield,cfield)*(*vp)[j];  
          }
        }
      }
    }       

  public:    

    TrackingObjective(ROL::Ptr<Discretization<Real>> disc, 
                      const ROL::Ptr<V> &target, 
                      const Real &gamma) : 
      disc_(disc), 
      numCells_(disc_->getNumCells()),
      numCubPts_(disc->getNumCubPts()),
      numFields_(disc->getNumFields()),
      spaceDim_(disc->getSpaceDim()),
      nDoF_(numCells_*(numFields_-1)+1),
      gamma_(gamma), 
      massMatrices_(disc->getMassMatrices()), 
      target_(target)  {}

    Real value(const V &u, const V &z, Real &tol) {

      ROL::Ptr<V> err = u.clone();
      err->set(u);
      err->axpy(-1.0,*target_);
      ROL::Ptr<V> Merr = u.clone();
      applyMass(*Merr,*err);
      
      ROL::Ptr<V> y = z.clone();
      y->set(z);
      ROL::Ptr<V> My = z.clone();   
      applyMass(*My,*y); 

      Real J = 0.5*(Merr->dot(*err)+gamma_*My->dot(*y));
      return J;       

    } // value()

    void gradient_1( V &g, const V &u, const V &z, Real &tol ) {

      ROL::Ptr<V> err = u.clone();
      err->set(u);
      err->axpy(-1.0,*target_);
      applyMass(g,*err);
     
      
    } // gradient_1()

    void gradient_2( V &g, const V &u, const V &z, Real &tol ) {

      ROL::Ptr<V> y = z.clone();
      y->set(z);
      applyMass(g,*y); 
      g.scale(gamma_);

    } // gradient_2()

    void hessVec_11( V &hv, const V &v, const V &u, const V &z, Real &tol) {

      ROL::Ptr<V> y = v.clone();
      y->set(v);
      applyMass(hv,*y);
         
    } // hessVec_11()  

    void hessVec_12( V &hv, const V &v, const V &u, const V &z, Real &tol) {

      hv.zero();

    } // hessVec_12()  

    void hessVec_21( V &hv, const V &v, const V &u, const V &z, Real &tol) {

      hv.zero();

    } // hessVec_21()  

    void hessVec_22( V &hv, const V &v, const V &u, const V &z, Real &tol) {

      ROL::Ptr<V> y = v.clone();
      y->set(v);
      applyMass(hv,*y);
      hv.scale(gamma_);
 
    } // hessVec_22()  

};



// BVP Constraint class 

template<class Real> 
class BVPConstraint : public Constraint_SimOpt<Real> {

  template <typename T> using ROL::Ptr  = ROL::Ptr<T>;
  template <typename T> using FC   = Intrepid::FieldContainer<T>;            

  typedef Teuchos::SerialDenseMatrix<int,Real> Matrix;
  typedef Teuchos::SerialDenseSolver<int,Real> Solver;

  typedef Intrepid::FunctionSpaceTools  FST;

  typedef Vector<Real>                  V;
  typedef StdVector<Real>               SV;
  typedef std::vector<Real>             vec;

  private:
    ROL::Ptr<Discretization<Real>> disc_;

    int numCells_;              // Number of cells (elements)
    int numCubPts_;             // Number of cubature points per cells
    int numFields_;             // Number of basis functions per cell
    int spaceDim_;              // Number of spatial dimensions (currently 1)
    int nDoF_;                  // Total number of degrees of freedom

    ROL::Ptr<FC<Real>> x_cub_;       // Physical cubature points
    ROL::Ptr<FC<Real>> tranVals_;    // Transformed values of basis functions
    ROL::Ptr<FC<Real>> tranGrad_;    // Transformed gradients of basis functions
    ROL::Ptr<FC<Real>> wtdTranVals_; // Weighted transformed values of basis functions
    ROL::Ptr<FC<Real>> wtdTranGrad_; // Weighted transformed gradients of basis functions

    ROL::Ptr<Matrix> Ju_;            // Constraint Jacobian w.r.t. u
    ROL::Ptr<Matrix> Jz_;            // Constraint Jacobian w.r.t. z
    ROL::Ptr<Matrix> JuFactors_;     // LU factors of the Sim Jacobian 

    vec dif_param_;             // Parameters passed to coefficient functions. Currently unused
    vec adv_param_;
    vec rea_param_; 

    enum var { sim, opt };      // Index for simulation and optimization variable

    // Write ROL vector into a one-column Teuchos::SerialDenseMatrix
    void vec2mat(ROL::Ptr<Matrix> &m, const V &v) {

      
      ROL::Ptr<const vec> vp = dynamic_cast<const SV&>(v).getVector();

      for(int i=0;i<nDoF_;++i) {
        (*m)(i,0) = (*vp)[i];  
      }       
    }

    // Write a one-column Teuchos::SerialDenseMatrix into a ROL vector 
    void mat2vec(V &v, const ROL::Ptr<Matrix> &m) {

       
       using ROL::constPtrCast;

       ROL::Ptr<vec> vp = (dynamic_cast<SV&>(v)).getVector(); 
       
       for(int i=0;i<nDoF_;++i) {
         (*vp)[i] = (*m)(i,0);   
       }
    }

    // Gather a ROL vector into a Intrepid Field Container
    template<class ScalarT>
    void gather(FC<ScalarT> &fc, const V& v) {

      
      ROL::Ptr<const vec> vp = dynamic_cast<const SV&>(v).getVector();

      for(int cell=0;cell<numCells_;++cell) {
        for(int field=0;field<numFields_;++field) {
          int i = cell*(numFields_-1) + field;
          fc(cell,field) = (*vp)[i]; 
        }
      }
    }

    // Compute the residual given u and z
    template<class ScalarT>
    void evaluate_res(FC<ScalarT> &c_fc, FC<ScalarT> &u_fc, FC<ScalarT> &z_fc) {

       

       ROL::Ptr<Coefficient<Real,ScalarT>> coeff = ROL::makePtr<ExampleCoefficient<Real,ScalarT>>();

       // Evaluate on the cubature points 
       FC<ScalarT> u_vals_cub(numCells_,numCubPts_);
       FC<ScalarT> z_vals_cub(numCells_,numCubPts_);
       FC<ScalarT> u_grad_cub(numCells_,numCubPts_,spaceDim_);       

       FST::evaluate<ScalarT>(u_vals_cub,u_fc,*tranVals_);
       FST::evaluate<ScalarT>(z_vals_cub,z_fc,*tranVals_);
       FST::evaluate<ScalarT>(u_grad_cub,u_fc,*tranGrad_); 

       // Evaluate terms on the cubature points
       FC<ScalarT> react_cub(numCells_,numCubPts_); 
       FC<ScalarT> advec_cub(numCells_,numCubPts_,spaceDim_); 
       FC<ScalarT> diff_cub(numCells_,numCubPts_); 
 
       coeff->reaction(react_cub, *x_cub_,u_vals_cub,z_vals_cub,rea_param_);
       coeff->advection(advec_cub,*x_cub_,u_vals_cub,z_vals_cub,adv_param_);
       coeff->diffusion(diff_cub, *x_cub_,u_vals_cub,z_vals_cub,dif_param_);

       FC<ScalarT> advec_term(numCells_,numCubPts_);
       FC<ScalarT> diff_term(numCells_,numCubPts_,spaceDim_);
       
       FST::scalarMultiplyDataData<ScalarT>(diff_term,diff_cub,u_grad_cub); 
       FST::dotMultiplyDataData<ScalarT>(advec_term,advec_cub,u_grad_cub);
 
       // Add terms to residual
       c_fc.initialize();
       FST::integrate<ScalarT>(c_fc,diff_term, *wtdTranGrad_,Intrepid::COMP_CPP,false);  
       FST::integrate<ScalarT>(c_fc,advec_term,*wtdTranVals_,Intrepid::COMP_CPP,true);
       FST::integrate<ScalarT>(c_fc,react_cub, *wtdTranVals_,Intrepid::COMP_CPP,true);

    }

    void applyJac(V &jv, const V &v, var comp, bool transpose = false) {

       ROL::Ptr<Matrix> J = (comp==sim) ? Ju_ : Jz_;

       ROL::Ptr<Matrix> vmat = ROL::makePtr<Matrix>(nDoF_,1,true);
       ROL::Ptr<Matrix> jvmat = ROL::makePtr<Matrix>(nDoF_,1,true);

       vec2mat(vmat, v);

       if(transpose) {
         jvmat->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*J,*vmat,0);
       }
       else {
         jvmat->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*J,*vmat,0);
       }
       mat2vec(jv,jvmat);
    }

    void applyAdjointHessian(V &ahwv, const V &w, const V &v, 
                             const V  &u, const V &z, var one, var two )  {

      

      typedef Sacado::Fad::SFad<Real,1> SFad;
      typedef Sacado::Fad::DFad<SFad>   DSFad;

      ROL::Ptr<vec> ahwvp = dynamic_cast<SV&>(ahwv).getVector();

      std::fill(ahwvp->begin(),ahwvp->end(),0.0);

      ROL::Ptr<const vec> vp = (dynamic_cast<const SV&>(v)).getVector();
      ROL::Ptr<const vec> wp = (dynamic_cast<const SV&>(w)).getVector();
      ROL::Ptr<const vec> up = (dynamic_cast<const SV&>(u)).getVector();
      ROL::Ptr<const vec> zp = (dynamic_cast<const SV&>(z)).getVector();

      FC<DSFad> u_fc(numCells_,numFields_);
      FC<DSFad> z_fc(numCells_,numFields_);

      if(one == sim && two == sim) { // H11 block
        for(int cell=0;cell<numCells_;++cell) {
          for(int field=0;field<numFields_;++field) {
            int i = cell*(numFields_-1) + field; 
            SFad temp(1,(*up)[i]);
            temp.fastAccessDx(0) = (*vp)[i]; 
            u_fc(cell,field) = DSFad(numFields_,field,temp);
            z_fc(cell,field) = (*zp)[i];
          }
        } 
      }
      else if(one == sim && two == opt) { // H12 block
        for(int cell=0;cell<numCells_;++cell) {
          for(int field=0;field<numFields_;++field) {
            int i = cell*(numFields_-1) + field; 
            SFad temp(1,(*up)[i]);
            temp.fastAccessDx(0) = (*vp)[i]; 
            u_fc(cell,field) = temp;
            z_fc(cell,field) = DSFad(numFields_,field,(*zp)[i]);
          }
        } 
      }
      else if(one == opt && two == sim) { // H21 block
        for(int cell=0;cell<numCells_;++cell) {
          for(int field=0;field<numFields_;++field) {
            int i = cell*(numFields_-1) + field; 
            SFad temp(1,(*zp)[i]);
            temp.fastAccessDx(0) = (*vp)[i]; 
            z_fc(cell,field) = temp;
            u_fc(cell,field) = DSFad(numFields_,field,(*up)[i]);
          }
        } 
      }
      else { // H22 block
        for(int cell=0;cell<numCells_;++cell) {
          for(int field=0;field<numFields_;++field) {
            int i = cell*(numFields_-1) + field; 
            SFad temp(1,(*zp)[i]);
            temp.fastAccessDx(0) = (*vp)[i]; 
            z_fc(cell,field) = DSFad(numFields_,field,temp);
            u_fc(cell,field) = (*up)[i];
          }
        } 
      }

      FC<DSFad> c_fc(numCells_,numFields_);         

      evaluate_res(c_fc,u_fc,z_fc);

     // Compute the cellwise dot product of the constraint value and the Lagrange multiply
     for(int cell=0;cell<numCells_;++cell) {
       DSFad wdotc(SFad(0.0));
       for(int field=0;field<numFields_;++field) {
         int i = cell*(numFields_-1) + field; 
         wdotc += c_fc(cell,field)*(*wp)[i];
       }
       for(int field=0;field<numFields_;++field) {
         int i = cell*(numFields_-1) + field; 
         (*ahwvp)[i] += wdotc.dx(field).dx(0);
       }
      } 
    }


  public:

    // Constructor
    BVPConstraint(const ROL::Ptr<Discretization<Real>> disc) : 
      disc_(disc),
      numCells_(disc->getNumCells()),
      numCubPts_(disc->getNumCubPts()),
      numFields_(disc->getNumFields()),
      spaceDim_(disc->getSpaceDim()),
      nDoF_(numCells_*(numFields_-1)+1),
      x_cub_(disc->getPhysCubPts()),
      tranVals_(disc_->getTransformedVals()),
      tranGrad_(disc_->getTransformedGrad()),
      wtdTranVals_(disc_->getWeightedTransformedVals()),
      wtdTranGrad_(disc_->getWeightedTransformedGrad()),
      Ju_(ROL::makePtr<Matrix>(nDoF_,nDoF_,true)),
      Jz_(ROL::makePtr<Matrix>(nDoF_,nDoF_,true)),
      JuFactors_(ROL::makePtr<Matrix>(nDoF_,nDoF_,true)) {
     
    }


    void update( const V& u, const V &z, bool flag, int iter = -1 ) {

        using Sacado::Fad::DFad;        

        
        typedef DFad<Real> FadType;

        ROL::Ptr<const vec> up = dynamic_cast<const SV&>(u).getVector();
        ROL::Ptr<const vec> zp = dynamic_cast<const SV&>(z).getVector();

        FC<FadType> u_fc1(numCells_,numFields_); 
        FC<FadType> z_fc1(numCells_,numFields_);
        FC<FadType> u_fc2(numCells_,numFields_); 
        FC<FadType> z_fc2(numCells_,numFields_);

        gather(z_fc1,z);
        gather(u_fc2,u);

        // Fill in field containers from vectors (gather)
        for(int cell=0;cell<numCells_;++cell) {
          for(int field=0;field<numFields_;++field) {
            int i = cell*(numFields_-1) + field;
            u_fc1(cell,field) = FadType(numFields_,field,(*up)[i]);
            z_fc2(cell,field) = FadType(numFields_,field,(*zp)[i]);
          } 
        }       

        FC<FadType> c_fc1(numCells_,numFields_); 
        FC<FadType> c_fc2(numCells_,numFields_); 

        evaluate_res(c_fc1,u_fc1,z_fc1); 
        evaluate_res(c_fc2,u_fc2,z_fc2); 

        Ju_->putScalar(0.0); // Zero out the matrix
        Jz_->putScalar(0.0); // Zero out the matrix

        for(int cell=0;cell<numCells_;++cell) {
          for(int rfield=0;rfield<numFields_;++rfield) {
            int row = cell*(numFields_-1) + rfield;
            for(int cfield=0;cfield<numFields_;++cfield) {
              int col = cell*(numFields_-1) + cfield;
              (*Ju_)(row,col) += c_fc1(cell,rfield).dx(cfield);
              (*Jz_)(row,col) += c_fc2(cell,rfield).dx(cfield);
            }
          } 
        }       

        JuFactors_->assign(*Ju_);
    }


    void value(V &c, const V &u, const V &z, Real &tol=0) {

      
      // Downcast and extract ROL::Ptrs to std::vectors
      ROL::Ptr<vec> cp = dynamic_cast<SV&>(c).getVector(); 

      std::fill(cp->begin(),cp->end(),0.0); 

      FC<Real> u_fc(numCells_,numFields_); 
      FC<Real> z_fc(numCells_,numFields_);

      gather(u_fc,u);
      gather(z_fc,z);

      FC<Real> c_fc(numCells_,numFields_);

      evaluate_res(c_fc,u_fc,z_fc);
 
      // Scatter residual back into ROL vector
      for(int cell=0;cell<numCells_;++cell) {
       for(int field=0;field<numFields_;++field) {
          int i = cell*(numFields_-1) + field;
          (*cp)[i] += c_fc(cell,field);
        } 
      }       
   } // value()


    void applyJacobian_1(V &jv, const V &v, const V &u, const V &z, Real &tol) {
      applyJac(jv,v,sim,false);

    } // applyJacobian_1()

   
    void applyJacobian_2(V &jv, const V &v, const V &u, const V &z, Real &tol) {
      applyJac(jv,v,opt,false);

    }  // applyJacobian_2()

      
    void applyAdjointJacobian_1(V &jv, const V &v, const V &u, const V &z, Real &tol) {
      applyJac(jv,v,sim,true);

    } // applyAdjointJacobian_1()

   
    void applyAdjointJacobian_2(V &jv, const V &v, const V &u, const V &z, Real &tol) {
      applyJac(jv,v,opt,true);

    }  // applyAdjointJacobian_2()


    void applyInverseJacobian_1(V &ijv, const V &v, const V &u, const V &z, Real &tol) {
      Teuchos::SerialDenseSolver<int,Real> solver;

      solver.setMatrix(JuFactors_);
      solver.factorWithEquilibration(true);
      solver.factor();

      ROL::Ptr<Matrix> rhs = ROL::makePtr<Matrix>(nDoF_,1);
      ROL::Ptr<Matrix> sol = ROL::makePtr<Matrix>(nDoF_,1);
       
      vec2mat(rhs,v);             // Write the ROL vector into the rhs 
      solver.setVectors(sol,rhs); 
      solver.solveWithTranspose(true);
      solver.solve();             // Solve the system
      solver.solveWithTranspose(false);
      mat2vec(ijv,sol);           // Write the solution into the ROL vector
      

    } // applyInverseJacobian_1()


    void applyInverseAdjointJacobian_1(V &iajv, const V &v, const V &u, const V &z, Real &tol) {
      Teuchos::SerialDenseSolver<int,Real> solver;

      solver.setMatrix(JuFactors_);
      solver.factorWithEquilibration(true);
      solver.factor();

      ROL::Ptr<Matrix> rhs = ROL::makePtr<Matrix>(nDoF_,1);
      ROL::Ptr<Matrix> sol = ROL::makePtr<Matrix>(nDoF_,1);
       
      vec2mat(rhs,v);                     // Write the ROL vector into the rhs 
      solver.setVectors(sol,rhs); 
      solver.solveWithTranspose(true);
      solver.solve();                     // Solve the system
      solver.solveWithTranspose(false);
      mat2vec(iajv,sol);                  // Write the solution into the ROL vector


    } // applyInverseAdjointJacobian_1()

    

    void applyAdjointHessian_11(V &ahwv, const V &w, const V &v, 
                                const V  &u, const V &z, Real &tol ) {
      applyAdjointHessian(ahwv, w, v, u, z, sim, sim);

    } // applyAdjointHessian_11()


    void applyAdjointHessian_12(V &ahwv, const V &w, const V &v, 
                                const V  &u, const V &z, Real &tol ) {
      applyAdjointHessian(ahwv, w, v, u, z, sim, opt);
     
    } // applyAdjointHessian_11()

    void applyAdjointHessian_21(V &ahwv, const V &w, const V &v, 
                                const V  &u, const V &z, Real &tol ) {
      applyAdjointHessian(ahwv, w, v, u, z, opt, sim);
    } // applyAdjointHessian_22()

    void applyAdjointHessian_22(V &ahwv, const V &w, const V &v, 
                                const V  &u, const V &z, Real &tol ) {
      applyAdjointHessian(ahwv, w, v, u, z, opt, opt);
    } // applyAdjointHessian_22()

};



