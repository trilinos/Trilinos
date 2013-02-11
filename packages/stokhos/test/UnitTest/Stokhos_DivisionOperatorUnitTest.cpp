// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Stokhos.hpp"
#include "Stokhos_UnitTestHelpers.hpp"
#include "Stokhos_DenseDirectDivisionExpansionStrategy.hpp"
#include "Stokhos_GMRESDivisionExpansionStrategy.hpp"
#include "Stokhos_CGDivisionExpansionStrategy.hpp"
#include "Stokhos_StandardStorage.hpp"
#include "Stokhos_OrthogPolyApprox.hpp"

namespace DivisionOperatorUnitTest {

  // Common setup for unit tests
  template <typename OrdinalType, typename ValueType>
  struct UnitTestSetup {
    ValueType rtol, atol;
    ValueType crtol, catol;
    OrdinalType sz;
    Teuchos::RCP<const Stokhos::CompletePolynomialBasis<OrdinalType,ValueType> > basis;
    Teuchos::RCP<const Stokhos::Quadrature<OrdinalType,ValueType> > quad;
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk, Cijk_linear;
    Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType> > exp, exp_linear;
    Teuchos::RCP< Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType> > qexp;
    Stokhos::OrthogPolyApprox<OrdinalType,ValueType> x, y, u, u2, cx, cu, cu2, sx, su, su2, c1;
    ValueType a;
    Teuchos::RCP< Stokhos::DenseDirectDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > direct_division_strategy;
      
    
    UnitTestSetup() {
      rtol = 1e-10;//4
      atol = 1e-10;//5
      crtol = 1e-12;
      catol = 1e-12;
      a = 3.1;
      const OrdinalType d = 2;
      const OrdinalType p = 7;
      
      // Create product basis
      Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<OrdinalType,ValueType> > > bases(d);
      for (OrdinalType i=0; i<d; i++)
	bases[i] = 
	  Teuchos::rcp(new Stokhos::LegendreBasis<OrdinalType,ValueType>(p));
      basis =
	Teuchos::rcp(new Stokhos::CompletePolynomialBasis<OrdinalType,ValueType>(bases));

      // Tensor product quadrature
      quad = 
	Teuchos::rcp(new Stokhos::TensorProductQuadrature<OrdinalType,ValueType>(basis));

      // Triple product tensor
      Cijk = basis->computeTripleProductTensor();
      Cijk_linear = basis->computeLinearTripleProductTensor();
      
      // Algebraic expansion
      exp = 
	Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType>(basis, Cijk, quad));
      exp_linear = 
	Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType>(basis, Cijk_linear, quad));

      // Quadrature expansion
      qexp = 
	Teuchos::rcp(new Stokhos::QuadOrthogPolyExpansion<OrdinalType,ValueType>(basis, Cijk, quad));
      //Dense Direct Division Operator
      direct_division_strategy =
         Teuchos::rcp(new Stokhos::DenseDirectDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(basis, Cijk));

      
      // Create approximation
//      sz = basis->size();
      x.reset(basis);
      y.reset(basis);
      u.reset(basis); 
      u2.reset(basis);
      cx.reset(basis, 1);
      x.term(0, 0) = 1.0;
//      y.term(0, 0) = 1.0:
      cx.term(0, 0) = a;
      cu.reset(basis);
//      cu2.reset(basis, 1);
//      sx.reset(basis, d+1);
//      su.reset(basis, d+1);
//      su2.reset(basis, d+1);
      for (OrdinalType i=0; i<d; i++) {
	x.term(i, 1) = 1.0;
//	y.term(i, 1) = 0.1;
      }
//      y.term(0, 0) = 2.0;
//      for (OrdinalType i=0; i<d; i++)
//	y.term(i, 1) = 0.25;
      
      c1.reset(basis);
      c1.term(0,0)=1;
      exp->exp(cu, x);

     }
};

  UnitTestSetup<int,double> setup;



  TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, CG_Divide ) {
    Teuchos::RCP< Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > cg_division_strategy =
       Teuchos::rcp(new Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 0, 100, 0, 0,1));
    cg_division_strategy->divide(setup.u, 1.0, setup.c1, setup.x, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.x, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2", 
				   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, CG_Jacobi_Divide ) {
    Teuchos::RCP< Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > cg_diag_division_strategy =
       Teuchos::rcp(new Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 2, 100, 0, 0,1));
    cg_diag_division_strategy->divide(setup.u, 1.0, setup.c1, setup.x, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.x, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, CG_SymGaussSeidel_Divide ) {
    Teuchos::RCP< Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > cg_jacobi_division_strategy =
       Teuchos::rcp(new Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 3, 100, 0, 0,1));
    cg_jacobi_division_strategy->divide(setup.u, 1.0, setup.c1, setup.x, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.x, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, CG_Schur_Divide ) {
    Teuchos::RCP< Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > cg_schur_division_strategy =
       Teuchos::rcp(new Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 0, 1e-12, 4, 100, 0, 0,1));
    cg_schur_division_strategy->divide(setup.u, 1.0, setup.c1, setup.x, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.x, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, CG_Nonlin_Divide ) {
    Teuchos::RCP< Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > cg_nonlin_division_strategy =
       Teuchos::rcp(new Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 0, 100, 0, 0,1));  
    cg_nonlin_division_strategy->divide(setup.u, 1.0, setup.c1, setup.cu, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.cu, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }
   TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, CG_Nonlin_Jacobi_Divide ) {
    Teuchos::RCP< Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > cg_nonlin_division_strategy =
       Teuchos::rcp(new Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 2, 100, 0, 0,1));
    cg_nonlin_division_strategy->divide(setup.u, 1.0, setup.c1, setup.cu, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.cu, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }
   TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, CG_Nonlin_SymGaussSeidel_Divide ) {
    Teuchos::RCP< Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > cg_nonlin_division_strategy =
       Teuchos::rcp(new Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 3, 100, 0, 0,1));
    cg_nonlin_division_strategy->divide(setup.u, 1.0, setup.c1, setup.cu, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.cu, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }

   TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, CG_Nonlin_Schur_Divide ) {
    Teuchos::RCP< Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > cg_nonlin_division_strategy =
       Teuchos::rcp(new Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 0, 1e-12, 4, 100, 0, 0,1));
    cg_nonlin_division_strategy->divide(setup.u, 1.0, setup.c1, setup.cu, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.cu, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }
   TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, CG_Nonlin_Schur_linearprec_Divide ) {
    Teuchos::RCP< Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > cg_nonlin_division_strategy =
       Teuchos::rcp(new Stokhos::CGDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 0, 1e-12, 4, 100, 1, 0,1));
    cg_nonlin_division_strategy->divide(setup.u, 1.0, setup.c1, setup.cu, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.cu, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }
  TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, GMRES_Divide ) {
    Teuchos::RCP< Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > gmres_division_strategy =
       Teuchos::rcp(new Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 0, 100, 0, 0,1));
    gmres_division_strategy->divide(setup.u, 1.0, setup.c1, setup.x, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.x, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }

  TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, GMRES_Jacobi_Divide ) {
    Teuchos::RCP< Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > gmres_division_strategy =
       Teuchos::rcp(new Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 2, 100, 0, 0,1));
    gmres_division_strategy->divide(setup.u, 1.0, setup.c1, setup.x, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.x, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }
 TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, GMRES_GaussSeidel_Divide ) {
    Teuchos::RCP< Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > gmres_division_strategy =
       Teuchos::rcp(new Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 3, 100, 0, 0,1));
    gmres_division_strategy->divide(setup.u, 1.0, setup.c1, setup.x, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.x, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }

TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, GMRES_Schur_Divide ) {
    Teuchos::RCP< Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > gmres_division_strategy =
       Teuchos::rcp(new Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 0, 1e-12, 4, 100, 0, 0,1));
    gmres_division_strategy->divide(setup.u, 1.0, setup.c1, setup.x, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.x, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }
TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, GMRES_Nonlin_Divide ) {
    Teuchos::RCP< Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > gmres_division_strategy =
       Teuchos::rcp(new Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 0, 100, 0, 0,1));
    gmres_division_strategy->divide(setup.u, 1.0, setup.c1, setup.cu, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.cu, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }
TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, GMRES_Nonlin_Jacobi_Divide ) {
    Teuchos::RCP< Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > gmres_division_strategy =
       Teuchos::rcp(new Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 2, 100, 0, 0,1));
    gmres_division_strategy->divide(setup.u, 1.0, setup.c1, setup.cu, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.cu, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);

  }
TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, GMRES_Nonlin_GaussSeidel_Divide ) {
    Teuchos::RCP< Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > gmres_division_strategy =
       Teuchos::rcp(new Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 1, 1e-12, 3, 100, 0, 0,1));
    gmres_division_strategy->divide(setup.u, 1.0, setup.c1, setup.cu, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.cu, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  
  }


TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, GMRES_Nonlin_Schur_Divide ) {
    Teuchos::RCP< Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > gmres_division_strategy =
       Teuchos::rcp(new Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 0, 1e-12, 4, 100, 0, 0,1));
    gmres_division_strategy->divide(setup.u, 1.0, setup.c1, setup.cu, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.cu, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }
TEUCHOS_UNIT_TEST( Stokhos_DivisionOperator, GMRES_Nonlin_Schur_linearprec_Divide ) {
    Teuchos::RCP< Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> > > gmres_division_strategy =
       Teuchos::rcp(new Stokhos::GMRESDivisionExpansionStrategy<int,double,Stokhos::StandardStorage<int, double> >(setup.basis, setup.Cijk, 0, 1e-12, 4, 100, 1, 0,1));
    gmres_division_strategy->divide(setup.u, 1.0, setup.c1, setup.cu, 0.0);
    setup.direct_division_strategy->divide(setup.u2, 1.0, setup.c1, setup.cu, 0.0);
    success = Stokhos::comparePCEs(setup.u, "u", setup.u2, "u2",
                                   setup.rtol, setup.atol, out);
  }


  


}
int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}


