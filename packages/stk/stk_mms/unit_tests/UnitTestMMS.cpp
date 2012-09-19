/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_mms/stk_mms.h>

namespace stk {
namespace mms {
namespace unit_tests {

//=============================================================================
//=============================================================================
//=============================================================================

STKUNIT_UNIT_TEST(derivatives, firstOrderPartials)
{
//   EXCEPTWATCH;

  const double xb=1.234, yb=2.345, zb=2.3, tb = 23.4;

  FAD_Type x(4, 0, xb), y(4, 1, yb), z(4, 2, zb), t(4, 3, tb);

  FAD_Type expr = (1.0+(x-t)*y*z)+(1.0+sin(x+y-t*z))+(2.0+exp(-x*x-t))+(3.0+pow(y-z*t,3));

  const double dx=DX(expr), dy=DY(expr), dz=DZ(expr), dt=DT(expr);

  // test against output from Mathematica
  STKUNIT_EXPECT_DOUBLE_EQ(expr.val(), -136504.5806378632);
  STKUNIT_EXPECT_DOUBLE_EQ(dx, 6.393200319571155);
  STKUNIT_EXPECT_DOUBLE_EQ(dy, 7899.044775319607);
  STKUNIT_EXPECT_DOUBLE_EQ(dz, -186082.60113247877);
  STKUNIT_EXPECT_DOUBLE_EQ(dt, -18290.45462323511);
}

STKUNIT_UNIT_TEST(exampleMMS, euler3D)
{
  const double xb=1.234, yb=2.345, zb=2.3, tb = 23.4;

  FAD_Type x(4, 0, xb), y(4, 1, yb), z(4, 2, zb), t(4, 3, tb);

  FAD_Type rho = 1.0+(x-t)*y*z;
  FAD_Type u = 1.0+sin(x+y-t*z);
  FAD_Type v = 2.0+exp(-x*x-t);
  FAD_Type w = 3.0+pow(y-z*t,3);
  FAD_Type T = 1.0+cos(z+x-t*y);

  const double Rgas = 287.0;
  FAD_Type p = Rgas*rho*T;

  const double gamma = 1.4;
  const double Cv = Rgas/(gamma-1.0);
  FAD_Type E = rho*(Cv*T + 0.5*(u*u+v*v+w*w));

  const double mass_src =   DT(rho)   + DX(rho*u)   + DY(rho*v)   + DZ(rho*w);
  const double momx_src =   DT(rho*u) + DX(rho*u*u) + DY(rho*u*v) + DZ(rho*u*w) + DX(p);
  const double momy_src =   DT(rho*v) + DX(rho*v*u) + DY(rho*v*v) + DZ(rho*v*w) + DY(p);
  const double momz_src =   DT(rho*w) + DX(rho*w*u) + DY(rho*w*v) + DZ(rho*w*w) + DZ(p);
  const double energy_src = DT(E)     + DX(u*(E+p)) + DY(v*(E+p)) + DZ(w*(E+p));

  // test against output from Mathematica
  STKUNIT_EXPECT_DOUBLE_EQ(mass_src, 2.914077175792222e7);
  STKUNIT_EXPECT_DOUBLE_EQ(momx_src, -3.4842036495958334e8);
  STKUNIT_EXPECT_DOUBLE_EQ(momy_src, 5.895967611187038e7);
  STKUNIT_EXPECT_DOUBLE_EQ(momz_src, -6.98207732332361e12);
  STKUNIT_EXPECT_DOUBLE_EQ(energy_src, 6.81241027548468e17);
}

STKUNIT_UNIT_TEST(exampleMMS, solidMechanics3D)
{
  // NOTE: some of the code below should be moved into library calls

  // point in ref config = unit cube
  const double Xb=0.8124, Yb=0.763, Zb=0.659, tb = 1.0;

  FAD2_Type X(4,0,FAD_Type(4,0,Xb));
  FAD2_Type Y(4,1,FAD_Type(4,0,Yb));
  FAD2_Type Z(4,2,FAD_Type(4,0,Zb));
  FAD2_Type t(4,3,FAD_Type(4,0,tb));

  // displacement function - user specified
  FAD2_Type phi[3];
  phi[0] = X + 0.2*sin(0.5*M_PI*t)*sin(M_PI*X)*sin(M_PI*Y)*sin(M_PI*Z);
  phi[1] = Y + 0.3*sin(0.5*M_PI*t)*sin(M_PI*X)*sin(M_PI*Y)*sin(M_PI*Z);
  phi[2] = Z + 0.4*sin(0.5*M_PI*t)*sin(M_PI*X)*sin(M_PI*Y)*sin(M_PI*Z);

  // deformation gradient - library call
  FAD_Type F[3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      F[i][j] = phi[i].dx(j);
    }
  }

  // determinant - library call
  FAD_Type J = 
     F[0][0]*(F[1][1]*F[2][2]-F[1][2]*F[2][1])
    -F[0][1]*(F[1][0]*F[2][2]-F[1][2]*F[2][0])
    +F[0][2]*(F[1][0]*F[2][1]-F[1][1]*F[2][0]);

  // left Cauchy-Green tensor (b = F*F^t) - library call
  FAD_Type b[3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      b[i][j] = 0;
      for (int k=0; k<3; k++) {
	b[i][j] += F[i][k]*F[j][k];
      }
    }
  }      

  const double mu  = 1.0;
  const double lambda = 1.0;
  const double bulk = lambda + mu*(2.0/3.0);

  // Cauchy stress tensor - material dependent.  
  // Here we use neo-Hookean
  FAD_Type sigma[3][3];
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      sigma[i][j] = 0;
    }
  }

  FAD_Type pressure = 0.5*bulk*(J-1.0/J);

  FAD_Type factor = pow(J, -5.0/3.0);
  FAD_Type trace_b = b[0][0]+b[1][1]+b[2][2];

  for (int i=0; i<3; i++) {
    sigma[i][i] += pressure - mu*factor*trace_b/3.0;
    for (int j=0; j<3; j++) {
      sigma[i][j] += mu*factor*b[i][j];
    }
  }

  // inverse deformation gradient - library call
  FAD_Type inv_F[3][3];
  
  inv_F[0][0] =  (F[2][2]*F[1][1]-F[2][1]*F[1][2])/J;
  inv_F[1][0] = -(F[2][2]*F[1][0]-F[2][0]*F[1][2])/J;
  inv_F[2][0] =  (F[2][1]*F[1][0]-F[2][0]*F[1][1])/J;

  inv_F[0][1] = -(F[2][2]*F[0][1]-F[2][1]*F[0][2])/J;
  inv_F[1][1] =  (F[2][2]*F[0][0]-F[2][0]*F[0][2])/J;
  inv_F[2][1] = -(F[2][1]*F[0][0]-F[2][0]*F[0][1])/J;
  
  inv_F[0][2] =  (F[1][2]*F[0][1]-F[1][1]*F[0][2])/J;
  inv_F[1][2] = -(F[1][2]*F[0][0]-F[1][0]*F[0][2])/J;
  inv_F[2][2] =  (F[1][1]*F[0][0]-F[1][0]*F[0][1])/J;
 
  // source term = -div(sigma) + rho*ddot(u) ...
  double mom_src[3];

  const double rho = 1.0;

  // compute div(sigma) in current coordinates using chain rule

  // this is because we compute everything based on a point 
  // in the original configuration, but need the divergence 
  // w.r.t. the current configuration

  for (int i=0; i<3; i++) {
    mom_src[i] = 0;
    for (int j=0; j<3; j++) {
      for (int k=0; k<3; k++) {
	// mom_i = - D_x_j sigma_ij = - D_X_k sigma_ij D_X_k/D_x_j
	mom_src[i] -= (sigma[i][j]).dx(k)*inv_F[k][j].val();
      }
    }
    mom_src[i] += rho*phi[i].dx(_T).dx(_T);
  }
  
  // TODO: test against output from Mathematica 

  // DEBUG
  std::cout << "Time: " << t.val().val() << std::endl;
  std::cout << "Initial point: " 
	    << X.val().val() << " " 
	    << Y.val().val() << " " 
	    << Z.val().val() << std::endl;
  std::cout << "Deformed point: " 
	    << phi[0].val().val() << " " 
	    << phi[1].val().val() << " " 
	    << phi[2].val().val() << std::endl;
  std::cout << "Deformation gradient: " << std::endl;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      std::cout << F[i][j].val() << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Inverse deformation gradient: " << std::endl;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      std::cout << inv_F[i][j].val() << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Determinant: " << J.val() << std::endl;
  std::cout << "Left Cauchy-Green tensor: " << std::endl;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      std::cout << b[i][j].val() << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "Pressure: " << pressure.val() << std::endl;
  std::cout << "Cauchy stress tensor: " << std::endl;
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      std::cout << sigma[i][j].val() << " ";
    }
    std::cout << std::endl;
  }
  std::cout << "momentum source: " << std::endl;
  for (int i=0; i<3; i++) {
    std::cout << mom_src[i] << " ";
  }
  std::cout << std::endl;
}

STKUNIT_UNIT_TEST(exampleMMS, viscousNavierStokes2D)
{
  const int ND = 2;

  const double gamma = 1.4;
  const double Rgas = 287.097384766765;

  const double suth_c1 = 0.01458;
  const double suth_c2 = 110.4;
  const double Pr = 0.72;

  const double Cv = Rgas/(gamma-1.0);
  
  const double M0 = 2.5;
  const double p0 = 35651.28116;
  const double T0 = 236.215;
  
  const double rho0 = p0/(Rgas*T0);
  const double u0 = M0*sqrt(gamma*p0/rho0);

  const double eps = 0.05;
  const double Cp = gamma*Cv;

  const double xb=1.34, yb=0.356, tb = 1.0;
  FAD2_Type x(4, 0, xb), y(4, 1, yb), t(4, 3, tb);

  x.val() = FAD_Type(4, 0, xb);
  y.val() = FAD_Type(4, 1, yb);
  t.val() = FAD_Type(4, 3, tb);

  FAD2_Type u = u0*(1.0+eps*pow(sin(0.5*M_PI*x),2)*pow(sin(M_PI*y),2));
  FAD2_Type v = u0*eps*pow(sin(0.5*M_PI*x),2)*pow(sin(M_PI*y),2);
  FAD2_Type p = p0*(1.0 + eps*sin(0.5*M_PI*x)*cos(M_PI*y));
  FAD2_Type T = T0*(1.0 + eps*sin(0.25*M_PI*x)*cos(M_PI*y));

  FAD2_Type velocity[3];
  velocity[0] = u;
  velocity[1] = v;
  velocity[2] = 0.0;

  FAD_Type viscosity = ( suth_c1 * T * sqrt(T)/(T + suth_c2) ).val();

  FAD_Type thermal_conductivity = viscosity*Cp/Pr;
  
  // form stress tensor (no pressure)
  FAD_Type tau[3][3];
  const double oneThird = 1./3.;
  
  FAD_Type div_vel = 0.0;
  for (int i = 0; i < ND; ++i ) {
    div_vel += D2(velocity[i], i);
  }

  for (int i = 0; i < ND; ++i ) {
    tau[i][i] = 2.0*viscosity*(D2(velocity[i], i) - oneThird*div_vel);
    
    for (int j = 0; j < i; ++j ) {
      tau[i][j] = viscosity*(D2(velocity[i], j) + D2(velocity[j], i));
      tau[j][i] = tau[i][j];
    }
  }
  
  // form total thermal heat flux
  FAD_Type thermal_flux[3];
  
  for (int i = 0; i < ND; ++i ) {
    thermal_flux[i] = -thermal_conductivity*D2(T, i);
  }
  
  for (int i = 0; i < ND; ++i ) {
    for (int j = 0; j < ND; ++j ) {
      thermal_flux[i] -= velocity[j].val()*tau[j][i];
    }
  }
  
  double mom_src[3];
  for (int i = 0; i < ND; ++i ) {
    mom_src[i] = 0.0;
  }

  double energy_src = 0.0;

  // sum terms into src (using divg)
  for (int i = 0; i < ND; ++i ) {
    energy_src += D(thermal_flux[i], i);
    for (int j = 0; j < ND; ++j ) {
      mom_src[i] -= D(tau[j][i], j);
    }
  }
  
  // test against output from Mathematica
  STKUNIT_EXPECT_NEAR(mom_src[0], 76.58413133799213, 1e-10);
  STKUNIT_EXPECT_NEAR(mom_src[1], 92.05081838663132, 1e-10);
  STKUNIT_EXPECT_NEAR(energy_src, 70295.4770236427,  1e-10);
}

} // namespace unit_tests
} // namespace mms
} // namespace stk

