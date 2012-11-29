/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_percept/Percept.hpp>
#include <stk_percept/Util.hpp>
#include <stk_percept/math/Math.hpp>
#include <stk_percept/math/DenseMatrix.hpp>
#include <stk_percept/fixtures/Fixture.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/diag/PrintTable.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <Teuchos_ScalarTraits.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>

namespace stk {
namespace percept {
namespace unit_tests {

//=============================================================================
//=============================================================================
//=============================================================================

static void test_eigen(double A[3][3], double eig_expected[3])
{
  double eig[3];
  DenseMatrix<3,3> M(A);
  eigen_3x3(M, eig);

  const double tol = 1.e-4;

  STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eig[0], eig_expected[0], tol);
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eig[1], eig_expected[1], tol);
  STKUNIT_EXPECT_DOUBLE_EQ_APPROX_TOL(eig[2], eig_expected[2], tol);
}

/// Code generated from Mathematica notebook @see eigen3x3.1.nb 
STKUNIT_UNIT_TEST(unit_math, eigen)
{
  EXCEPTWATCH;
  MPI_Barrier( MPI_COMM_WORLD );

double m0[3][3] = {{0.755447, 0.111342, 0.37761}, {0.111342, 0.651947, 0.267815}, {0.37761, 0.267815, 0.269684}}; double ei0[3] = {1.0880529453139673, 0.5813510432810213, 0.007673854192603399};

double m1[3][3] = {{0.16826, 0.789665, 0.00973783}, {0.789665, 0.711058, 0.382209}, {0.00973783, 0.382209, 0.949962}}; double ei1[3] = {1.4750846518274794, 0.7842389812485816, -0.4300438040795849};

double m2[3][3] = {{0.773226, 0.452952, 0.27961}, {0.452952, 0.865696, 0.022551}, {0.27961, 0.022551, 0.646788}}; double ei2[3] = {1.3385147226626053, 0.6894027149103933, 0.257792593153505};

double m3[3][3] = {{0.909289, 0.501889, 0.45783}, {0.501889, 0.238991, 0.900382}, {0.45783, 0.900382, 0.469206}}; double ei3[3] = {1.784328235733481, 0.39129509787679345, -0.5581366853170395};

double m4[3][3] = {{0.702383, 0.127649, 0.522772}, {0.127649, 0.81726, 0.434568}, {0.522772, 0.434568, 0.857965}}; double ei4[3] = {1.5498268486692401, 0.6480407900710551, 0.17973967795170437};

double m5[3][3] = {{0.354512, 0.0275941, 0.42483}, {0.0275941, 0.146907, 0.972302}, {0.42483, 0.972302, 0.077632}}; double ei5[3] = {1.2071403752044616, 0.29917540311243573, -0.927265061302539};

double m6[3][3] = {{0.651604, 0.693955, 0.692692}, {0.693955, 0.211936, 0.629053}, {0.692692, 0.629053, 0.0471666}}; double ei6[3] = {1.6887408967270607, -0.26366985417330047, -0.514363859422579};

double m7[3][3] = {{0.783403, 0.710047, 0.171223}, {0.710047, 0.808176, 0.883021}, {0.171223, 0.883021, 0.240841}}; double ei7[3] = {1.8802888922785086, 0.43140165441588174, -0.4792708055096206};

double m8[3][3] = {{0.46884, 0.680527, 0.360249}, {0.680527, 0.423581, 0.034272}, {0.360249, 0.034272, 0.822562}}; double ei8[3] = {1.3050430327295293, 0.6932431548874134, -0.2833022831435616};

double m9[3][3] = {{0.00573793, 0.395987, 0.609442}, {0.395987, 0.675655, 0.0334358}, {0.609442, 0.0334358, 0.318355}}; double ei9[3] = {1.032679556235905, 0.5077523809442646, -0.5406837254017315};

double m10[3][3] = {{0.957837, 0, 0}, {0, 0.9817, 0}, {0, 0, 0.340743}}; double ei10[3] = {0.9817003436069414, 0.9578372391289818, 0.3407434051765223};

double m11[3][3] = {{0.106419, 0, 0}, {0, 0.328784, 0}, {0, 0, 0.934534}}; double ei11[3] = {0.9345337227256809, 0.3287839495964591, 0.10641896032545504};

double m12[3][3] = {{0.55734, 0, 0}, {0, 0.396372, 0}, {0, 0, 0.157561}}; double ei12[3] = {0.5573403962034714, 0.3963716192114425, 0.15756074821379618};

double m13[3][3] = {{0.126358, 0, 0}, {0, 0.674319, 0}, {0, 0, 0.155531}}; double ei13[3] = {0.6743193624291033, 0.15553063699815822, 0.12635797272724567};

double m14[3][3] = {{0.68872, 0, 0}, {0, 0.445831, 0}, {0, 0, 0.31407}}; double ei14[3] = {0.6887204656460029, 0.44583093515781363, 0.3140699233021876};

double m15[3][3] = {{0.731949, 0, 0}, {0, 0.654448, 0}, {0, 0, 0.623269}}; double ei15[3] = {0.7319492335383265, 0.6544484721945408, 0.6232687167120586};

double m16[3][3] = {{0.308332, 0, 0}, {0, 0.335962, 0}, {0, 0, 0.0450069}}; double ei16[3] = {0.33596196231820047, 0.3083319964015657, 0.04500694087060534};

double m17[3][3] = {{0.947614, 0, 0}, {0, 0.274896, 0}, {0, 0, 0.0176067}}; double ei17[3] = {0.9476136622146641, 0.27489616374231096, 0.017606731937778692};

double m18[3][3] = {{0.0871697, 0, 0}, {0, 0.965913, 0}, {0, 0, 0.934153}}; double ei18[3] = {0.9659133186077228, 0.9341527585657887, 0.08716970174162356};

double m19[3][3] = {{0.911188, 0, 0}, {0, 0.758386, 0}, {0, 0, 0.0313796}}; double ei19[3] = {0.9111877716123237, 0.7583857521451645, 0.03137959588204193};

test_eigen( m0, ei0);

test_eigen( m1, ei1);

test_eigen( m2, ei2);

test_eigen( m3, ei3);

test_eigen( m4, ei4);

test_eigen( m5, ei5);

test_eigen( m6, ei6);

test_eigen( m7, ei7);

test_eigen( m8, ei8);

test_eigen( m9, ei9);

test_eigen( m10, ei10);

test_eigen( m11, ei11);

test_eigen( m12, ei12);

test_eigen( m13, ei13);

test_eigen( m14, ei14);

test_eigen( m15, ei15);

test_eigen( m16, ei16);

test_eigen( m17, ei17);

test_eigen( m18, ei18);

test_eigen( m19, ei19);

  
}

}
}
}
