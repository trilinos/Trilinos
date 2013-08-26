/*--------------------------------------------------------------------*/
/*    Copyright 2009, 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_percept/mesh/geometry/stk_geom/SplineFit.hpp>
#include <stk_percept/mesh/geometry/stk_geom/LocalCubicSplineFit.hpp>
#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>
#include <math.h>

namespace stk
{
  namespace geom
  {
    namespace unit_tests
    {

      //=============================================================================
      //=============================================================================
      //=============================================================================

      STKUNIT_UNIT_TEST(unit_stk_geom, test_1)
      {
        const bool debug_print = true;

        /// fit a quadratic
        LocalCubicSplineFit cf;
        int n = 3;
        Vectors2D Q(n);
        for (int i=0; i < n; i++)
          {
            double x = double(i)/double(n-1);
            //double y = x*x;
            double y = x;
            if (i==1) {
              x=.2;
              y=.2;
            }
            double xy[] = {x, y};
            Q[i] = Vector2D(xy);
          }
        if (debug_print) std::cout << "Fitting quadratic with cubic spline - points = " << std::endl;
        DPRINTLN(Q);
        ON_Curve *curve = cf.fit(Q);
        if (debug_print) 
          cf.print();

        double t0=0,t1=1;
        //ON_BOOL32 nerr = curve->GetDomain( &t0, &t1);
        curve->GetDomain( &t0, &t1);
        DPRINTLN2(t0,t1);
        DPRINTLN(curve->SpanCount());
        std::vector<double> spanKnots(curve->SpanCount()+1);
        curve->GetSpanVector(&spanKnots[0]);
        DPRINTLN(spanKnots);
        DPRINTLN(curve->HasNurbForm());
        DPRINTLN(curve->PointAtStart());
        DPRINTLN(curve->PointAtEnd());
        for (size_t i=0; i < cf.m_U.size(); i++)
          {
            DPRINTLN2(i, curve->PointAt(cf.m_U[i]));
          }
        // this tolerance is too high - why is this evaluation not more accurate?  FIXME
        double tol = 1.e-4;
        int kk=0;
        for (size_t i=1; i < cf.m_U.size()-1; i += 2)
          {
            Point3D pt = curve->PointAt(cf.m_U[i]);
            Point3D Qkk = Q[kk++];
            DPRINTLN2(pt, Qkk);
            double dist = pt.DistanceTo(Qkk);
            STKUNIT_EXPECT_NEAR(dist, 0, tol);
          }
        for (size_t i=0; i < cf.m_U.size(); i++)
          {
            DPRINTLN2(i, curve->DerivativeAt(cf.m_U[i]));
          }
        for (size_t i=0; i < cf.m_U.size(); i++)
          {
            DPRINTLN2(i, curve->CurvatureAt(cf.m_U[i]));
          }

        for (size_t i=0; i < cf.m_U.size(); i++)
          {
            ON_3dPoint point;
            ON_3dVector first_derivative;
            ON_3dVector second_derivative;
            curve->Ev2Der(cf.m_U[i], point, first_derivative, second_derivative);
            if (debug_print)
              std::cout << PR(i) << PR(cf.m_U[i]) << PR(point)
                        << PR(first_derivative) << PR(second_derivative) << std::endl;
          }

        {
          Point3D p(0,.01,0);
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          STKUNIT_EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN(closest_point);
          double dist = closest_point.DistanceTo(p);
          STKUNIT_EXPECT_NEAR(dist, 0.01, tol);
          STKUNIT_EXPECT_NEAR(closest_point[0], 0.0, tol);
          STKUNIT_EXPECT_NEAR(closest_point[1], 0.0, tol);
        }

        {
          Point3D p(0.5,0.25,0);
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          STKUNIT_EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN2(u,closest_point);
          double dist = closest_point.DistanceTo(p);
          STKUNIT_EXPECT_NEAR(dist, 0,tol);
        }

        {
          Point3D p(1,1,0);
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          STKUNIT_EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN2(u,closest_point);
          double dist = closest_point.DistanceTo(p);
          STKUNIT_EXPECT_NEAR(dist, 0,tol);
        }

        {
          int nn = 101;
          for (int i = 0; i < nn; i++)
            {
              double t = double(i)/double(nn-1);
              Point3D p = Q[0] + t*(Q[n-1]-Q[0]);
              double u=0;
              bool success = curve->GetClosestPoint(p, &u);
              STKUNIT_EXPECT_TRUE(success);
              Point3D closest_point = curve->PointAt(u);
              DPRINTLN2(u,closest_point);
              double dist = closest_point.DistanceTo(p);
              STKUNIT_EXPECT_NEAR(dist, 0,tol);
            }
        }

        {
          Point3D p(1e-6,1e-6,0);
          double u=0;
          bool success = curve->GetClosestPoint(p, &u);
          STKUNIT_EXPECT_TRUE(success);
          Point3D closest_point = curve->PointAt(u);
          DPRINTLN2(u,closest_point);
          double dist = closest_point.DistanceTo(p);
          STKUNIT_EXPECT_NEAR(dist, 0,tol);
        }


        delete curve;
      }

      STKUNIT_UNIT_TEST(unit_stk_geom, test_2)
      {
        const bool debug_print = true;
        if (1)
          {
            // write a wiggly cubic curve on the "green NURBS wiggle" layer
            ON_NurbsCurve* wiggle = new ON_NurbsCurve(
                                                      3, // dimension
                                                      false, // true if rational
                                                      4,     // order = degree+1
                                                      6      // number of control vertices
                                                      );
            int i;
            Points3D pts;
            for ( i = 0; i < wiggle->CVCount(); i++ ) {
              Point3D pt( 2*i, -i, (i-3)*(i-3) ); // pt = some 3d point
              DPRINTLN2(i,pt);
              wiggle->SetCV( i, pt );
              pts.push_back(pt);
            }

            // ON_NurbsCurve's have order+cv_count-2 knots.
            wiggle->SetKnot(0, 0.0);
            wiggle->SetKnot(1, 0.0);
            wiggle->SetKnot(2, 0.0);
            wiggle->SetKnot(3, 1.5);
            wiggle->SetKnot(4, 2.3);
            wiggle->SetKnot(5, 4.0);
            wiggle->SetKnot(6, 4.0);
            wiggle->SetKnot(7, 4.0);

            DPRINTLN(wiggle->IsValid());
            DPRINTLN(wiggle->IsRational());
            DPRINTLN(wiggle->CVCount());
            DPRINTLN(wiggle->CVSize());
            DPRINTLN(wiggle->Order());
            DPRINTLN(wiggle->KnotCount());

            DPRINTLN(wiggle->HasNurbForm());
            DPRINTLN(wiggle->PointAtStart());
            DPRINTLN(wiggle->PointAtEnd());

            int kn=wiggle->KnotCount();
            std::vector<double> knots;
            for (int k = 0; k < kn; k++)
              {
                DPRINTLN2(k, wiggle->PointAt(wiggle->Knot(k)));
                knots.push_back(wiggle->Knot(k));
              }
            ON_TextLog log;
            wiggle->Dump(log);

            DPRINTLN(pts);
            DPRINTLN(knots);

            delete wiggle;
          }

        if (0)
          {
            double pts[][2] = {{0, 1},     {0.2, 1.3}, {1.5, 2},
                               {1.75, 2},  {2, 2}, {2.1, 0.5},
                               {2.5, .5},  {2.9, .5}, {3.3, .9},
                               {4, 2}};
            //  Piegl & Tiller convention - start/end with order number of repeated knots
            //double[] knots =  {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7};
            // OpenNURBS convention - start/end with (order-1) repeated knots
            //double knots[] =  {0, 0,   0, 1, 2, 3, 4, 5, 6, 7,   7, 7 };
            double knots[] =  {0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4 };

            // n = 3, Q[0]..Q[n] --> n = 3;
            int np=10;
            //int nk=12;

            ON_NurbsCurve* wiggle = new ON_NurbsCurve(
                                                      2, // dimension
                                                      false, // true if rational
                                                      4,     // order = degree+1
                                                      np      // number of control vertices
                                                      );
            int i;
            for ( i = 0; i < wiggle->CVCount(); i++ ) {
              Point3D pt( pts[i][0], pts[i][1], 0);
              DPRINTLN2(i,pt);
              wiggle->SetCV( i, pt );
            }

            // ON_NurbsCurve's have order+cv_count-2 knots.
            int kn=wiggle->KnotCount();
            for (int k = 0; k < kn; k++)
              {
                wiggle->SetKnot(k, knots[k]);
              }

            DPRINTLN(wiggle->IsValid());
            DPRINTLN(wiggle->IsRational());
            DPRINTLN(wiggle->CVCount());
            DPRINTLN(wiggle->CVSize());
            DPRINTLN(wiggle->Order());
            DPRINTLN(wiggle->KnotCount());

            DPRINTLN(wiggle->HasNurbForm());
            DPRINTLN(wiggle->PointAtStart());
            DPRINTLN(wiggle->PointAtEnd());

            for (int k = 0; k < kn; k++)
              {
                DPRINTLN2(k, wiggle->PointAt(wiggle->Knot(k)));
              }
            delete wiggle;
          }

      }

    }
  }
}
