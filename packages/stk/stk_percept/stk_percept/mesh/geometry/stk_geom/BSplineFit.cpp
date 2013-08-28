#include <stk_percept/mesh/geometry/stk_geom/BSplineFit.hpp>

namespace stk {
  namespace geom {


    enum Option { ThreePoint = 3, FivePoint = 5 };

    ON_Curve * BSplineFit::
    fit(Vectors2D& points_in)
    {
      const bool debug_print = SplineFit::s_debug_print;
      int n_i = points_in.size();
      Option option = ThreePoint; // FivePoint
      //int n_stencil = option; // option.n_stencil;
      int n = n_i - 1; // 0...n define the points

      /// using Piegl and Tiller "The NURBS Book" notation and equation numbers
      //Vectors2D& Q = *(&tmp_points[n_stencil]);
      Vectors2D Q = points_in;  // FIXME
      m_Q = Q;

      double arclen = 0.0;
      for (int k=1; k <= n; k++)
        {
          // eq (9.3)
          //u[k] = double(k)/double(n+1);
          // eq (9.4,5)
          arclen += (Q[k] - Q[k-1]).Length();
        }

      // decide if points are closed
      m_isClosed = false;
      double reltol = 1.e-8;
      double tol = reltol*arclen;
      Vector2D diff = Q[0];
      diff -= Q[n];
      if (diff.Length() < tol)
        {
          m_isClosed = true;
        }

      Vectors2D q(n+1), d(n+1), D(n+1), V(n+1), T(n+1);
      std::vector<double> u(n+1), alpha(n+1);
      u[0] = 0.0;
      u[n] = 1.0;

      for (int k=1; k <= n-1; k++)
        {
          // eq (9.4,5)
          u[k] = u[k-1] + (Q[k] - Q[k-1]).Length() / arclen;
          DPRINTLN2(k,u[k]);
        }

      for (int k=1; k <= n; k++)
        {
          // eq (9.28)
          q[k] = Q[k] - Q[k-1];
          d[k] = q[k]/(u[k] - u[k-1]);
        }
      if (option == ThreePoint)
        {
          int N = n-1;
          if (m_isPeriodic)
            {
              double unp1 = u[n] + (u[1]-u[0]);
              u.push_back(unp1);
              q.push_back(q[1]);
              d.push_back(d[1]);
              N = n;
            }
          DPRINTLN(u);
          for (int k=1; k <= N; k++)
            {
              // eq (9.30)
              alpha[k] = (u[k] - u[k-1]) / (u[k+1] - u[k-1]);
              double akm1 = (u[k+1] - u[k]) / (u[k+1] - u[k-1]);
              // eq (9.28)
              D[k] = akm1*d[k] + alpha[k]*d[k+1];
              V[k] = akm1*q[k] + alpha[k]*q[k+1];
              double vkl =  V[k].Length();
              if (vkl < 1.e-12) throw std::runtime_error("can't normalize T vector");
              DPRINTLN2(k,alpha[k]);
              // eq (9.29)
              T[k] = V[k] / vkl;
              DPRINTLN2(k,T[k]);
            }
          if (m_isPeriodic)
            {
              T[0] = T[n];
            }
          else
            {
              // endpoints
              // eq (9.32)
              D[0] = d[1] + (d[1] - D[1]);
              D[n] = d[n] + (d[n] - D[n-1]);
              V[0] = 2*q[1] - V[1];
              V[n] = 2*q[n] - V[n-1];
              double v0l = V[0].Length();
              double vnl = V[n].Length();
              if (v0l < 1.e-12) throw std::runtime_error("can't normalize T vector[0]");
              if (vnl < 1.e-12) throw std::runtime_error("can't normalize T vector[n]");
              T[0] = V[0]/v0l;
              T[n] = V[n]/vnl;
              DPRINTLN(d[1]);
              DPRINTLN(D[1]);
              DPRINTLN(T[0]);
              DPRINTLN(T[n]);
            }
        }
      else if (option == FivePoint)
        {
          throw std::runtime_error("BSplineFit FivePoint option not yet implemented");
        }
      if (0 && m_isPeriodic && m_isClosed)
        {
          Vector2D Tave = 0.5*(T[0] + T[n]);
          T[0] = Tave;
          T[n] = Tave;
        }
      fit_internal(n, Q, T);
      return create();
    }

    ON_Curve * BSplineFit::create()
    {
      const bool debug_print = SplineFit::s_debug_print;
      const int dimension = 2;
      const int order = 4;

      ON_NurbsCurve* pNurbsCurve =
        ON_NurbsCurve::New(
                           dimension,         // int dimension,
                           false,             // ON_BOOL32 bIsRational,
                           order,             // int order = degree+1
                           (int)m_CV.size()   // int cv_count
                           );
      if (!pNurbsCurve) throw std::runtime_error("failed in ON_NurbsCurve");

      for (int k = 0; k < (int)m_CV.size(); k++)
        {
          Point3D v;
          v.Zero();
          v[0] = m_CV[k][0];
          v[1] = m_CV[k][1];
          ON_BOOL32 success = pNurbsCurve->SetCV(k, v);
          if (!success) {
            PRINTLN(k);
            throw std::runtime_error("failed in SetCV");
          }
        }
      DPRINTLN(m_CV.size());
      DPRINTLN(m_U.size());
      DPRINTLN(pNurbsCurve->KnotCount());

      for (int k = 0; k < (int)m_U.size(); k++)
        {
          ON_BOOL32 success = pNurbsCurve->SetKnot(k, m_U[k]);
          if (!success) {
            PRINTLN(k);
            throw std::runtime_error("failed in SetKnot");
          }
        }
      m_curve = pNurbsCurve;
      return pNurbsCurve;
    }


    void BSplineFit::print()
    {
      std::cout << "BSplineFit::print: getIsPeriodic= " << getIsPeriodic() << " getIsClosed= " << getIsClosed() << std::endl;
      std::cout << "n = " << m_n << " cv count= " << m_CV.size() << " U count = " << m_U.size() << std::endl;
      std::cout << "Q= " << m_Q << std::endl;
      std::cout << "CV= " << m_CV << std::endl;
      std::cout << "U= " << m_U << std::endl;
      if (m_curve)
        {
          PRINTLN(m_curve->IsValid());
          PRINTLN(m_curve->IsRational());
          PRINTLN(m_curve->CVCount());
          PRINTLN(m_curve->CVSize());
          PRINTLN(m_curve->Order());
          PRINTLN(m_curve->KnotCount());
          ON_TextLog log;
          m_curve->Dump(log);
        }
    }

  }
}

