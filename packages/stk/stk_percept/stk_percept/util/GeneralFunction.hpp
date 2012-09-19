/*--------------------------------------------------------------------*/
/*    Copyright 2010 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef stk_percept_util_generalfunction_hpp
#define stk_percept_util_generalfunction_hpp

#include <vector>

using namespace std;

#if 0
template<typename T> void push_back( vector<T>& dst, const vector<T>& src)
{
  dst.insert(dst.end(), src.begin(), src.end());
}
#endif

namespace stk {

  namespace percept {

    namespace util {

      /** example usage
       *
       *   GeneralFunction<double, vector<double> > f;
       *   GeneralFunction<vector<double>, int > g;
       *   vector<double> y;
       *   double x;
       *   y = f(x);
       *   f(x, y);
       *
       *   // multiple points
       *   vector<double> xv;
       *   vector< vector<double> > yv;
       *   f(xv, yv);
       *
       *   // function composition
       *   GeneralFunction<double, int> g_of_f = g(f);
       *   int iy;
       *   iy = g_of_f(x);
       *
       *   // gradient
       *   GeneralFunction<double, vector<vector<double> > > grad_f = f.grad();
       *
       */
  
      template<typename domain_f, typename codomain_f_and_domain_g, typename codomain_g> 
      class CompositeGeneralFunction;

      template<typename domain, typename codomain> 
      class GeneralFunction
      {
      public:
        GeneralFunction() {}
  
        /// single value (default impl is identity op)
        virtual codomain operator()(const domain& x) { return x; }  // return value or reference?  
        virtual void operator()(const domain& x, codomain& y) { y = x; }

        // multiple values
        virtual vector<codomain> operator()(const vector<domain>& x) 
        {
          // inefficient default impl
          int n = x.size();
          vector<codomain> y(n);
          for(int i = 0; i < n; i++)
            {
              y[i] = (*this)(x[i]);
            }
          return y;
        };  // return value or reference?
        virtual void operator()(const vector<domain>& x, vector<codomain>& y) 
        {
          // inefficient default impl
          int n = x.size();
          for(int i = 0; i < n; i++)
            {
              y[i] = (*this)(x[i]);
            }
        }

        // composition
        template< typename codomain_f_and_domain_g>
        GeneralFunction<domain, codomain> operator()(const GeneralFunction<domain, codomain_f_and_domain_g>& f) 
        {
          return CompositeGeneralFunction<domain, codomain_f_and_domain_g, codomain>(*this, f);
        }

      };

      template<typename domain, typename codomain> 
      class GeneralFunctionWithGrad : public GeneralFunction<domain, codomain>
      {
        // return a function that computes the gradient of this
        virtual GeneralFunction<domain, vector<codomain> > grad()=0;  // return GeneralFunction<domain, vector<codomain> >(); }
      };

      template<typename domain_f, typename codomain_f_and_domain_g, typename codomain_g> 
      class CompositeGeneralFunction : public GeneralFunction<domain_f, codomain_g>
      {
      private:
        GeneralFunction<codomain_f_and_domain_g, codomain_g> m_g;
        GeneralFunction<domain_f, codomain_f_and_domain_g> m_f;
        codomain_f_and_domain_g m_tmp_fv;

      public:
        CompositeGeneralFunction(GeneralFunction<codomain_f_and_domain_g, codomain_g>& g, GeneralFunction<domain_f, codomain_f_and_domain_g>& f, 
                          codomain_f_and_domain_g& tmp_fv = 0) : m_g(g), m_f(f), m_tmp_fv(tmp_fv? tmp_fv : codomain_f_and_domain_g())
        {
        }
        virtual codomain_g operator()(domain_f& x) { return m_g(m_f(x)); }
        virtual void operator()(domain_f& x, codomain_g& y) { 
          codomain_f_and_domain_g& fv =  m_tmp_fv;
          m_f(x, fv); 
          m_g(fv, y); 
        }
      };

    } // namespace util
  }//namespace percept
}//namespace stk

#endif
