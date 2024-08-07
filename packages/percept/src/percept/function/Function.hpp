// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_Function_hpp
#define stk_encr_Function_hpp


#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>

#include <percept/function/internal/Dimensions.hpp>
#include <percept/function/internal/GenericFunction.hpp>
#include <percept/Util.hpp>
#include <stk_mesh/base/Entity.hpp>

#include <Teuchos_RCP.hpp>

  namespace percept
  {

    class Function : public GenericFunction // , HasGradient
    {
    public:
      typedef std::map<std::string , Function * > NameToFunctionMap;

      /**
       * Create a function with the given name, domain dimensions and codomain dimensions, and integration order.
       *  If domain_dimensions and codomain_dimensions are defaulted, then they are dimensioned using the
       *  defaults of domain_dimensions = {s_spatialDimDefault} and 
       *              codomain_dimensions = {s_codomainDimDefault}
       *  If integration_order is not specified, it defaults to s_integration_order_default.
       *  The defaults of the defaults is s_spatialDimDefault=3, s_codomainDimDefault=1, and s_integration_order_default=1
       *
       *  [DEPRECATED START] NOTE: we assume that input arrays have either {x,y,t} 2D+time or {x,y,z,t} for 3D+time.  There is no separate
       *  argument for time.  This means arrays should be dimensioned with length 3 (2D+time), or 4 (3D+time) respectively.
       *  [DEPRECATED END]: now we have time in the operator() arg list
       *
       *  A Function has a {\em core} dimensioning of its domain and codamain, as specified in its constructor (or by
       *      setting with accessors setDomainDimensions() and setCodomainDimensions().  This means that all Function's
       *      implementations of  operator()(MDArray& input, MDArray& output) expect the rightmost dimensions of input
       *      to be equal to the rightmost dimensions of its domain dimensions, and similar for output and codomain_dimensions
       *   
       *  Conventions:
       *  1. core dimensions given by getDomainDimensions() and getCodomainDimensions() properties
       *  2. rightmost dimensions of input array must match rightmost dimensions of Function's domain
       *  3. rightmost dimensions of output array must match rightmost dimensions of Function's codomain
       *  4. domain and codomain can have different ranks and dimensions
       *  5. usage of conventions 1-4 are in scenarios where the input and output arrays contain multiple points
       *       that the client is requesting the Function to evaluate.  For example:
       *
       *         int numElements   = bucket.size();
       *         int numQuadPoints = quadratureFactory.getNumPoints();
       *         int spaceDim      = 3;
       *         int numValuesAtQP = 1;
       *
       *         MDArray quadrature_points(numElements, numQuadPoints, spaceDim);
       *         MDArray function_value_at_quadrature_points(numElements, numQuadPoints, numValuesAtQP);
       *
       *         StringFunction sf("<... some definition ...>", "<... some name ...>", Dimensions(spaceDim), Dimensions(numValuesAtQP) );
       *
       *         MDArray single_point(spaceDim);
       *         MDArray value_at_single_point(numValuesAtQP);
       *
       *         // note that this same sf can be evaluated on either of these pairs of in/out arrays
       *         string_function(quadrature_points, function_value_at_quadrature_points);
       *         string_function(single_point, value_at_single_point);
       *
       */
      Function() {throw new std::runtime_error("shouldn't be invoked"); }
      Function(const char *name, 
               Dimensions domain_dimensions = Dimensions(),
               Dimensions codomain_dimensions = Dimensions(),
               unsigned integration_order = 0);

      /// this version uses the MDOutVal argument as a predefined array to ensure the python returned value is properly sized
      void value(MDArray& domain, MDArray& MDOutVal, double time_value_optional=0.0)
      {
        this->operator()(domain, MDOutVal, time_value_optional);
      }

      /// this version creates a dummy output and returns it based on the codomain-dimensions specified at construction
      MDArray value(MDArray& domain, double time_value_optional=0.0)
      {
        MDArray output = getNewCodomain();
        this->operator()(domain, output, time_value_optional);
        return output;
      }

      virtual void operator()(MDArray& domain, MDArray& codomain, double time = 0.0)=0;
      //using GenericFunction::operator();

      // FIXME make protected
      virtual void operator()(MDArray& in, MDArray& out, const stk::mesh::Entity element, const MDArray& parametric_coords, double time_value_optional=0.0)
      {
        EXCEPTWATCH;
        throw std::runtime_error("Not implemented");
      }

      // FIXME make protected
      virtual void operator()(MDArray& in, MDArray& out, const stk::mesh::Bucket& bucket, const MDArray& parametric_coords, double time_value_optional=0.0)
      {
        EXCEPTWATCH;
        throw std::runtime_error("Not implemented");
      }


      /*!
       *  Return a function that is the derivative of this function.  The derivative is specified as a rank-2 array
       *    of strings that specify what derivative to take and how many derivatives.
       *    For example,
       *
       *    \code
       *    Function func(....);
       *    MDArray point(3);
       *
       *    MDArray value(1);
       *    MDArrayString spec_x(1,1);
       *    spec_x(0,0) = "x";
       *    Function deriv = func.derivative(spec_x);
       *    deriv(point, value);
       *
       *    // mixed derivative d^2/dx dy
       *    MDArray value(1);
       *    MDArrayString spec_xy(1,2);
       *    spec_xy(0,0) = "x";
       *    spec_xy(0,1) = "y";
       *    deriv = func.derivative(spec_xy);
       *    deriv(point, value);
       *
       *    // gradient (returns a 3-vector)
       *    MDArray gradient_value(3);
       *    MDArrayString spec_grad(3,1);
       *    spec_grad(0,0) = "x"; 
       *    spec_grad(1,0) = "y"; 
       *    spec_grad(2,0) = "z"; 
       *    deriv = func.derivative(spec_grad);
       *    deriv(point, gradient_value);
       *
       *    // time deriv
       *    MDArrayString spec_t(1,1);
       *    spec_t(0,0) = "t";
       *    deriv = func.derivative(spec_t);
       *    deriv(point, value);
       *  
       *    // Derivative w.r.t. a parameter
       *    ValueFunction vf_A ( "A", 1.234);  // setup a "function" with parameter name "A" with current value 1.234
       *    StringFunction sf_param("A * sin(2.0 * x)");
       *    MDArrayString spec_param(1,1);
       *    spec_param(0,0) = "A";
       *    deriv = func.derivative(spec_param);
       *    deriv(point, value);
       *    // now change A and recompute
       *    vf_A.setValue(2.456);
       *    deriv(point, value);
       *
       *    \endcode
       *
       */

      virtual Teuchos::RCP<Function > derivative(MDArrayString& deriv_spec)
      //virtual Function& derivative(MDArrayString& deriv_spec)
      {
        throw std::runtime_error("not implemented");
      }

      virtual Teuchos::RCP<Function > gradient(int spatialDim=3)
      {
        throw std::runtime_error("not implemented");
      }

      void derivativeAtPoint(MDArrayString& deriv_spec, MDArray& domain, MDArray& codomain, double time = 0.0)
      {
        derivative(deriv_spec)->operator()(domain, codomain, time);
      }


      static void setIntegrationOrderDefault(unsigned integration_order) { s_integration_order_default = integration_order; }
      static void setSpatialDimDefault(unsigned spatialDim) { s_spatialDimDefault = spatialDim; }
      static void setCodomainDimDefault(unsigned codomainDim) { s_codomainDimDefault = codomainDim; }

      /// allow this function to have one or more aliases 
      //    FIXME add ability to delete an alias
      Function * add_alias(const char *alias);

      std::string& getName() { return m_name; }

      void setIntegrationOrder(unsigned iord) { m_integration_order=iord; }
      unsigned getIntegrationOrder(void) { return m_integration_order; }

      void setDomainDimensions(const Dimensions dims);
      void setCodomainDimensions(const Dimensions dims);

      static const Function& Identity;
      //static Function Zero;

      /// Verify that the last dimensions of @param in and @param out are the same; this allows Functions
      ///   to be invoked at multiple points where the first M indices represent an M-d array of points
      ///   to evaluate the function, while the last N indices should match the Functions domain and
      ///   codomain dimensions.
      bool argsAreValid(const MDArray& in, const MDArray& out);

    protected:

      static int last_dimension(MDArray& arr) { return arr.extent_int(arr.rank()-1); }
      static NameToFunctionMap& getNameToFunctionMap();

      std::string m_name;
      unsigned m_integration_order;

    private:

      static NameToFunctionMap s_nameToFunctionMap;
      static unsigned s_integration_order_default;
      static unsigned s_spatialDimDefault;
      static unsigned s_codomainDimDefault;

    };

#ifndef SWIG
    std::ostream &operator<<(std::ostream& out,  Function& func);
#endif

    double eval(double x, double y, double z, double t, Function& func);
    double eval(double x, double y, double z, double t, Teuchos::RCP<Function>& func);

    double eval2(double x, double y, double t, Function& func);
    double eval2(double x, double y, double t, Teuchos::RCP<Function>& func);

    void eval_print(double x, double y, double z, double t, Function& func);
    void eval_print(double x, double y, double z, double t, Teuchos::RCP<Function>& func);

    void eval_print2(double x, double y, double t, Function& func);
    void eval_print2(double x, double y, double t, Teuchos::RCP<Function>& func);

    MDArray eval_vec3(double x, double y, double z, double t, Function& func);
    MDArray eval_vec3(double x, double y, double z, double t, Teuchos::RCP<Function>& func);

    void eval_vec3_print(double x, double y, double z, double t, Function& func);
    void eval_vec3_print(double x, double y, double z, double t, Teuchos::RCP<Function>& func);


  }//namespace percept

#endif
