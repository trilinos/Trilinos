// see http://stackoverflow.com/questions/1579007/how-can-i-include-a-subset-of-a-cpp-file-in-a-doxygen-comment

/*! \mainpage STK Percept Documentation
 
  \image html intrepid.png
  \image latex intrepid.jpg "Reconnaissance balloon ``Intrepid''" width=1in

  \section intro_sec Introduction

  %STK Percept is a library and application built on the Sierra Tool Kit (STK).  Percept
  is designed to perform operations on results from Sierra mechanics codes, such as post-processing,
  norms, error estimation, adaptivity, and mesh transfer, thus the name "Percept".

  Percept is based around the concept of a generic function interface that takes multi-dimensional
  arrays as input (the functions' domain objects) and produce multi-dimensional arrays as output
  (codomain).  From this generic interface, users can apply functions to evaluate the fields on
  a mesh, such as pressure or temperature (\ref FieldFunction), or evaluate the gradient (\ref GradientFunction)
  or analytic functions defined by a string expression (\ref StringFunction).
  Multi-dimensional arrays used in Percept are based on mainly the Intrepid FieldContainer
   <a href="http://trilinos.sandia.gov/packages/intrepid/">Intrepid</a> array but also use
   the    <a href="http://trilinos.sandia.gov/packages/shards/">Shards</a> array package.
   Intrepid is used extensively in Percept for:
   - evaluation of mesh Jacobians for geometry verification
   - integration of Functions over mesh domains
   - evaluation of norms of Functions using numerical quadrature
   - evaluation (by interpolation) of field functions
   - support for manufactured solutions
   
  \section overview_sec Overview

  %Percept includes the following features:
   - Definition of analytic functions by providing a string expression, such as "x * sin(y * z)" (StringFunction)
   - Read an Exodus model and results and define FieldFunction for each field
   - add fields to a model (e.g. the gradient) - fields can be nodal or elemental (piecewise constant/linear/etc.)
   - evaluate string or field functions at a given point \c (x,y,z)
   - construct composite functions of StringFunction's by arithmetic operations
   - embedding references to field functions inside a string function
   - computing the norm (L2, H1 semi-norm, etc) of any Function (which includes the difference of two functions,
        for example a field function and a known/manufactured solution)
   - verifying mesh integrity by checking its topologic and geometric consistency/quality
   
  \section future_dev Future Capabilities

  STK Percept is a port and expansion from Sierra Percept (\ref xxx).  Sierra Percept has a wealth of capabilities
  for post-processing, error-estimation, adaptivity, mesh transfere, etc., to name a few.  STK Percept will incorporate
  these features over time, in the rough order of:
   -# uniform mesh refinement for convergence studies
   -# advanced stand-alone executable with new Python-based input syntax
   -# a suite of physics-independent error indicators and error estimators  
   -# tools for physics-dependent error estimation
   -# etc  FIXME
   
  \section quickstart_sec Quick Start

  Familiarity with with the following concepts, objects, and tools is required:
  - multi-dimensional arrays / Intrepid::FieldContainer / \ref md_array_page,
  - <a href="http://trilinos.sandia.gov/packages/shards/">Shards</a> arrays

  Examples of the C++ interface to STK::Percept are shown below.  We start with a simple example of an analytic
  function and evaluate it at a given point.  Then read in a model, construct a FieldFunction and evaluate it.
  Then construct a complex function of the pressure field from a mesh.  Finally, some other examples are shown.

  \subsection ex_string_function String Function
  

  The following can be seen in \see UnitTestStringFunction.cpp unit tests file.

   \clip{UnitTestStringFunction.cpp,start_demo_stringFunction_xy_basic}

  Now a more complicated example that builds a new StringFunction from two others (\see  TEST(function, stringFunction_arithmetic_ops) )

   \clip{UnitTestStringFunction.cpp,start_demo_stringFunction_arithmetic_ops}

  Here's a series of examples on using the derivative capability to evaluate gradient of a StringFunction (note: the methods
  with names derivative_test and derivative_test_fd are only for internal testing).  The first example is very simple and shows
  the derivative being defined for a simple linear function "x-y".  

   \clip{UnitTestStringFunction.cpp,start_demo_stringFunction_derivative}

  The next example shows the gradient of the same function
  
   \clip{UnitTestStringFunction.cpp,start_demo_stringFunction_derivative_1}
   
  while the last derivative example shows the gradient of a more complex function "sin(x*y*z^2)"

   \clip{UnitTestStringFunction.cpp,start_demo_stringFunction_derivative_2}
   
  Here is an example of setting up a FieldFunction and evaluating it at a point

    \clip{UnitTestFieldFunction.cpp, start_demo_fieldFunction_1}

  Here is the same example but now we add a nodal field to contain the magnitude of the coordinates field at each node,
  and use a StringFunction to define and compute the coordinate magnitude at any point.  Finally, we show the interpolation
  of that coordinate magnitude function to the field we created.

    \clip{UnitTestFieldFunction.cpp, start_demo_fieldFunction_2}

  Here is an example following on from the last one that shows a StringFunction referring to a FieldFunction's values:
    \clip{UnitTestFieldFunction.cpp, start_demo_fieldFunction_3}
  \code
     ...
  \endcode
    
  An example of evaluating the difference between a computed pressure field and an analytic function

    \clip{RegressionTestPerceptMeshFieldFunction.cpp, start_demo_open_new_close_PerceptMesh}

  Continuing, we can now compute the norm of the error.  But, first a remark: now that we have introduced Functions, we 
  introduce the concept of operators on functions, under the base class of FunctionOperator.  One type of operator is a norm.  Norms take
  a field function and return a value of the norm.  We consider this to be a special type of more general operation, such
  as a gradient operation which takes a field function and returns another one.  So, we consider the single value returned
  by Norm to be a special type of function, ConstantFunction.

  \code
     Norm<2> l2Norm(bulkData, &metaData.universal_part());
     ConstantFunction norm_result(0.0, "norm_result");
     l2Norm (sf_error, norm_result);
     std::cout << "Norm of error is: " << norm_result.getValue() << std::endl;
  \endcode

  \section math_sec Mathematical Foundations

  Given a domain \f$ \Omega \in R^d \f$, discretized to represent a finite domain \f$ \Omega_h \f$
  with a finite element mesh \f$ T = {\Omega_e} s.t. \Omega_e \cap \Omega_f = \null
 \forall e,f \mbox{and} \cup_e \Omega_e = \Omega_h \f$.  We define \f$ x \f$ to be any point in \f$ \Omega_h \f$
and a function \f$ f : x \rightarrow y \f$ where \f$ y \in Y \f$  where \f$ Y \f$  is a space of results from
the operator of \f$ f \f$ such as a scalar, vector or tensor field.  For example, \f$ Y = R^d X R^d \f$.

  A function \c f has a \e domain and \e codomain.  The domain of \c f is the set of points \c x where it is defined
to be evaluatable.  The codomain is the set of all images \c y of the domain of \c f.  This is not to be confused with
the \e range of \c f, which is a superset of the codomain.  For example, the range of \f$ |x| \f$ is the real line, but
the codomain is just \f$ R^+ \f$.  In Percept the user explicitly specifies the domain and codomain of each function.
 Typically the domain is defined as \f$ R^3 \f$ and the codomain is \f$ R \f$.

  A function operator \c S has as its domain a subset of all functions and its codomain is a subset of all functions.  
  For example, the L2 norm \f$ L^2 \f$ is defined by 
   \f[  L^2 : f \rightarrow c = \sqrt{\int_{\Omega_h} |f(x)|^2 dx } \f]

  \page parallel_node_registry Parallel Node Registry

  The following demo shows use of the NodeRegistry class.  This class enables registration of the need for new nodes
  associated with sub-entities of the elements, for example edges and faces.  This is used for uniform (or other)
  refinement schemes where new nodes are inserted at the mid-points of edges, centroids of faces, etc., and then 
  become part of new elements that break the element into same-topology sub-elements.

  \clip{UnitTestNodeRegistry.cpp, start_demo_nodeRegistry_test_parallel_1}

*/
  
