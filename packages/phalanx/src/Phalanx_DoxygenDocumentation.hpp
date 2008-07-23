#ifndef PHX_DOXYGEN_DOCUMENTATION_HPP
#define PHX_DOXYGEN_DOCUMENTATION_HPP

/*!

\mainpage

*** WARNING: DOCUMENTATION IS UNDER CONSTRUCTION ***

This will be updated for release 9.0

\section contents Contents

  - \ref introduction

  - \ref overview

  - \ref user_guide

  - \ref developer_guide

  - \ref faq

  - \ref history

  - \ref authors

\section introduction Introduction

Phalanx is a library to handle arbitrary function evaluation with complex nonlinear dependency chains for discretized partial differential equation systems.  It provides a flexible and efficient mechanism for switching dependencies and computing sensitivities either analytically or via automatic differentiation.  It can be used with any cell-based discretization technique including finite element, finite volume and finite difference.

\section overview Overview

Phalanx is used to break down complex dependencies into manageable algorithmic blocks.  It was written to address a variety of difficulties encountered when writing a general PDE code that allows for flexible switching of dependencies at runtime.  

Definitions:
 - A "Field" is a concrete data type that has allocated storage space.  A field can be defined for a number of cells and is specific to a DataLayout.
 - A "Data Layout" is a description of a unique topologic entity used for defining where a field lives.  The API for a data layout contains a size (number of entities on a cell) and allows for a comparison operator.
 - An "Evaluator" will evaluate the values of fields and is in turn dependent on other fields.

A simple example (found in phalanx/example/energyFlux) is the construction of a Fourier energy flux for the heat equation.  Suppose that we want to solve the heat equation using finite elements:

\f[
  \nabla \cdot (\mathbf{q}) + s = 0
\f] 

where \f$\mathbf{q}\f$ is the heat flux, \f$ \mathbf{q} = -\rho C_p \nabla T \f$  and \f$s\f$ is a nonlinear source term.  Integrating by parts and using the test function \f$ \phi \f$ results in the following weak form equation:

\f[
  - \int_{\Omega} \nabla \phi \cdot \mathbf{q} d\Omega 
  + \int_{\Gamma} \phi \mathbf{n} \cdot \mathbf{q} 
  + \int_{\Omega} \phi s d\Omega = 0 
\f]




\f[
  \mathbf{q} = -\rho C_p \nabla T
\f]

Advantages of using Phalanx:

 - Increased flexibility because each simpler piece becomes an extension point that can be swapped out with different implementations.

 - Easier to implement new code because each piece is simpler, more focused and easier to test in isolation.  

 - Ensures that complex dependency chains are evaluated correctly.  

 - Determine the order of variable evaluation so that the correct
 evaluation is performed: For example, to compute the density at the
 quadrature points, it could be a constant, it could be a function of
 temperature or it could be a function of temperature and pressure.  A
 factory creates the density provider and from this, the manager
 automatically figures out the correct evaluation order - whether T
 and/or P will have to be evaluated before the density.  This is an
 easy example, but there are much more complex dependency chains.

Disadvantages of using Phalanx:

 - A potential performance loss due to fragmentation of the over-all algorithm (e.g., several small loops instead of one large loop) 

 - A potential loss of visibility of the original, composite problem (since the code is scattered into multiple places).  

Managing these trade-offs can result in application code that both performs well and supports rapid development and extensibility.  

\section user_guide User's Guide

Phalanx is distributed as a package in the <a href="http://trilinos.sandia.gov">Trilinos Framework</a>.  It has direct dependencies on the following third party libraries:

 - Requires the <a href="http://trilinos.sandia.gov/packages/teuchos">Teuchos</a> utilities library, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.
 
 - Requires the <a href="http://trilinos.sandia.gov/packages/sacado">Sacado Automatic Differentiation Library</a>, part of the <a href="http://trilinos.sandia.gov/">Trilinos Framework</a>.

 - Requires the <a href="http://www.boost.org">Boost Template Metaprogramming (MPL) Library</a>


\section developer_guide Developer's Guide

\section faq Frequently Asked Questions

  - Why name it Phalanx?  The phalanx was one of the most dominant military formations of the Greek armies in the classical period.  It was a strictly ordered formation.  The Phalanx software library figures out ordered dependencies of field evaluators.   A second more obscure reference relates to the US Navy.  The Phalanx software package was designed to provide nonlinear functionality for the <a href="http://trilinos.sandia.gov/packages/intrepid">Intrepid</a> itegration library.  Intrepid was the name of an aircraft carrier during in WW II.  Modern US aircraft carriers are protected by a Close-In Weapons System (CIWS) named the Phalanx.  Finally, the PI of this project is an avid strategy warfare gamer and leans towards military references.

  - Why is the scalar type not embedded in the FieldTag?  We designed the FieldTag to be scalar type independent.  The reason is that we could then use FieldTags as arguments for Evaluator objects of ANY scalar type.  We have a factory that can automatically build Evaluators for all scalar types.  This automation requires constructor arguments that are not dependent on the scalar type.  

\section history History

Phalanx grew out of the Variable Manger in the Charon code and the Expression Manager in the Aria code.  It is an attempt at merging the two capabilities and is slated to provide nonlinear function evaluation to the Intrepid discretiztion package.

\section authors Authors

  - Roger Pawlowski (PI), SNL 01414

*/

/* ************************************************************************ */
/* ************************************************************************ */

/*! \page developer_notes Developer Notes

\section dnotes_notz Description

This is a description Pat Notz called his elevator speech on Phalanx:

\verbatim
The {NAME} package is a library to help decompose a complex problem
into a set of simpler problems with managed dependencies. Two benefits
of decomposing a problem using {NAME} are (1) increased flexibility
because each simpler piece becomes an extension point that can be
swapped out with different implementations and (2) easier to craft
code because each piece is simpler, more focused and easier to test in
isolation.  The counter-benefits to these are (1) a potential
performance loss due to fragmentation of the over-all algorithm (e.g.,
several small loops instead of one large loop) and (2) a potential
loss of visibility of the original, composite problem (since the code
is scattered into multiple places).  Managing these trade-offs can
result in application code that both performs well and supports rapid
development and extensibility.  
\endverbatim

\section dnotes_contiguous_allocator Contiguous Allocator

Simple concept for chunk-layout of heterogeneous members with proper
data alignment.

A 'malloc'ed chunk of memory is always suitable to hold any data type,
or array of data types (do a 'man malloc').

\verbatim
        unsigned char * chunk = malloc( chunk_size );
\endverbatim

If a 'double' is to occupy a portion of the chunk then is must be
aligned such that:

\verbatim
        double * d = (double *)( chunk + offset_d );
\endverbatim

where:

\verbatim
        0 == offset_d % sizeof(double);
\endverbatim

The key point is for the offset of each member within the chunk to
have the correct alignment, and of course not overlap any other
member.

\section dnotes_old_description Old email with description of variable manager 

Email description of Phalanx (formerly known as the variable manager):

Phalanx is a tool for controlling the evaluation of
data used in the assembly of a residual/Jacobian fill for PDE solvers.
I implemented the one in Charon and one also exists in Aria (although it
is used for a different purpose).  The idea behind the VM is to allow
for fast changes to dependency trees.  The VM provides variable
quantities.  A variable consists of a name (string), a type (scalar,
vector, tensor - all templated for ad), and the place that it lives
(quadrature point, node, edge, face, element, etc...).  A variable is
evaluated using a provider.  The provider contains lists of variables it
provides (evlauates) and a list of variables it requires to perform the
evaluation (its dependencies).  The VM figures out the correct order to
call the providers to get the dependency chain correct.

I developed this for Charon to address the following issues:

1. Determine the order of variable evaluation so that the correct
evaluation is performed: For example, to compute the density, it could
be a constant, it could be a function of temperature or it could be a
function of temperature and pressure.  A factory creates the density
provider and from this, the manager automatically figures out the
correct evaluation order - whether T and/or P will have to be
evaluated before the density.  This is an easy example, but there are
much more complex dependency chains.

2. Provides a single place to store all evaluation variables so that the
procedure by which the variable is evaluated is not known by the objects
using the variable.  For example, in Charon a quantity of interest is
the velocity components.  The velocity could be a degree of freedom or
it could be a user specified constant, or it could be a user specified
field (read from an exodus file).  Before we had the variable manager,
everywhere in Charon that used velocity, there were repeated loops that
looked like:

\verbatim
if (velocity is DOF)
  // get values from dof structure
else if (velocity is auxiliary data)
  // go get from aux data structure
else if (...)
\endverbatim

There were hard coded loops in multiple places that had to be updated
every time we added a new way to specify velocity.   Now, when velocity
is required, you just go to the VM and get it.  No questions asked.

3. By having a single place to store the variables, we eliminated the
reevaluation of temporary quantities that are used in different PDE
operators and different equations.  Prior to using the variable manager,
the reuse of variables like density was ignored. It was recomputed for
each operator evaluation.  After implementing the variable manager, my
residual evaluation times were reduced by a factor of 7 for reacting
flow problems!

4. By having the evaluation order, one can perform automatic
differentiation without any library like saccado.  When the user writes
a provider for a function, they can also write a local Jacobian for the
individual function.  Then by using the chain rule, the local element
Jacobian can be assembled by traversing the dependency tree generated by
the variable manager.  Ariawrote their VM to do this, since Sacado was
not around when they started.  Aria does use Sacado for user defined
functions, but uses the VM for internal functions.

There are other benefits, but you get the idea.  A primary goal of my
current ASC funding is to generalize the variable manager in Charon
(there are many things I can do better now to make it more efficient and
easier to use) as a trilinos package to be used with Intrepid to aid in
element assembly for nonlinear equation sets.

  */

/* ************************************************************************ */
/* ************************************************************************ */

/*! \page junk Junk

  \todo Add a configure check for BOOST.

  \todo Add a configure check for TVMET.

  \todo Add documentation for tvmet.

*/

/* ************************************************************************ */
/* ************************************************************************ */

#endif
