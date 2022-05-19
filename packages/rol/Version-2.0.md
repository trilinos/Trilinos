# ROL Version 2.0

-----------------



## Introduction

__Rapid Optimization Library (ROL) Version 2.0__ offers new interfaces to define
optimization problems and new algorithms to solve them.  The old ROL interfaces and
algorithms will be maintained for approximately 18 months after the initial writing
of this document, until October 31, 2022.
The purpose of this document is to assist existing users of ROL in transitioning
their code from the "old" version, herein denoted by ___Version 1.0___,
to ___Version 2.0___.

#### Contents

1. [Optimization problem](#optimization-problem)
2. [Optimization solver](#optimization-solver)
3. [Algorithms](#algorithms)
4. [Input XML files](#input-xml-files)
5. [Output to stream and string](#output-to-stream-and-string)
6. [Stochastic optimization](#stochastic-optimization)
7. [Avoiding recomputations](#avoiding-recomputations)



## Optimization problem

For user convenience, ROL supports the definition of a "problem," which is subsequently sent to a "solver."  ___Version 2.0___ allows the user to modify the problem through several convenience functions.  Additionally, ___Version 2.0___ enables an explicit specification of **linear constraints**, which can be handled more efficiently by ROL's algorithms.  In contrast, ___Version 1.0___ only supports special treatment of bound constraints.

**Key change**: `ROL::OptimizationProblem` becomes `ROL::Problem`.

#### Basic syntax in Version 1.0

```cpp
    // Instantiate objective function and initial guess vector.
    ROL::Ptr<ROL::Objective<double>> obj = ROL::makePtr<MyObjective<double>>();
    ROL::Ptr<ROL::Vector<double>>      x = ROL::makePtr<MyOptimizationVector<double>>();
    // Instantiate OptimizationProblem.
    ROL::OptimizationProblem<double> problem(obj,x);
```

#### Basic syntax in Version 2.0

```cpp
    // Instantiate objective function and initial guess vector.
    ROL::Ptr<ROL::Objective<double>> obj = ROL::makePtr<MyObjective<double>>();
    ROL::Ptr<ROL::Vector<double>>      x = ROL::makePtr<MyOptimizationVector<double>>();
    // Instantiate Problem.
    ROL::Problem<double> problem(obj,x);
```

#### Added features in Version 2.0

In this example the optimization problem types (TypeU, TypeB, TypeE and TypeG) are as described in [Algorithms](#algorithms).

```cpp
    // TypeU (unconstrained) specification
    ROL::Ptr<ROL::Objective<double>> obj = ROL::makePtr<MyObjective<double>>();
    ROL::Ptr<ROL::Vector<double>>      x = ROL::makePtr<MyOptimizationVector<double>>();
    ROL::Problem<double> problem(obj,x);
    // TypeU can now handle linear equality constraints.  If a linear
    // equality constraint is added to the problem, the TypeU problem
    // will eliminite the equality constraint using the change of variables
    //     x = x_0 + N.apply(y)
    // where x_0 is a feasible point for the linear equality constraint
    // and N is a LinearOperator that projects onto the null space of the
    // constraint Jacobian.
    ROL::Ptr<ROL::Constraint<double>> lin_econ = ROL::makePtr<MyLinearEqualityConstraint<double>>();
    ROL::Ptr<ROL::Vector<double>      lin_emul = ROL::makePtr<MyLinearEqualityConstraintMultiplier<double>>();
    problem.addLinearConstraint("Linear Equality Constraint",lin_econ,lin_mul);

    // TypeB (bound constrained) specification
    ROL::Ptr<ROL::BoundConstraint<double>> bnd = ROL::makePtr<MyBoundConstraint<double>>();
    problem.addBoundConstraint(bnd)    
    // TypeB can now handle polyhedral constraints specified by
    // a bound constraint and linear equality/inequality constraints.
    // If a linear equality/inequality constraint is added to the problem,
    // the TypeB problem will create a PolyhedralProjection object to
    // handle the linear constraints.
    ROL::Ptr<ROL::Constraint<double>>     lin_icon = ROL::makePtr<MyLinearInequalityConstraint<double>>();
    ROL::Ptr<ROL::Vector<double>>         lin_imul = ROL::makePtr<MyLinearInequalityConstraintMultiplier<double>>();
    ROL::Ptr<ROL:BoundConstraint<double>> lin_ibnd = ROL::makePtr<MyLinearInequalityConstraintBound<double>>();
    problem.addLinearConstraint("Linear Inequality Constraint",lin_icon,lin_imul,lin_ibnd);
    // You can set the PolyhedralProjection algorithmic parameters,
    // by calling setProjectionAlgorithm.
    ROL::ParameterList polyProjList;
    ... // fill parameter list with desired polyhedral projection options
    problem.setProjectionAlgorithm(polyProjList);
    
    // TypeG (generally constrained) specification
    ROL::Ptr<ROL::Constraint<double>> econ = ROL::makePtr<MyEqualityConstraint<double>>();
    ROL::Ptr<ROL::Vector<double>>     emul = ROL::makePtr<MyEqualityConstraintMultiplier<double>>();
    problem.addConstraint("Equality Constraint",econ,emul);
    ROL::Ptr<ROL::Constraint<double>>      icon = ROL::makePtr<MyInequalityConstraint<double>>();
    ROL::Ptr<ROL:Vector<double>>          imul = ROL::makePtr<MyInequalityConstraintMultiplier<double>>();
    ROL::Ptr<ROL::BoundConstraint<double>> ibnd = ROL::makePtr<MyInequalityConstraintBound<double>>();
    problem.addConstraint("Inequality Constraint",icon,imul,ibnd);
    
    // Finalize problem (not required, but prints a cool problem
    // summary if called).
    bool lumpConstraints = false;
    bool printToStream   = true;
    std::ostream outStream;
    // If lumpConstraints = false, then any linear constraints will be
    // treated using the TypeU reduction or the TypeB polyhedral projection.
    // Otherwise, the linear constraints will be treated as nonlinear constraints.
    problem.finalize(lumpConstraints,printToStream,outStream)
    
    // Check problem the vector implementations, the derivatives,
    // and the linearity of linear constraints (not required).
    problem.check(printToStream,outStream);
    
    // If the problem has been finalized, but you need to add or
    // remove constraints call the edit function.
    problem.edit();
    
    // TypeE (equality constrained) specification
    // Since we already have equality constraints added, we can remove
    // the inequality constraints to create a TypeE problem.  Like
    // TypeU, Problem will eliminate the linear equality constraints
    // using the change of variables described in the TypeU section.
    problem.removeConstraint("Inequality Constraint");
    problem.removeLinearConstraint("Linear Inequality Constraint");
    
    // Finalize problem again (not required, but prints a cool problem
    // summary if called).
    problem.finalize(lumpConstraints,printToStream,outStream);
```



## Optimization solver

**Key change**: `ROL::OptimizationSolver` becomes `ROL::Solver`.

#### Basic syntax in Version 1.0

```cpp
    // Instantiate Problem.
    ROL::Ptr<ROL::Objective<double>> obj = ROL::makePtr<MyObjective<double>>();
    ROL::Ptr<ROL::Vector<double>>      x = ROL::makePtr<MyOptimizationVector<double>>();
    ROL::OptimizationProblem<double> problem(obj,x);
    ... // add constraints if needed
    
    // Instantiate Solver.
    ROL::ParameterList parlist;
    ... // fill parameter list with desired algorithmic options
    ROL::OptimizationSolver<double> solver(problem,parlist);
    
    // Solve optimization problem.
    std::ostream outStream;
    solver.solve(outStream);
```

#### Basic syntax in Version 2.0

```cpp
    // Instantiate Problem.
    ROL::Ptr<ROL::Objective<double>> obj = ROL::makePtr<MyObjective<double>>();
    ROL::Ptr<ROL::Vector<double>>      x = ROL::makePtr<MyOptimizationVector<double>>();
    ROL::Problem<double> problem(obj,x);
    ... // add constraints if needed
    
    // Finalize Problem (not required, but prints a cool problem
    // summary if you call it).
    bool lumpConstraints = false;
    bool printToStream   = true;
    std::ostream outStream;
    problem.finalize(lumpConstraints,printToStream,outStream);
    
    // Instantiate Solver.
    ROL::ParameterList parlist;
    ... // fill parameter list with desired algorithmic options
    ROL::Solver<double> solver(problem,parlist);
    
    // Solve optimization problem.
    std::ostream outStream;
    solver.solve(outStream);
```



## Algorithms

ROL ___Version 2.0___ maintains a fine-grained interface to directly use specific algorithmic objects.
This is done through a variety of `ROL::Type___::___Algorithm` classes, which replace the `ROL::Step` and
`ROL::Algorithm` classes.
ROL ___Version 2.0___ explicitly categorizes all algorithms into four
groups, based on the problem type:
1. `TypeU`: algorithms for _unconstrained_ problems, e.g., `min f(x)`;
2. `TypeB`: algorithms for _bound-constrained_ problems, e.g., `min f(x) s.t. a <= x <= b`;
3. `TypeE`: algorithms for _equality-constrained_ problems, e.g., `min f(x) s.t. c(x) = 0`; and
4. `TypeG`: algorithms for problems with _general constraints_, e.g., `min f(x) s.t. c(x) >= 0`.

**Key change**: `ROL::Step` and `ROL::Algorithm` are rewritten as `ROL::Type___::___Algorithm`.
For instance, the `ROL::TrustRegionStep` based on the Lin-More method corresponds to `ROL::TypeB::LinMoreAlgorithm`.

Here is an example of ___Version 1.0___ and ___Version 2.0___ usage.

#### Example: Trust-region algorithm for bound-constrained problems

##### Version 1.0

```cpp
    // Instantiate parameter list.
    ROL::ParameterList parlist;
    ... // fill parameter list with desired algorithmic options (e.g., selecting the Lin-More step)
    // Instantiate objective function, initial guess vector and bound constraint.
    ROL::MyObjective<double>        obj;
    ROL::MyVector<double>           x;
    ROL::MyBoundConstraint<double>  bnd;
    // Instantiate step and status test.
    ROL::Ptr<ROL::Step<double>>         step = ROL::makePtr<ROL::TrustRegionStep<double>>(parlist);
    ROL::Ptr<ROL::StatusTest<double>> status = ROL::makePtr<ROL::StatusTest<double>>(parlist);
    // Instantiate algorithm, and run it.
    ROL::Algorithm<double> algo(step,status,false);
    std::outstream outStream;
    algo.run(x, obj, bnd, true, outStream);
```

##### Version 2.0

```cpp
    // Instantiate parameter list.
    ROL::ParameterList parlist;
    ... // fill parameter list with desired algorithmic options
    // Instantiate objective function, initial guess vector and bound constraint.
    ROL::MyObjective<double>        obj;
    ROL::MyVector<double>           x;
    ROL::MyBoundConstraint<double>  bnd;
    // Instantiate the Lin-More trust-region algorithm.
    ROL::TypeB::LinMoreAlgorithm<double> algo(parlist);
    // OPTIONAL: Instantiate and set user-defined status test; otherwise, algo will
    // use a default status test with parameters provided in parlist.
    ROL::Ptr<ROL::StatusTest<double>> status = ROL::makePtr<ROL::MySpecialStatusTest<double>>();
    algo.setStatusTest(status);
    // Run algorithm.
    std::outstream outStream;
    algo.run(x, obj, bnd, outStream);
```

## Output to stream and string

The interface for ___Version 2.0___ can print to an `std::ostream`, but does not return
an `std::vector<std::string>` as in ___Version 1.0___.  To reproduce ___Version 1.0___
behavior, the user can inherit from `std::streambuf`, as follows.

```cpp
#include <iostream>
#include <sstream>

template<class charT=char, class Traits=std::char_traits<charT>>
class mybuffer : public std::basic_streambuf<charT,Traits> {
public:
  const std::stringstream& getStringStream() const { return ss; }

protected:
  inline virtual int overflow(int c = Traits::eof()) {
    if (c != Traits::eof())          ss << static_cast<charT>(c);
    if (putchar(c) == Traits::eof()) return Traits::eof();
    return c;
  }

private:
  std::stringstream ss;
};
```

The user can then build an `std::ostream` using this buffer, i.e.,

```cpp
  mybuffer<> buf;
  std::ostream mystream(&buf);
  ... // Print something to mystream
```
To print the output to file, the user simply grabs the `std::stringstream`
member, i.e.,

```cpp
  std::ofstream file;
  file.open("output.txt");
  if (file.is_open()) {
    file << buf.getStringStream().str();
    file.close();
  }
```


## Input XML files

There are a few minor changes in the naming of steps/algorithms and optimization models.



## Stochastic optimization



## Avoiding recomputations

For computationally expensive objective function and constraint evaluations,
including derivatives, ROL's performance can be enhanced significantly by
taking advantage of `update()` functions.

**Key change**: The second input argument to `update()` has changed.
In ___Version 1.0___, the second argument was a boolean.  In ___Version 2.0___,
the second argument is an enumeration type, which enables finer-grained
control of computations.  The enumeration values and their descriptions are
listed below.

```cpp
enum class UpdateType : std::uint8_t {
  Initial = 0, // Update has not been called before
  Accept,      // This is the new iterate, trial must be called before
  Revert,      // Revert to the previous iterate, trial must be called before
  Trial,       // This is a candidate for the next iterate
  Temp         // For temporary uses including finite difference computations
};
```

Here is a simple example demonstrating the use of the new `update()` function.

```cpp

// Evaluate f(S(x),x) where S(x)=u solves c(u,x) = 0

template<typename Real>
class MyObjective : public ROL::Objective<Real> {
private:
  ... // Member variables
  ROL::Ptr<Objective_SimOpt<Real>> obj_; // Full-space objective
  ROL::Ptr<Constraint_SimOpt<Real>> con_; // Eliminated constraint
  ROL::Ptr<Vector<Real>> u_, ucache_, utemp_; // Storage for state u
  ROL::Ptr<Vector<Real>> r_; // Storage for residual r = c(u,x);

public:
  void update(const ROL::Vector<Real> &x, ROL::UpdateType type, int iter) {
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    if (type == ROL::UpdateType::Initial)  {
      // This is the first call to update
      ucache_ = ROL::makePtr<MyStateVector<Real>>();
      utemp_  = ROL::makePtr<MyStateVector<Real>>();
      u_      = ucache_;
      con_->solve(*r_,*u_,x,tol);
    }
    else if (type == ROL::UpdateType::Accept) {
      // u_ was set to u=S(x) during a trial update
      // and has been accepted as the new iterate
      utemp_  = ucache_;
      ucache_ = u_;      // Cache the accepted value
    }
    else if (type == ROL::UpdateType::Revert) {
      // u_ was set to u=S(x) during a trial update
      // and has been rejected as the new iterate
      u_ = ucache_;      // Revert to cached value
    }
    else if (type == ROL::UpdateType::Trial) {
      // This is a new value of x
      u_ = utemp_;
      con_->solve(*r_,*u_,x,tol);
    }
    else { // ROL::UpdateType::Temp
      // This is a new value of x used for,
      // e.g., finite-difference checks
      u_ = utemp_;
      con_->solve(*r_,*u_,x,tol);
    }
  }
  
  Real value(const ROL::Vector<Real> &x, Real &tol) {
    return obj_->value(*u_,x,tol);
  }
  ... // Define remaining functions
};
```
