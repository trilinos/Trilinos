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
5. [Stochastic optimization](#stochastic-optimization)
6. [Avoiding recomputations](#avoiding-recomputations)



## Optimization problem

Describe the need for changes at a high level.

#### Basic syntax in Version 1.0

```cpp
    ROL::OptimizationProblem() problem;
```

#### Basic syntax in Version 2.0

```cpp
    ROL::Problem() problem;
```

#### Added features in Version 2.0

```cpp
    ROL::Problem() problem;
    problem.edit();
    problem.addConstraint();
    problem.finalize();
```



## Optimization solver

#### Basic syntax in Version 1.0

```cpp
    ROL::OptimizationSolver() solver;
```

#### Basic syntax in Version 2.0

```cpp
    ROL::Solver() solver;
```



## Algorithms

ROL Version 2.0 maintains a fine-grained interface to directly use
specific algorithmic objects.  This is done through a variety of
`ROL::Algorithm` classes, which replace the `ROL::Step` classes.
ROL Version 2.0 explicitly categorizes all algorithms into four
groups, based on the problem type:
1. `TypeU`: algorithms for _unconstrained_ problems, e.g., `min f(x)`;
2. `TypeB`: algorithms for _bound-constrained_ problems, e.g., `min f(x) s.t. a <= x <= b`;
3. `TypeE`: algorithms for _equality-constrained_ problems, e.g., `min f(x) s.t. c(x) = 0`; and
4. `TypeG`: algorithms for problems with _general constraints_, e.g., `min f(x) s.t. c(x) >= 0`.

Here are a few examples of Version 1.0 and Version 2.0 usage.

#### Example 1: Trust-region algorithm for unconstrained problems

##### Version 1.0

```cpp
    ROL::TrustRegionStep() step;
    more_detail();
```

##### Version 2.0

```cpp
    ROL::TypeU::TrustRegionAlgorithm() algorithm;
    more_detail();
```

#### Example 2: Trust-region algorithm for bound-constrained problems

##### Version 1.0

```cpp
    ROL::TrustRegionStep() step;
    parlist_selection();
    more_detail();
```

##### Version 2.0

```cpp
    ROL::ParameterList parlist;
    ... // fill parameter list with desired algorithmic options
    ROL::TypeB::LinMoreAlgorithm<double> algo(parlist);
    ... // define ROL::Vector 'x', ROL::Objective 'f' and ROL::BoundConstraint 'bnd'
    algo.run(x, f, bnd);
```



## Input XML files

There are minor changes in the naming of steps/algorithms and models.



## Stochastic optimization



## Avoiding recomputations

For computationally expensive objective function and constraint evaluations,
including derivatives, ROL's performance can be enhanced significantly by
taking advantage of `update()` functions.

