// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

# Belos Tpetra Preconditioning in Stratimikos

The adapters in this directory allow a (unpreconditioned) Belos solver using Tpetra linear algebra to be called from a Stratimikos input deck as a preconditioner to another Belos solver.  

This can be used, for example, to apply Belos' polynomial preconditioner using Stratimikos.
It can also be used with the fixed-point solver as the outer solver to perform iterative refinement around an inner solver. 
(See example in...)



## Usage

To enable the Belos-as-preconditioner functionality, simply include the following in your Stratimikos code: 

```
#include <Stratimikos_BelosTpetraPrecHelpers.hpp>
```

Then, after creating your Stratimikos linearSolverBuilder, enable Belos-as-preconditioner as follows:

```
// Create Stratimikos::LinearSolverBuilder
Stratimikos::LinearSolverBuilder<Scalar> linearSolverBuilder;

// Register Belos+Tpetra as preconditioner:
// (Note: One could use Tpetra::Operator instead of Tpetra::CrsMatrix.)
Stratimikos::enableBelosPrecTpetra<Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>>(linearSolverBuilder);
```

## Examples

One can find the following example parameter list inputs in `trilinos/packages/stratimikos/test/`

+ list1.xml
+ list2.xml
+ list3.xml
