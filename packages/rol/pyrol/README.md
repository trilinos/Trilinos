# PyROL

PyROL is the Python interface to the
[Rapid Optimization Library (ROL)](https://github.com/sandialabs/rol), a C++
library for large-scale numerical optimization. You write the objective,
constraints and, if needed, the vectors in Python (with NumPy, for example),
and ROL's algorithms solve the problem.

## Installation

```
pip install rol-python
```

Wheels are available for CPython 3.9 to 3.14 on Linux (x86_64) and macOS
(arm64), including free-threaded CPython 3.14.

## Example

Minimize the Rosenbrock function with a trust-region method:

```python
import numpy as np

from pyrol import Objective, Problem, Solver, getCout
from pyrol.pyrol.Teuchos import ParameterList
from pyrol.vectors import NumPyVector


class Rosenbrock(Objective):

    def value(self, x, tol):
        return 100 * (x[1] - x[0] ** 2) ** 2 + (1 - x[0]) ** 2

    def gradient(self, g, x, tol):
        g[0] = -400 * x[0] * (x[1] - x[0] ** 2) - 2 * (1 - x[0])
        g[1] = 200 * (x[1] - x[0] ** 2)


x = NumPyVector(np.array([-3.0, -4.0]))

parameters = ParameterList()
parameters["Step"] = ParameterList()
parameters["Step"]["Type"] = "Trust Region"

solver = Solver(Problem(Rosenbrock(), x), parameters)
solver.solve(getCout())
print(x[:])  # approximately [1. 1.]
```

ROL's algorithms are configured with a parameter list, which you can build in
Python as above or read with `getParametersFromXmlFile` or
`getParametersFromYamlFile`.

## What is available

- `pyrol` exports the commonly used classes for double precision:
  `Objective`, `Constraint`, `LinearOperator`, `LinearConstraint`, `Vector`,
  `Bounds`, `Problem` and `Solver`, plus `getCout`,
  `getParametersFromXmlFile` and `getParametersFromYamlFile`.
  `pyrol.<Name>_forScalar(scalar_type)` looks up the class for a scalar type
  (only `double` is built at the moment).
- `pyrol.unsupported` holds the rest of the wrapped interface, for example
  the SimOpt classes (`Objective_SimOpt`, `Constraint_SimOpt`,
  `Reduced_Objective_SimOpt`, `Vector_SimOpt`).
- `pyrol.vectors.NumPyVector` is a ROL vector backed by a NumPy array.
- Where ROL takes a `std::vector`, pass a Python list, tuple or NumPy array,
  for example `PartitionedVector([x, y])`. The `std::vector`s that ROL
  returns behave like Python lists, and those of numbers can be viewed as
  NumPy arrays without a copy (`numpy.asarray(v)`). A view shares the
  vector's storage, so it is valid only until the vector changes size, for
  example with `append`, `extend` or `resize`; after that, make a new one.

More examples (SimOpt, dynamic constraints, PyTorch tensors and networks) are
in [`pyrol/example`](https://github.com/sandialabs/rol/tree/develop/pyrol/example).

## Free-threaded Python

On free-threaded CPython, importing `pyrol` turns the GIL back on, and Python
prints a `RuntimeWarning` saying so, because PyROL does not yet declare that
it is safe to run without the GIL.

## License

BSD-3-Clause. See [LICENSE](https://github.com/sandialabs/rol/blob/develop/LICENSE).
