The classes KokkosUnaryFunction, KokkosBinaryFunction, and KokkosReductionOp (and their derived types),
exist to circumvent a particular difficulty with implementing generic elementwise functions operating
on ROL::Vector types in the context of GPU computing. Currently, ROL uses a virtual inheritance
design for both its Vector types and its elementwise function types. This allows the ROL::Vector
interface to be relatively simple as we need only three member functions to accomodate UnaryFunction, 
BinaryFunction, and ReductionOp type operations. If we need the elementwise logarithm of a vector,
for example, x(i) <- log(x(i)), the class ROL::Elementwise::Logarithm is used in conjunction with 
ROL::Vector::applyUnary, where the implementation details of applying the unary function are specified
in a concrete Vector class, such as ROL::StdVector. 

The problem arises that when one wants to perform an operation of this kind on a GPU using CUDA,
the elementwise function must be defined using the __device__ keyword. ROL is intended to be
agnostic to such details. The classes here have been created to provide a GPU-compute compatible
adapter for ROL through Kokkos without injecting Kokkos dependencies in to the core of ROL's 
source code. The only modification to rol/src needed was to allow ROL's elementwise functions
to accept an abstract Visitor class. This allows for double dispatch and runtime identification
of the specific elementwise functions being passed to ROL::Vector using vtable lookups. Factory
classes are derived from the Visitor types and dynamically allocate Kokkos::DualView compatible
elementwise functions that provide their own implementation details using Kokkos::parallel_for
and Kokkos::parallel_reduce.
