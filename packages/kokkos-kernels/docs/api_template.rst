namespace::function
###################

Defined in header: :code:`<namespace>_<function>.hpp`

.. code:: c++

  template <class execution_space, class input_type, class output_type>
  void function(const execution_space& space, const input_type& in, const output_type& out);

  template <class input_type, class output_type>
  void function(const input_type& in, const output_type& out);

Explain what the <function> does and what each overall specificities are.

Parameters
==========

:space: execution space instance

:in: input

:out: output

Type Requirements
-----------------

- `execution_space` must be a Kokkos `execution space <https://kokkos.org/kokkos-core-wiki/API/core/execution_spaces.html>`_

- `input_type` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_

  - ``Kokkos::SpaceAccessibility<execution_space, typename input_type::memory_space>::accessible``
  - ...

- `output_type` must be a Kokkos `View <https://kokkos.org/kokkos-core-wiki/API/core/view/view.html>`_ that satisfies

  - ``Kokkos::SpaceAccessibility<execution_space, typename output_type::memory_space>::accessible``
  - ...

Example
=======

.. code:: cpp

  #include <Kokkos_Core.hpp>
  #include <Namespace_function.hpp>

  using Scalar  = default_scalar;
  using Ordinal = default_lno_t;
  using Offset  = default_size_type;
  using Layout  = default_layout;

  int main(int argc, char* argv[]) {
    Kokkos::initialize(argc, argv);
    {
      ...
    }
    Kokkos::finalize();
    return 0;
  }
