// @HEADER
// *****************************************************************************
//          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
//
// Copyright 2022 NTESS and the PyTrilinos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PYTRILINOS2_TPETRA_CUSTOM
#define PYTRILINOS2_TPETRA_CUSTOM

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <Tpetra_CrsGraph_decl.hpp>
#include <Tpetra_CrsGraph_def.hpp>

#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_CrsMatrix_def.hpp>

namespace py = pybind11;

template<typename T>
Teuchos::ArrayView< T > copyNumPyToTeuchosArrayView(pybind11::array_t<T> array) {
    auto np_array = array.template mutable_unchecked<1>();
    int size = array.shape(0);
    T * vec = new T[size];
    for (int i=0; i < size; ++i)
      vec[i] = np_array(i);
    Teuchos::ArrayView< T > av(vec, size);
    return av;
}

template<typename T>
Teuchos::ArrayView< const T > copyNumPyToTeuchosConstArrayView(pybind11::array_t<T> array) {
    return copyNumPyToTeuchosArrayView<T>(array).getConst();
}

// ----------------

// The implementation of the conversion from numpy array to Kokkos view
// in both directions is based on:
// https://github.com/sandialabs/compadre/blob/master/pycompadre/pycompadre.cpp

namespace py = pybind11;

template< bool B, class T = void >
using enable_if_t = typename std::enable_if<B,T>::type;

template<typename T>
Teuchos::ArrayView< T > convert_np_to_ArrayView(pybind11::array_t<T> array) {

    int size = array.shape(0);
    Teuchos::ArrayView< T > av(array.mutable_data(0), size);

    return av;
}

// conversion of numpy arrays to kokkos arrays

template<typename ViewType>
void convert_np_to_kokkos_1d(pybind11::array_t<typename ViewType::non_const_value_type> array,  ViewType kokkos_array_device) {

    auto np_array = array.template unchecked<1>();

    auto kokkos_array_host = Kokkos::create_mirror_view(kokkos_array_device);
    //Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, array.shape(0)), [&](int i) {
      for (int i=0; i<array.shape(0); ++i) {
        kokkos_array_host(i) = np_array(i);
      }
    //});
    Kokkos::fence();
    Kokkos::deep_copy(kokkos_array_device, kokkos_array_host);
}

template<typename ViewType>
void convert_np_to_kokkos_2d(pybind11::array_t<typename ViewType::non_const_value_type> array,  ViewType kokkos_array_device) {

    auto np_array = array.template unchecked<2>();

    auto kokkos_array_host = Kokkos::create_mirror_view(kokkos_array_device);
    //Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, array.shape(0)), [&](int i) {
      for (int i=0; i<array.shape(0); ++i) {
        for (int j=0; j<array.shape(1); ++j) {
            kokkos_array_host(i,j) = np_array(i,j);
        }
      }
    //});
    Kokkos::fence();
    Kokkos::deep_copy(kokkos_array_device, kokkos_array_host);
}

// conversion of kokkos arrays to numpy arrays

template<typename T, typename T2=void>
struct cknp1d {
    pybind11::array_t<typename T::value_type> result;
    cknp1d (T kokkos_array_host) {

        const int dim_out_0 = kokkos_array_host.extent(0);
        result = pybind11::array_t<typename T::value_type>(dim_out_0);
        auto data = result.template mutable_unchecked<1>();
        //Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0,dim_out_0), [&](int i) {
          for (int i=0; i<dim_out_0; ++i) {
            data(i) = kokkos_array_host(i);
          }
        //});
        Kokkos::fence();

    }
    pybind11::array_t<typename T::value_type> convert() { return result; }
};

template<typename T>
struct cknp1d<T, enable_if_t<(T::rank!=1)> > {
    pybind11::array_t<typename T::value_type> result;
    cknp1d (T kokkos_array_host) {
        result = pybind11::array_t<typename T::value_type>(0);
    }
    pybind11::array_t<typename T::value_type> convert() { return result; }
};

template<typename T, typename T2=void>
struct cknp2d {
    pybind11::array_t<typename T::value_type> result;
    cknp2d (T kokkos_array_host) {

        const int dim_out_0 = kokkos_array_host.extent(0);
        const int dim_out_1 = kokkos_array_host.extent(1);

        result = pybind11::array_t<typename T::value_type>(dim_out_0*dim_out_1);
        result.resize({dim_out_0,dim_out_1});
        auto data = result.template mutable_unchecked<T::rank>();
        //Kokkos::parallel_for(Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0,dim_out_0), [&](int i) {
          for (int i=0; i<dim_out_0; ++i) {
            for (int j=0; j<dim_out_1; ++j) {
                data(i,j) = kokkos_array_host(i,j);
            }
          }
        //});
        Kokkos::fence();

    }
    pybind11::array_t<typename T::value_type> convert() { return result; }
};

template<typename T>
struct cknp2d<T, enable_if_t<(T::rank!=2)> > {
    pybind11::array_t<typename T::value_type> result;
    cknp2d (T kokkos_array_host) {
        result = pybind11::array_t<typename T::value_type>(0);
    }
    pybind11::array_t<typename T::value_type> convert() { return result; }
};


template<typename T>
pybind11::array_t<typename T::value_type> convert_kokkos_to_np(T kokkos_array_device) {

    // ensure data is accessible
    auto kokkos_array_host =
        Kokkos::create_mirror_view(kokkos_array_device);
    Kokkos::deep_copy(kokkos_array_host, kokkos_array_device);

    pybind11::array_t<typename T::value_type> result;
    if (T::rank==1) {
        result = cknp1d<decltype(kokkos_array_host)>(kokkos_array_host).convert();
    } else if (T::rank==2) {
        result = cknp2d<decltype(kokkos_array_host)>(kokkos_array_host).convert();
    } else {
        result = pybind11::array_t<typename T::value_type>(0);
    }
    return result;

}

// ----------------

template<typename SCALAR, typename LO, typename GO, typename NODE>
pybind11::array_t<SCALAR> getLocalViewHost(Teuchos::RCP<Tpetra::Vector<SCALAR,LO,GO,NODE>> &vector) {
    return convert_kokkos_to_np(Kokkos::subview(vector->getLocalViewDevice(Tpetra::Access::ReadOnly), Kokkos::ALL, 0));
}

template<typename SCALAR, typename LO, typename GO, typename NODE>
pybind11::array_t<SCALAR> getLocalViewHost(Teuchos::RCP<Tpetra::MultiVector<SCALAR,LO,GO,NODE>> &mvector) {
    return convert_kokkos_to_np(mvector->getLocalViewDevice(Tpetra::Access::ReadOnly));
}

template<typename SCALAR, typename LO, typename GO, typename NODE>
void setLocalViewHost(Teuchos::RCP<Tpetra::Vector<SCALAR,LO,GO,NODE>> &vector, pybind11::array_t<SCALAR> input) {
    auto view = Kokkos::subview(vector->getLocalViewDevice(Tpetra::Access::ReadWrite), Kokkos::ALL, 0);
    convert_np_to_kokkos_1d(input, view);
}

template<typename SCALAR, typename LO, typename GO, typename NODE>
void setLocalViewHost(Teuchos::RCP<Tpetra::MultiVector<SCALAR,LO,GO,NODE>> &mvector, pybind11::array_t<SCALAR> input) {
    auto view = mvector->getLocalViewDevice(Tpetra::Access::ReadWrite);
    convert_np_to_kokkos_2d(input, view);
}

template<typename T>
void define_CrsGraph_member_functions(T cl) {
  using LO = typename T::type::local_ordinal_type;
  using GO = typename T::type::global_ordinal_type;
  using NODE = typename T::type::node_type;
  cl.def("insertGlobalIndices", [](Teuchos::RCP<Tpetra::CrsGraph<LO,GO,NODE>> &m,
  const GO row,
  pybind11::array_t<GO> cols) {
    m->insertGlobalIndices(row, copyNumPyToTeuchosArrayView(cols));
  }, "Insert global indices into the graph.\n\n \n  is a valid index in the row Map.  It need\n   not be owned by the calling process.\n \n\n isLocallyIndexed() == false\n \n\n isStorageOptimized() == false\n\n \n indicesAreAllocated() == true\n \n\n isGloballyIndexed() == true\n\n If  does not belong to the graph on this process,\n then it will be communicated to the appropriate process when\n globalAssemble() is called.  (That method will be called\n automatically during the next call to fillComplete().)\n Otherwise, the entries will be inserted into the part of the\n graph owned by the calling process.\n\n If the graph row already contains entries at the indices\n corresponding to values in  then the redundant\n indices will be eliminated.  This may happen either at\n insertion or during the next call to fillComplete().\n\nC++: Tpetra::CrsGraph<int, long long, Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> >::insertGlobalIndices(const long long, const class Teuchos::ArrayView<const long long> &) --> void", pybind11::arg("gblRow"), pybind11::arg("inputGblColInds"));
  cl.def("insertLocalIndices", [](Teuchos::RCP<Tpetra::CrsGraph<LO,GO,NODE>> &m,
  const LO row,
  pybind11::array_t<LO> cols) {
    m->insertLocalIndices(row, copyNumPyToTeuchosArrayView(cols));
  }, "Insert local indices into the graph.\n\n       \n  is a local row belonging to the graph on this process.\n       \n\n isGloballyIndexed() == false\n       \n\n isStorageOptimized() == false\n       \n\n hasColMap() == true\n\n       \n indicesAreAllocated() == true\n       \n\n isLocallyIndexed() == true\n\n       \n If the graph row already contains entries at the indices\n         corresponding to values in  then the redundant\n         indices will be eliminated; this may happen at insertion or\n         during the next call to fillComplete().\n\nC++: Tpetra::CrsGraph<int, long long, Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> >::insertLocalIndices(const int, const class Teuchos::ArrayView<const int> &) --> void", pybind11::arg("localRow"), pybind11::arg("indices"));
}

template<typename T>
void define_CrsMatrix_member_functions(T cl) {
  using SCALAR = typename T::type::scalar_type;
  using LO = typename T::type::local_ordinal_type;
  using GO = typename T::type::global_ordinal_type;
  using NODE = typename T::type::node_type;

  cl.def("insertGlobalValues", [](Teuchos::RCP<Tpetra::CrsMatrix<SCALAR,LO,GO,NODE>> &m,
  const GO row,
  pybind11::array_t<GO> cols,
  pybind11::array_t<SCALAR> vals) {
    m->insertGlobalValues(row, copyNumPyToTeuchosArrayView(cols), copyNumPyToTeuchosArrayView(vals));
  }, "If the matrix has a column Map (hasColMap() == true),\n and if globalRow is owned by process p, then it is forbidden\n to insert column indices that are not in the column Map on\n process p.  Tpetra will test the input column indices to\n ensure this is the case, but if  is not owned by\n the calling process, the test will be deferred until the next\n call to globalAssemble() or fillComplete().\n\n \n The behavior described in the above paragraph differs\n   from that of Epetra.  If the matrix has a column Map,\n   Epetra_CrsMatrix \"filters\" column indices not in the column\n   Map.  Many users found this confusing, so we changed it so\n   that nonowned column indices are forbidden.\n\n It is legal to call this method whether the matrix's column\n indices are globally or locally indexed.  If the matrix's\n column indices are locally indexed (isLocallyIndexed() ==\n true), then this method will convert the input global\n column indices to local column indices.\n\n For better performance when filling entries into a sparse\n matrix, consider the following tips:\n \n Use local indices (e.g., insertLocalValues()) if you know\n   the column Map in advance.  Converting global indices to\n   local indices is expensive.  Of course, if you don't know\n   the column Map in advance, you must use global indices.\n When invoking the CrsMatrix constructor, give the best\n   possible upper bounds on the number of entries in each row\n   of the matrix.  This will avoid expensive reallocation if\n   your bound was not large enough.\n If you plan to reuse a matrix's graph structure, but\n   change its values, in repeated fillComplete() / resumeFill()\n   cycles, you can get the best performance by creating the\n   matrix with a const CrsGraph.  Do this by using the\n   CrsMatrix constructor that accepts an RCP of a const\n   CrsGraph.  If you do this, you must use the \"replace\" or\n   \"sumInto\" methods to change the values of the matrix; you\n   may not use insertGlobalValues() or\n   insertLocalValues().\n \n\nC++: Tpetra::CrsMatrix<double, int, long long, Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> >::insertGlobalValues(const long long, const class Teuchos::ArrayView<const long long> &, const class Teuchos::ArrayView<const double> &) --> void", pybind11::arg("gblRow"), pybind11::arg("indices"), pybind11::arg("values"));
  cl.def("insertLocalValues", [](Teuchos::RCP<Tpetra::CrsMatrix<SCALAR,LO,GO,NODE>> &m,
  const LO row,
  pybind11::array_t<LO> cols,
  pybind11::array_t<SCALAR> vals) {
    m->insertLocalValues(row, copyNumPyToTeuchosArrayView(cols), copyNumPyToTeuchosArrayView(vals));
  }, "", pybind11::arg("lclRow"), pybind11::arg("indices"), pybind11::arg("values"));
  cl.def("insertLocalValues", [](Teuchos::RCP<Tpetra::CrsMatrix<SCALAR,LO,GO,NODE>> &m,
  const LO row,
  pybind11::array_t<LO> cols,
  pybind11::array_t<SCALAR> vals,
  const enum Tpetra::CombineMode mode) {
    m->insertLocalValues(row, copyNumPyToTeuchosArrayView(cols), copyNumPyToTeuchosArrayView(vals), mode);
  }, "Insert one or more entries into the matrix, using local\n   column indices.\n\n \n [in] Local index of the row into which to\n   insert the entries.  It must be owned by the row Map on the\n   calling process.\n \n\n [in] Local indices of the columns into which to\n   insert the entries.  All of the column indices must be owned\n   by the column Map on the calling process.\n \n\n [in] Values to insert into the above columns.\n \n\n [in] How values should be inserted. Valid options\n   are: ADD (default) inserts values that are not yet in the\n   matrix graph, and sums values that are already present. INSERT\n   inserts values that are not yet in the matrix graph, and\n   replaces values that are already present.\n\n For all k in 0, ..., cols.size()-1, insert the value\n values[k] into entry (globalRow, cols[k]) of\n the matrix.  If that entry already exists, add the new value\n to the old value, if CM=ADD, otherwise replace\n the old value.\n\n In order to call this method, the matrix must be locally\n indexed, and it must have a column Map.\n\n For better performance when filling entries into a sparse\n matrix, consider the following tips:\n \n When invoking the CrsMatrix constructor, give the best\n   possible upper bounds on the number of entries in each row\n   of the matrix.  This will avoid expensive reallocation if\n   your bound was not large enough.\n If you plan to reuse a matrix's graph structure, but\n   change its values, in repeated fillComplete() / resumeFill()\n   cycles, you can get the best performance by creating the\n   matrix with a const CrsGraph.  Do this by using the\n   CrsMatrix constructor that accepts an RCP of a const\n   CrsGraph.  If you do this, you must use the \"replace\" or\n   \"sumInto\" methods to change the values of the matrix; you\n   may not use insertGlobalValues() or\n   insertLocalValues().\n \n\nC++: Tpetra::CrsMatrix<double, int, long long, Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> >::insertLocalValues(const int, const class Teuchos::ArrayView<const int> &, const class Teuchos::ArrayView<const double> &, const enum Tpetra::CombineMode) --> void", pybind11::arg("lclRow"), pybind11::arg("indices"), pybind11::arg("values"), pybind11::arg("CM"));
  cl.def("replaceGlobalValues", [](Teuchos::RCP<Tpetra::CrsMatrix<SCALAR,LO,GO,NODE>> &m,
  const GO row,
  pybind11::array_t<GO> cols,
  pybind11::array_t<SCALAR> vals) {
    m->replaceGlobalValues(row, copyNumPyToTeuchosArrayView(cols), copyNumPyToTeuchosArrayView(vals));
  }, "Overload of replaceGlobalValues (see above), that takes\n   Teuchos::ArrayView (host pointers) instead of Kokkos::View.\n\nC++: Tpetra::CrsMatrix<double, int, long long, Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> >::replaceGlobalValues(const long long, const class Teuchos::ArrayView<const long long> &, const class Teuchos::ArrayView<const double> &) --> int", pybind11::arg("globalRow"), pybind11::arg("inputGblColInds"), pybind11::arg("inputVals"));
  cl.def("replaceLocalValues", [](Teuchos::RCP<Tpetra::CrsMatrix<SCALAR,LO,GO,NODE>> &m,
  const LO row,
  pybind11::array_t<LO> cols,
  pybind11::array_t<SCALAR> vals) {
    m->replaceLocalValues(row, copyNumPyToTeuchosArrayView(cols), copyNumPyToTeuchosArrayView(vals));
  }, "Backwards compatibility version of replaceLocalValues\n   (see above), that takes Teuchos::ArrayView (host pointers)\n   instead of Kokkos::View.\n\nC++: Tpetra::CrsMatrix<double, int, long long, Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> >::replaceLocalValues(const int, const class Teuchos::ArrayView<const int> &, const class Teuchos::ArrayView<const double> &) --> int", pybind11::arg("localRow"), pybind11::arg("lclCols"), pybind11::arg("vals"));
  cl.def("sumIntoGlobalValues", [](Teuchos::RCP<Tpetra::CrsMatrix<SCALAR,LO,GO,NODE>> &m,
  const GO row,
  pybind11::array_t<GO> cols,
  pybind11::array_t<SCALAR> vals) {
    m->sumIntoGlobalValues(row, copyNumPyToTeuchosArrayView(cols), copyNumPyToTeuchosArrayView(vals));
  }, "", pybind11::arg("gblRow"), pybind11::arg("inputGblColInds"), pybind11::arg("inputVals"));
  cl.def("sumIntoGlobalValues", [](Teuchos::RCP<Tpetra::CrsMatrix<SCALAR,LO,GO,NODE>> &m,
  const GO row,
  pybind11::array_t<GO> cols,
  pybind11::array_t<SCALAR> vals,
  const bool atomic) {
    m->sumIntoGlobalValues(row, copyNumPyToTeuchosArrayView(cols), copyNumPyToTeuchosArrayView(vals), atomic);
  }, "Sum into one or more sparse matrix entries, using\n   global indices.\n\n This is a local operation; it does not involve communication.\n However, if you sum into rows not owned by the calling\n process, it may result in future communication in\n globalAssemble() (which is called by fillComplete()).\n\n If  is owned by the calling process, then this\n method performs the sum-into operation right away.  Otherwise,\n if the row is not owned by the calling process, this\n method defers the sum-into operation until globalAssemble().\n That method communicates data for nonowned rows to the\n processes that own those rows.  Then, globalAssemble() does\n one of the following:\n \n  It calls insertGlobalValues() for that data if the matrix\n      has a dynamic graph. \n  It calls sumIntoGlobalValues() for that data if the matrix\n      has a static graph.  The matrix silently ignores\n      (row,column) pairs that do not exist in the graph.\n \n\n \n [in] The global index of the row in which to\n   sum into the matrix entries.\n \n\n [in] One or more column indices.\n \n\n [in] One or more values corresponding to those\n   column indices.  vals[k] corresponds to\n   cols[k].\n \n\n [in] Whether to use atomic updates.\n\n \n The number of indices for which values were actually\n   modified; the number of \"correct\" indices.\n\n This method has the same preconditions and return value\n meaning as replaceGlobalValues() (which see).\n\nC++: Tpetra::CrsMatrix<double, int, long long, Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> >::sumIntoGlobalValues(const long long, const class Teuchos::ArrayView<const long long> &, const class Teuchos::ArrayView<const double> &, const bool) --> int", pybind11::arg("gblRow"), pybind11::arg("inputGblColInds"), pybind11::arg("inputVals"), pybind11::arg("atomic"));
  cl.def("sumIntoLocalValues", [](Teuchos::RCP<Tpetra::CrsMatrix<SCALAR,LO,GO,NODE>> &m,
  const LO row,
  pybind11::array_t<LO> cols,
  pybind11::array_t<SCALAR> vals) {
    m->sumIntoLocalValues(row, copyNumPyToTeuchosArrayView(cols), copyNumPyToTeuchosArrayView(vals));
  }, "", pybind11::arg("localRow"), pybind11::arg("indices"), pybind11::arg("values"));
  cl.def("sumIntoLocalValues", [](Teuchos::RCP<Tpetra::CrsMatrix<SCALAR,LO,GO,NODE>> &m,
  const LO row,
  pybind11::array_t<LO> cols,
  pybind11::array_t<SCALAR> vals,
  const bool atomic) {
    m->sumIntoLocalValues(row, copyNumPyToTeuchosArrayView(cols), copyNumPyToTeuchosArrayView(vals), atomic);
  }, "Sum into one or more sparse matrix entries, using local\n   row and column indices.\n\n For local row index  and local column indices\n cols, perform the update A(localRow, cols[k]) +=\n vals[k].  The row index and column indices must be valid\n on the calling process, and all matrix entries A(localRow,\n cols[k]) must already exist.  (This method does\n not change the matrix's structure.)  If the row index\n is valid, any invalid column indices are ignored, but counted\n in the return value.\n\n This overload of the method takes the column indices and\n values as Teuchos::ArrayView.  See above for an overload that\n takes Kokkos::View instead.\n\n \n [in] Local index of a row.  This row\n   must be owned by the calling process.\n \n\n [in] Local indices of the columns whose entries we\n   want to modify.\n \n\n [in] Values corresponding to the above column\n   indices.  vals[k] corresponds to cols[k].\n \n\n [in] Whether to use atomic updates.\n\n \n The number of indices for which values were actually\n   modified; the number of \"correct\" indices.\n\n This method has the same preconditions and return value\n meaning as replaceLocalValues() (which see).\n\nC++: Tpetra::CrsMatrix<double, int, long long, Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Serial, Kokkos::HostSpace> >::sumIntoLocalValues(const int, const class Teuchos::ArrayView<const int> &, const class Teuchos::ArrayView<const double> &, const bool) --> int", pybind11::arg("localRow"), pybind11::arg("indices"), pybind11::arg("values"), pybind11::arg("atomic"));
}

template<typename T>
void define_Vector_member_functions(T cl) {
  using SCALAR = typename T::type::scalar_type;
  using LO = typename T::type::local_ordinal_type;
  using GO = typename T::type::global_ordinal_type;
  using NODE = typename T::type::node_type;

  cl.def("getLocalViewHost",[](Teuchos::RCP<Tpetra::Vector<SCALAR,LO,GO,NODE>> &m){
      return getLocalViewHost(m);
  });
  cl.def("setLocalViewHost",[](Teuchos::RCP<Tpetra::Vector<SCALAR,LO,GO,NODE>> &m, py::array_t<SCALAR> input){
      return setLocalViewHost(m, input);
  });
}

template<typename T>
void define_MultiVector_member_functions(T cl) {
  using SCALAR = typename T::type::scalar_type;
  using LO = typename T::type::local_ordinal_type;
  using GO = typename T::type::global_ordinal_type;
  using NODE = typename T::type::node_type;

  cl.def("getLocalViewHost",[](Teuchos::RCP<Tpetra::MultiVector<SCALAR,LO,GO,NODE>> &m){
      return getLocalViewHost(m);
  });
  cl.def("setLocalViewHost",[](Teuchos::RCP<Tpetra::MultiVector<SCALAR,LO,GO,NODE>> &m, py::array_t<SCALAR> input){
      return setLocalViewHost(m, input);
  });
}

template <typename T>
void def_initialize_Kokkos(T m) {
  m.def("initialize_Kokkos",[](int num_threads,
                               int num_devices,
                               int device_id){
        if(!Kokkos::is_initialized()) {
          Kokkos::InitializationSettings args;
          args.set_num_threads(num_threads);
          args.set_num_devices(num_devices);
          args.set_device_id(device_id);
          Kokkos::initialize(args);
        }
      },
      py::arg("num_threads") = -1,
      py::arg("num_devices") = -1,
      py::arg("device_id") = -1
    );
  m.def("finalize_Kokkos",[](){
      if(Kokkos::is_initialized())
          Kokkos::finalize();
      }
    );
}

#endif // PYTRILINOS2_TPETRA_CUSTOM
