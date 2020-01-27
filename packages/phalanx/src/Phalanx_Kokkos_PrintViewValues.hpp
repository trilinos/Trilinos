#ifndef PHALANX_KOKKOS_PRINT_VIEW_VALUES_HPP
#define PHALANX_KOKKOS_PRINT_VIEW_VALUES_HPP

namespace PHX {

  template<typename ViewType,int Rank> struct PrintViewValues;

  // Dynamic rank view specialization
  template<typename ViewType>
  struct PrintViewValues<ViewType,0>  {
    void print(const ViewType& vd, std::ostream& os)
    {
      typename ViewType::HostMirror v = Kokkos::create_mirror(vd);
      Kokkos::deep_copy(v,vd);
      const auto rank = v.rank();
      if (rank == 1) {
        for (int i=0; i < v.extent(0); ++i)
          os << i << " " << v(i) << std::endl;
      }
      else if (rank == 2) {
        for (int i=0; i < v.extent(0); ++i)
          for (int j=0; j < v.extent(1); ++j)
            os << i << " " << j << " " << v(i,j) << std::endl;
      }
      else if (rank == 3) {
        for (int i=0; i < v.extent(0); ++i)
          for (int j=0; j < v.extent(1); ++j)
            for (int k=0; k < v.extent(2); ++k)
            os << i << " " << j << " " << k << " " << v(i,j,k) << std::endl;
      }
      else if (rank == 4) {
        for (int i=0; i < v.extent(0); ++i)
          for (int j=0; j < v.extent(1); ++j)
            for (int k=0; k < v.extent(2); ++k)
              for (int m=0; m < v.extent(3); ++m)
                os << i << " " << j << " " << k << " " << m << " " 
                   << v(i,j,k,m) << std::endl;
      }
      else if (rank == 5) {
        for (int i=0; i < v.extent(0); ++i)
          for (int j=0; j < v.extent(1); ++j)
            for (int k=0; k < v.extent(2); ++k)
              for (int m=0; m < v.extent(3); ++m)
                for (int n=0; n < v.extent(4); ++n)
                  os << i << " " << j << " " << k << " " << m << " " 
                     << n << " " 
                     << v(i,j,k,m,n) << std::endl;
      }
      else if (rank == 6) {
        for (int i=0; i < v.extent(0); ++i)
          for (int j=0; j < v.extent(1); ++j)
            for (int k=0; k < v.extent(2); ++k)
              for (int m=0; m < v.extent(3); ++m)
                for (int n=0; n < v.extent(4); ++n)
                for (int p=0; p < v.extent(5); ++p)
                  os << i << " " << j << " " << k << " " << m << " " 
                     << n << " " << p << " "
                     << v(i,j,k,m,n,p) << std::endl;
      }
      else if (rank == 7) {
        for (int i=0; i < v.extent(0); ++i)
          for (int j=0; j < v.extent(1); ++j)
            for (int k=0; k < v.extent(2); ++k)
              for (int m=0; m < v.extent(3); ++m)
                for (int n=0; n < v.extent(4); ++n)
                  for (int p=0; p < v.extent(5); ++p)
                    for (int q=0; q < v.extent(6); ++q)
                      os << i << " " << j << " " << k << " " << m << " " 
                         << n << " " << p << " " << q << " "
                         << v(i,j,k,m,n,p,q) << std::endl;
      }
    }
  };

  template<typename ViewType>
  struct PrintViewValues<ViewType,1>  {
    void print(const ViewType& vd, std::ostream& os)
    {
      typename ViewType::HostMirror v = Kokkos::create_mirror(vd);
      Kokkos::deep_copy(v,vd);
      for (int i=0; i < v.extent(0); ++i)
        os << i << " " << v(i) << std::endl;
    }
  };

  template<typename ViewType>
  struct PrintViewValues<ViewType,2>  {
    void print(const ViewType& vd, std::ostream& os)
    {
      typename ViewType::HostMirror v = Kokkos::create_mirror(vd);
      Kokkos::deep_copy(v,vd);
      for (int i=0; i < v.extent(0); ++i)
        for (int j=0; j < v.extent(1); ++j)
          os << i << " " << j << " " << v(i,j) << std::endl;
    }
  };

  template<typename ViewType>
  struct PrintViewValues<ViewType,3>  {
    void print(const ViewType& vd, std::ostream& os)
    {
      typename ViewType::HostMirror v = Kokkos::create_mirror(vd);
      Kokkos::deep_copy(v,vd);
      for (int i=0; i < v.extent(0); ++i)
        for (int j=0; j < v.extent(1); ++j)
          for (int k=0; k < v.extent(2); ++k)
            os << i << " " << j << " " << k << " " << v(i,j,k) << std::endl;
    }
  };

  template<typename ViewType>
  struct PrintViewValues<ViewType,4>  {
    void print(const ViewType& vd, std::ostream& os)
    {
      typename ViewType::HostMirror v = Kokkos::create_mirror(vd);
      Kokkos::deep_copy(v,vd);
      for (int i=0; i < v.extent(0); ++i)
        for (int j=0; j < v.extent(1); ++j)
          for (int k=0; k < v.extent(2); ++k)
            for (int m=0; m < v.extent(3); ++m)
              os << i << " " << j << " " << k << " " << m << " " 
                 << v(i,j,k,m) << std::endl;
    }
  };

  template<typename ViewType>
  struct PrintViewValues<ViewType,5>  {
    void print(const ViewType& vd, std::ostream& os)
    {
      typename ViewType::HostMirror v = Kokkos::create_mirror(vd);
      Kokkos::deep_copy(v,vd);
      for (int i=0; i < v.extent(0); ++i)
        for (int j=0; j < v.extent(1); ++j)
          for (int k=0; k < v.extent(2); ++k)
            for (int m=0; m < v.extent(3); ++m)
              for (int n=0; n < v.extent(4); ++n)
                os << i << " " << j << " " << k << " " << m << " " 
                   << n << " " 
                   << v(i,j,k,m,n) << std::endl;
    }
  };

  template<typename ViewType>
  struct PrintViewValues<ViewType,6>  {
    void print(const ViewType& vd, std::ostream& os)
    {
      typename ViewType::HostMirror v = Kokkos::create_mirror(vd);
      Kokkos::deep_copy(v,vd);
      for (int i=0; i < v.extent(0); ++i)
        for (int j=0; j < v.extent(1); ++j)
          for (int k=0; k < v.extent(2); ++k)
            for (int m=0; m < v.extent(3); ++m)
              for (int n=0; n < v.extent(4); ++n)
                for (int p=0; p < v.extent(5); ++p)
                  os << i << " " << j << " " << k << " " << m << " " 
                     << n << " " << p << " "
                     << v(i,j,k,m,n,p) << std::endl;
    }
  };

  template<typename ViewType>
  struct PrintViewValues<ViewType,7>  {
    void print(const ViewType& vd, std::ostream& os)
    {
      typename ViewType::HostMirror v = Kokkos::create_mirror(vd);
      Kokkos::deep_copy(v,vd);
      for (int i=0; i < v.extent(0); ++i)
        for (int j=0; j < v.extent(1); ++j)
          for (int k=0; k < v.extent(2); ++k)
            for (int m=0; m < v.extent(3); ++m)
              for (int n=0; n < v.extent(4); ++n)
                for (int p=0; p < v.extent(5); ++p)
                  for (int q=0; q < v.extent(6); ++q)
                    os << i << " " << j << " " << k << " " << m << " " 
                       << n << " " << p << " " << q << " "
                       << v(i,j,k,m,n,p,q) << std::endl;
    }
  };

}

#endif
