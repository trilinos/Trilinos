namespace Kokkos {
namespace Example {

//Get mean values matrix for mean-based preconditioning
template <class ViewType>
class GetMeanValsFunc {

  public:

  typedef typename ViewType::device_type device_type;
  GetMeanValsFunc(ViewType mean_vals, ViewType vals) : mean_vals_(mean_vals), vals_(vals) {}

  private:
  ViewType mean_vals_;
  ViewType vals_;

  public:

  KOKKOS_INLINE_FUNCTION
  void operator() (const int i) const {
    mean_vals_[i] = vals_[i].fastAccessCoeff(0);
  }

};


} // namespace Example
} // namespace Kokkos

