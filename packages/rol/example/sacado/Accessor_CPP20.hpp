#include <concepts>
#include <utility>

// Helper to check if T is a specialization of Accessor (from previous answer)
template <typename T>
struct is_accessor_specialization_trait_cpp20 : std::false_type {};

template <typename VectorType>
struct is_accessor_specialization_trait_cpp20<Accessor<VectorType>> : std::true_type {};

template <typename T>
concept IsAccessorSpecialization = is_accessor_specialization_trait_cpp20<T>::value;

// Concept to check for full compliance
template<typename Acc>
concept CompliantAccessor = IsAccessorSpecialization<Acc> && requires (
    const Acc& accessor,
    // Use placeholders that will be deduced based on Acc's member types
    // If Acc::container_type doesn't exist, this part of the requires clause will fail.
    typename Acc::container_type& mutable_container,
    const typename Acc::container_type& const_container,
    typename Acc::size_type index,
    typename Acc::value_type val
) {
    // 1. Check for required member types
    typename Acc::container_type;
    typename Acc::value_type;
    typename Acc::size_type;

    // 2. Check for get_value member function
    //    - const-correct
    //    - takes const container_type&, size_type
    //    - returns value_type
    { accessor.get_value(const_container, index) } -> std::same_as<typename Acc::value_type>;

    // 3. Check for set_value member function
    //    - const-correct (the method itself is const)
    //    - takes container_type&, size_type, value_type
    //    - returns void
    { accessor.set_value(mutable_container, index, val) } -> std::same_as<void>;
};
