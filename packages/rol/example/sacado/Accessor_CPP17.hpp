#ifndef ACCESSOR_HPP
#define ACCESSOR_HPP

#include <type_traits>
#include <utility>

// --- C++17 Trait Definitions ---

// 1. Trait to check if T is a specialization of Accessor (from previous answer)
template <typename T>
struct is_accessor_specialization_trait_cpp17 : std::false_type {};
template <typename VectorType>
struct is_accessor_specialization_trait_cpp17<Accessor<VectorType>> : std::true_type {};

// 2. Trait to check for member types
template <typename T, typename = void>
struct has_required_member_types : std::false_type {};

template <typename T>

                                        typename T::value_type,
                                        typename T::size_type
                                    >> : std::true_type {};

// 3. Trait to check for get_value method (if member types exist)
template <typename T, bool HasTypes, typename = void> // bool HasTypes is a precondition
struct has_correct_get_value_method : std::false_type {};

template <typename T>
struct has_correct_get_value_method<T, /* HasTypes = */ true,
    std::void_t<decltype(
        std::declval<const T&>().get_value(
            std::declval<const typename T::container_type&>(),
            std::declval<typename T::size_type>()
        )
    )>
> : std::is_same<
        decltype(
            std::declval<const T&>().get_value(
                std::declval<const typename T::container_type&>(),
                std::declval<typename T::size_type>()
            )
        ),
        typename T::value_type
    > {};

// 4. Trait to check for set_value method (if member types exist)
template <typename T, bool HasTypes, typename = void> // bool HasTypes is a precondition
struct has_correct_set_value_method : std::false_type {};

template <typename T>
struct has_correct_set_value_method<T, /* HasTypes = */ true,
    std::void_t<decltype(
        std::declval<const T&>().set_value(
            std::declval<typename T::container_type&>(), // Note: non-const container
            std::declval<typename T::size_type>(),
            std::declval<typename T::value_type>()
        )
    )>
> : std::is_same< // Check return type is void
        decltype(
            std::declval<const T&>().set_value(
                std::declval<typename T::container_type&>(),
                std::declval<typename T::size_type>(),
                std::declval<typename T::value_type>()
            )
        ),
        void
    > {};

// 5. Combine all checks for C++17
template<typename T>
struct is_compliant_accessor_trait {
private:
    static constexpr bool is_specialization = is_accessor_specialization_trait_cpp17<T>::value;
    static constexpr bool has_types = has_required_member_types<T>::value;
    // Only check methods if it's a specialization AND has the types, to avoid hard errors with declval on non-existent types
    static constexpr bool has_get = has_correct_get_value_method<T, is_specialization && has_types>::value;
    static constexpr bool has_set = has_correct_set_value_method<T, is_specialization && has_types>::value;

public:
    static constexpr bool value = is_specialization && has_types && has_get && has_set;
};


template <typename T>
inline constexpr bool is_compliant_accessor_v = is_compliant_accessor_trait<T>::value;

#endif // ACCESSOR_HPP
