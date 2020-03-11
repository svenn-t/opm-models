// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Properties
 * \ingroup TypeTraits
 * \author Timo Koch
 * \brief The Opm property system, traits with inheritance
 */
#ifndef OPM_PROPERTY_SYSTEM_HH
#define OPM_PROPERTY_SYSTEM_HH

#include <tuple>
#include <type_traits>
#include <iostream>

namespace Opm {
namespace Properties {

//! a tag to mark properties as undefined
struct UndefinedProperty {};

template <class TypeTag, class MyTypeTag>
struct Splices
{
     using tuple = std::tuple<>;
};
    

//! implementation details for template meta programming
namespace Detail {

//! check if a property P is defined
template<class P>
constexpr auto isDefinedProperty(int)
-> decltype(std::integral_constant<bool, !std::is_same<typename P::type, UndefinedProperty>::value>{})
{ return {}; }

//! fall back if a Property is defined
template<class P>
constexpr std::true_type isDefinedProperty(...) { return {}; }

//! check if a TypeTag inherits from other TypeTags
//! the enable_if portion of decltype is only needed for the macro hack to work, if no macros are in use anymore it can be removed,
//! i.e. then trailing return type is then -> decltype(std::declval<typename T::InheritsFrom>(), std::true_type{})
template<class T>
constexpr auto hasParentTypeTag(int)
-> decltype(std::declval<typename T::InheritsFrom>(), std::enable_if_t<!std::is_same<typename T::InheritsFrom, void>::value, int>{}, std::true_type{})
{ return {}; }

//! fall back if a TypeTag doesn't inherit
template<class T>
constexpr std::false_type hasParentTypeTag(...) { return {}; }

//! helper alias to concatenate multiple tuples
template<class ...Tuples>
using ConCatTuples = decltype(std::tuple_cat(std::declval<Tuples>()...));

//! helper struct to get the first property that is defined in the TypeTag hierarchy
template<class TypeTag, template<class,class> class Property, class TTagList>
struct GetDefined;

//! helper struct to iterate over the TypeTag hierarchy
template<class TypeTag, template<class,class> class Property, class TTagList, class Enable>
struct GetNextTypeTag;

template<class TypeTag, template<class,class> class Property, class LastTypeTag>
struct GetNextTypeTag<TypeTag, Property, std::tuple<LastTypeTag>, std::enable_if_t<hasParentTypeTag<LastTypeTag>(int{}), void>>
{ using type = typename GetDefined<TypeTag, Property, typename LastTypeTag::InheritsFrom>::type; };

template<class TypeTag, template<class,class> class Property, class LastTypeTag>
struct GetNextTypeTag<TypeTag, Property, std::tuple<LastTypeTag>, std::enable_if_t<!hasParentTypeTag<LastTypeTag>(int{}), void>>
{ using type = UndefinedProperty; };

template<class TypeTag, template<class,class> class Property, class FirstTypeTag, class ...Args>
struct GetNextTypeTag<TypeTag, Property, std::tuple<FirstTypeTag, Args...>, std::enable_if_t<hasParentTypeTag<FirstTypeTag>(int{}), void>>
{ using type = typename GetDefined<TypeTag, Property, ConCatTuples<typename FirstTypeTag::InheritsFrom, std::tuple<Args...>>>::type; };

template<class TypeTag, template<class,class> class Property, class FirstTypeTag, class ...Args>
struct GetNextTypeTag<TypeTag, Property, std::tuple<FirstTypeTag, Args...>, std::enable_if_t<!hasParentTypeTag<FirstTypeTag>(int{}), void>>
{ using type = typename GetDefined<TypeTag, Property, std::tuple<Args...>>::type; };

template<class TypeTag, template<class,class> class Property, class LastTypeTag>
struct GetDefined<TypeTag, Property, std::tuple<LastTypeTag>>
{
// For clang, the following alias triggers compiler warnings if instantiated
// from something like `GetPropType<..., DeprecatedProperty>`, even if that is
// contained in a diagnostic pragma construct that should prevent these warnings.
// As a workaround, also add the pragmas around this line.
// See the discussion in MR 1647 for more details.
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif
     using LastType = Property<TypeTag, LastTypeTag>;
#ifdef __clang__
#pragma clang diagnostic pop
#endif
     using type = std::conditional_t<isDefinedProperty<LastType>(int{}), LastType,
                                     typename GetNextTypeTag<TypeTag, Property, std::tuple<LastTypeTag>, void>::type>;
};

template<class TypeTag, template<class,class> class Property, class FirstTypeTag, class ...Args>
struct GetDefined<TypeTag, Property, std::tuple<FirstTypeTag, Args...>>
{
// See the comment above.
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif
     using FirstType = Property<TypeTag, FirstTypeTag>;
#ifdef __clang__
#pragma clang diagnostic pop
#endif
     using type = std::conditional_t<isDefinedProperty<FirstType>(int{}), FirstType,
                                     typename GetNextTypeTag<TypeTag, Property, std::tuple<FirstTypeTag, Args...>, void>::type>;
};


//! check if a splice S is defined
template<class S>
constexpr auto isDefinedSplice(int)
-> decltype(std::integral_constant<bool, !std::is_same<typename S::tuple, std::tuple<>>::value>{})
{ return {}; }

//! fall back if a splice is defined
template<class S>
constexpr std::true_type isDefinedSplice(...) { return {}; }

template<class TypeTag, class TTagList>
struct GetSplicesTypeTags;

template<class TypeTag, class TTagList, class Enable>
struct GetNextSplicesTypeTag;

template<class TypeTag, class LastTypeTag>
struct GetNextSplicesTypeTag<TypeTag, std::tuple<LastTypeTag>, std::enable_if_t<hasParentTypeTag<LastTypeTag>(int{}), void>>
{ using tuple = typename GetSplicesTypeTags<TypeTag, typename LastTypeTag::InheritsFrom>::tuple; };

template<class TypeTag, class LastTypeTag>
struct GetNextSplicesTypeTag<TypeTag, std::tuple<LastTypeTag>, std::enable_if_t<!hasParentTypeTag<LastTypeTag>(int{}), void>>
{ using tuple = std::tuple<>; };

template<class TypeTag, class FirstTypeTag, class ...Args>
struct GetNextSplicesTypeTag<TypeTag, std::tuple<FirstTypeTag, Args...>, std::enable_if_t<hasParentTypeTag<FirstTypeTag>(int{}), void>>
{ using tuple = typename GetSplicesTypeTags<TypeTag, ConCatTuples<typename FirstTypeTag::InheritsFrom, std::tuple<Args...>>>::tuple; };

template<class TypeTag, class FirstTypeTag, class ...Args>
struct GetNextSplicesTypeTag<TypeTag, std::tuple<FirstTypeTag, Args...>, std::enable_if_t<!hasParentTypeTag<FirstTypeTag>(int{}), void>>
{ using tuple = typename GetSplicesTypeTags<TypeTag, std::tuple<Args...>>::tuple; };

template<class TypeTag, class LastTypeTag>
struct GetSplicesTypeTags<TypeTag, std::tuple<LastTypeTag>>
{
     using LastSplices = Splices<TypeTag, LastTypeTag>;
     using nexttuple = typename GetNextSplicesTypeTag<TypeTag, std::tuple<LastTypeTag>, void>::tuple;
     // originally intended
//      using tuple = std::conditional_t<isDefinedSplice<LastSplices>(int{}),
//                                       typename ConCatTuples<nexttuple, typename LastSplices::tuple>::tuple,
//                                       nexttuple>;
     using tuple = std::conditional_t<isDefinedSplice<LastSplices>(int{}),
                                      typename LastSplices::tuple,
                                      typename GetNextSplicesTypeTag<TypeTag, std::tuple<LastTypeTag>, void>::tuple>;
};

template<class TypeTag, class FirstTypeTag, class ...Args>
struct GetSplicesTypeTags<TypeTag, std::tuple<FirstTypeTag, Args...>>
{
     using FirstSplices = Splices<TypeTag, FirstTypeTag>;
     using nexttuple = typename GetNextSplicesTypeTag<TypeTag, std::tuple<FirstTypeTag, Args...>, void>::tuple;
     // originally intended
//      using tuple = std::conditional_t<isDefinedSplice<FirstSplices>(int{}),
//                                       typename ConCatTuples<typename FirstSplices::tuple, nexttuple>::tuple,
//                                       nexttuple>;
     using tuple = std::conditional_t<isDefinedSplice<FirstSplices>(int{}),
                                      typename FirstSplices::tuple,
                                      typename GetNextSplicesTypeTag<TypeTag, std::tuple<FirstTypeTag, Args...>, void>::tuple>;
};

} // end namespace Detail

template<class TypeTag, class MyTypeTag>
struct SpatialDiscretizationSplice { using type = UndefinedProperty; };

namespace TTag {
struct MyTTag {};
struct MySDTTag {};
}

template<class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::MyTTag>
{ using type = TTag::MySDTTag; };

template<class TypeTag, class MyTypeTag>
struct VtkOutputFormat { using type = UndefinedProperty; };

template<class TypeTag>
struct VtkOutputFormat<TypeTag, TTag::MySDTTag>
{ static constexpr auto value = 2; };

namespace TTag {
struct FvBaseDiscretization;
struct VcfvDiscretization;
struct MultiPhaseBaseModel;
}

namespace Detail {

    //! helper struct to extract get the Property specilization given a TypeTag, asserts that the property is defined
template<class TypeTag, template<class,class> class Property>
struct GetPropImpl
{
    using PType = Properties::SpatialDiscretizationSplice<TypeTag, TTag::MyTTag>;
    using testtype = std::conditional_t<isDefinedProperty<PType>(int{}), PType, UndefinedProperty>;
    static_assert(!std::is_same<testtype, UndefinedProperty>::value, "SpatialDiscretizationSplice is undefined in MyTTag!");

    using QType = Properties::VtkOutputFormat<TypeTag, TTag::MySDTTag>;
    using testtype2 = std::conditional_t<isDefinedProperty<QType>(int{}), QType, UndefinedProperty>;
    static_assert(!std::is_same<testtype2, UndefinedProperty>::value, "VtkOutputFormat is undefined in MySDTTag!");

    using RType = Properties::VtkOutputFormat<TypeTag, TTag::FvBaseDiscretization>;
    using testtype3 = std::conditional_t<isDefinedProperty<RType>(int{}), RType, UndefinedProperty>;
    static_assert(!std::is_same<testtype3, UndefinedProperty>::value, "VtkOutputFormat is undefined in FvBaseDiscretization!");

    // works:
//     using type = typename Detail::GetDefined<TypeTag,
//                                              Property,
//                                              std::tuple<TypeTag, TTag::VcfvDiscretization>
//                                             >::type;
    using tuple = typename GetSplicesTypeTags<TypeTag, std::tuple<TypeTag>>::tuple;
    using type = typename Detail::GetDefined<TypeTag,
                                             Property,
                                             typename ConCatTuples<std::tuple<TypeTag>, tuple>::type
                                            >::type;
    static_assert(!std::is_same<type, UndefinedProperty>::value, "Property is undefined!");
};

} // end namespace Detail
} // end namespace Property

//! get the type of a property (equivalent to old macro GET_PROP(...))
template<class TypeTag, template<class,class> class Property>
using GetProp = typename Properties::Detail::GetPropImpl<TypeTag, Property>::type;

// See the comment above.
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif
//! get the type alias defined in the property (equivalent to old macro GET_PROP_TYPE(...))
template<class TypeTag, template<class,class> class Property>
using GetPropType = typename Properties::Detail::GetPropImpl<TypeTag, Property>::type::type;

//! get the value data member of a property
template<class TypeTag, template<class,class> class Property>
constexpr auto getPropValue() { return Properties::Detail::GetPropImpl<TypeTag, Property>::type::value; }
#ifdef __clang__
#pragma clang diagnostic pop
#endif

namespace Properties {
template <class TypeTag>
void printValues(std::ostream& os = std::cout)
{
    os <<
    "The eWoms property system was compiled with the macro\n"
    "NO_PROPERTY_INTROSPECTION defined.\n"
    "No diagnostic messages this time, sorry.\n";
}
}

/*!
 * \ingroup Properties
 * \brief Indicates that property definitions follow
 */
#define BEGIN_PROPERTIES namespace Opm { namespace Properties {

/*!
 * \ingroup Properties
 * \brief Indicates that all properties have been specified (for now)
 */
#define END_PROPERTIES }}

} // end namespace Opm

#endif
