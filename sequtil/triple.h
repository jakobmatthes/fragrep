#ifndef __GLIBCPP_INTERNAL_TRIPLE_H
#define __GLIBCPP_INTERNAL_TRIPLE_H

namespace std
{

/// Triple implementation for a triple class in C++.  Axel Mosig,
/// University of Leipzig, Departement of Computer Science,
/// Bioinformatics Group, 2004.
/// This class has
/// been obtained by straightforward modifications of the Pair
/// implementation of the gcc c++ stl implementation, which is
/// copyrighted (C) 2001 by the Free Software Foundation, Inc.
/// For further documnentation of the pair implementation, see
/// stl_pair.h in the stl implementation of pair.
/// Although not part of the stl, the implementation is completely
/// canonical, so that I impudently put it in the std namespace.
///
/// triple holds three objects of arbitrary type.
template <class _T1, class _T2, class _T3>
struct triple {
  typedef _T1 first_type;    ///<  @c first_type is the first bound type
  typedef _T2 second_type;   ///<  @c second_type is the second bound type
  typedef _T3 third_type;    ///<  @c third_type is the third bound type

  _T1 first;                 ///< @c first is a copy of the first object
  _T2 second;                ///< @c second is a copy of the second object
  _T3 third;                ///< @c third is a copy of the third object
#ifdef _GLIBCPP_RESOLVE_LIB_DEFECTS
//265.  std::triple::triple() effects overly restrictive
  /** The default constructor creates @c first, @c second and @c third using their
   *  respective default constructors.  */
  triple() : first(), second(), third() {}
#else
  triple() : first(_T1()), second(_T2()), third (_T3()) {}
#endif
  /** Three objects may be passed to a @c triple constructor to be copied.  */
  triple(const _T1& __a, const _T2& __b, const _T3& __c) : first(__a), second(__b), third(__c) {}

  /** There is also a templated copy ctor for the @c triple class itself.  */
  template <class _U1, class _U2, class _U3>
  triple(const triple<_U1, _U2, _U3>& __p) : first(__p.first), second(__p.second), third(__p.third) {}
};

/// Two triples of the same type are equal iff their members are equal.
template <class _T1, class _T2, class _T3>
inline bool operator==(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y)
{ 
	return __x.first == __y.first && __x.second == __y.second && __x.third == __y.third;
}

/// <http://gcc.gnu.org/onlinedocs/libstdc++/20_util/howto.html#triplelt>
template <class _T1, class _T2, class _T3>
inline bool operator<(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y)
{ 
  return __x.first < __y.first || 
         (!(__y.first < __x.first) && __x.second < __y.second) ||
	  ((!(__y.first < __x.first) && !(__y.second < __x.second) && __x.third < __y.third)); 
}

/// Uses @c operator== to find the result.
template <class _T1, class _T2, class _T3>
inline bool operator!=(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y) {
  return !(__x == __y);
}

/// Uses @c operator< to find the result.
template <class _T1, class _T2, class _T3>
inline bool operator>(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y) {
  return __y < __x;
}

/// Uses @c operator< to find the result.
template <class _T1, class _T2, class _T3>
inline bool operator<=(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y) {
  return !(__y < __x);
}

/// Uses @c operator< to find the result.
template <class _T1, class _T2, class _T3>
inline bool operator>=(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y) {
  return !(__x < __y);
}

/**
 *  @brief A convenience wrapper for creating a triple from two objects.
 *  @param  x  The first object.
 *  @param  y  The second object.
 *  @param  z  The third object.
 *  @return   A newly-constructed triple<> object of the appropriate type.
 *
 *  The standard requires that the objects be passed by reference-to-const,
 *  but LWG issue #181 says they should be passed by const value.  We follow
 *  the LWG by default.
*/
template <class _T1, class _T2, class _T3>
#ifdef _GLIBCPP_RESOLVE_LIB_DEFECTS
//181.  make_triple() unintended behavior
inline triple<_T1, _T2, _T3> make_triple(_T1 __x, _T2 __y, _T3 __z)
#else
inline triple<_T1, _T2, _T3> make_triple(const _T1& __x, const _T2& __y, const _T3& __z)
#endif
{
  return triple<_T1, _T2, _T3>(__x, __y, __z);
}

} // namespace std

#endif /* __GLIBCPP_INTERNAL_TRIPLE_H */

// Local Variables:
// mode:C++
// End:
