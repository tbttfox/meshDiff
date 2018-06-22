#pragma once

#include "math.h"
#include "vec.h"

template <typename R>
struct range2 {
	vec2<R> min, max;
};

template <typename R>
struct range3 {
	vec3<R> min, max; // min and max corners
};

template <typename R>
inline bool isvalid(const range3<R>& r) { return r.max.x >= r.min.x && r.max.y >= r.min.y && r.max.z >= r.min.z; }

template <typename R>
inline range2<R> makerange2(const vec2<R>& min, const vec2<R>& max) {
	range2<R> r; r.min = min; r.max = max; return r;
}
template <typename R>
inline range3<R> makerange3(const vec3<R>& min, const vec3<R>& max) {
	range3<R> r; r.min = min; r.max = max; return r;
}
template <typename R>
inline range2<R> invalidrange2() {
	return makerange2<R>(
		makevec2<R>(std::numeric_limits<R>::max(),
			std::numeric_limits<R>::max()),
		makevec2<R>(-std::numeric_limits<R>::max(),
			-std::numeric_limits<R>::max()));
}
template <typename R>
inline range3<R> invalidrange3() {
	return makerange3<R>(
		makevec3<R>(std::numeric_limits<R>::max(),
			std::numeric_limits<R>::max(),
			std::numeric_limits<R>::max()),
		makevec3<R>(-std::numeric_limits<R>::max(),
			-std::numeric_limits<R>::max(),
			-std::numeric_limits<R>::max()));
}

template<typename R>
inline range2<R> runion(const range2<R>& range, const vec2<R> &v) {
	return makerange2<R>(min_component(range.min, v), max_component(range.max, v));
}

template<typename R>
inline range2<R> runion(const range2<R>& range, const range2<R>& r) {
	return makerange2<R>(min_component(range.min, r.min), max_component(range.max, r.max));
}

template<typename R>
inline range3<R> runion(const range3<R>& range, const vec3<R> &v) {
	return makerange3<R>(min_component(range.min, v), max_component(range.max, v));
}

template<typename R>
inline range3<R> runion(const range3<R>& range, const range3<R>& r) {
	return makerange3<R>(min_component(range.min, r.min), max_component(range.max, r.max));
}

template<typename R>
inline vec3<R> center(const range3<R>& range) {
	error_if_not(isvalid(range), "center: invalid range");
	return (range.min + range.max) * 0.5f;
}

template<typename R>
inline vec3<R> extent(const range3<R>& range) {
	error_if_not(isvalid(range), "extent: invalid range");
	return (range.max - range.min);
}

//template <typename R>
//inline R diagonal(const range3<R>& range) {
//    return length(extent(range));
//}

//template <typename R>
//inline R diagonalSqr(const range3<R>& range) {
//    return lengthSqr(extent(range));
//}

template<typename R>
inline range3<R> scale(const range3<R>& range, R s) {
	error(isvalid(range), "scale: invalid range");
	vec3<R> rc = center(range);
	vec3<R> rs = size(range);
	return makerange3<R>(rc - rs * ((R)0.5*s), rc + rs * ((R)0.5*s));
}

typedef range2<int>     range2i;
typedef range2<float>   range2f;
typedef range2<double>  range2d;

typedef range3<int>     range3i;
typedef range3<float>   range3f;
typedef range3<double>  range3d;

inline range2i makerange2i(const vec2i& min, const vec2i& max) { return makerange2<int>(min, max); }
inline range2f makerange2f(const vec2f& min, const vec2f& max) { return makerange2<float>(min, max); }
inline range3i makerange3i(const vec3i& min, const vec3i& max) { return makerange3<int>(min, max); }
inline range3f makerange3f(const vec3f& min, const vec3f& max) { return makerange3<float>(min, max); }
const range3f invalidrange3f = invalidrange3<float>();
const range2f invalidrange2f = invalidrange2<float>();
const range3i invalidrange3i = invalidrange3<int>();
const range2i invalidrange2i = invalidrange2<int>();

