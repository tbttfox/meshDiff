#pragma once

#include "vec.h"

template <typename R>
struct ray3 {
	static R epsilon;
	static R infinity;

	vec3<R> o; // origin
	vec3<R> d; // direction
	R mint, maxt; // bounds
};

template<typename R> R ray3<R>::epsilon = R(1e-5);
template<typename R> R ray3<R>::infinity = numeric_limits<R>::infinity();

template<typename R>
inline ray3<R> makeray3(const vec3<R>& o, const vec3<R>& d,
	R mint = ray3<R>::epsilon, R maxt = ray3<R>::infinity) {
	ray3<R> r; r.o = o; r.d = d; r.mint = mint; r.maxt = maxt; return r;
}

template<typename R>
inline ray3<R> makeraysegment3(const vec3<R>& p1, const vec3<R>& p2,
	R mint = ray3<R>::epsilon) {
	return makeray3<R>(p1, normalize(p2 - p1), mint, length(p2 - p1) - mint);
}

template<typename R>
inline vec3<R> evalray3(const ray3<R>& r, R t) { return r.o + r.d*t; }

typedef ray3<float> ray3f;
typedef ray3<double> ray3d;

inline ray3f makeray3f(const vec3f& o, const vec3f& d,
	float mint = ray3f::epsilon, float maxt = ray3f::infinity) {
	return makeray3<float>(o, d, mint, maxt);
}

inline ray3f makeraysegment3f(const vec3f& o, const vec3f& d, float mint = ray3f::epsilon) {
	return makeraysegment3<float>(o, d, mint);
}

