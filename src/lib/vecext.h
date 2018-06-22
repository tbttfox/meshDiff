#pragma once

#include "math.h"

template <typename R>
struct vec2 {
	R x, y;
};

template <typename R>
struct vec3 {
	R x, y, z;

	inline R operator[](const int i) const {
		return (i == 0) ? x : ((i == 1) ? y : z);
	}
};

template <typename R>
struct vec4 {
	R x, y, z, w;
};

template<typename R>
inline vec2<R> makevec2(const R& x, const R& y) {
	vec2<R> v; v.x = x; v.y = y; return v;
}

template <typename R>
inline vec2<R> operator+(const vec2<R>& a, const vec2<R>& b) {
	return makevec2<R>(a.x + b.x, a.y + b.y);
}

template <typename R>
inline vec2<R> operator-(const vec2<R>& a, const vec2<R>& b) {
	return makevec2<R>(a.x - b.x, a.y - b.y);
}

template <typename R>
inline vec2<R> operator*(const vec2<R>& a, const R& b) {
	return makevec2<R>(a.x*b, a.y*b);
}

template <typename R>
inline vec2<R> operator*(const vec2<R>& a, const vec2<R>& b) {
	return makevec2<R>(a.x*b.x, a.y*b.y);
}

template <typename R>
inline vec2<R> operator/(const vec2<R>& a, const R& b) {
	return makevec2<R>(a.x / b, a.y / b);
}

template <typename R>
inline vec2<R> operator/(const R& a, const vec2<R>& b) {
	return makevec2<R>(a / b.x, a / b.y);
}

template <typename R>
inline vec2<R> operator/(const vec2<R>& a, const vec2<R>& b) {
	return makevec2<R>(a.x / b.x, a.y / b.y);
}

template<typename R>
inline vec3<R> makevec3(const R& x, const R& y, const R& z) {
	vec3<R> v; v.x = x; v.y = y; v.z = z; return v;
}

template <typename R>
inline vec3<R> operator-(const vec3<R>& a) {
	return makevec3<R>(-a.x, -a.y, -a.z);
}

template <typename R>
inline vec3<R> operator+(const vec3<R>& a, const vec3<R>& b) {
	return makevec3<R>(a.x + b.x, a.y + b.y, a.z + b.z);
}

template <typename R>
inline vec3<R> operator-(const vec3<R>& a, const vec3<R>& b) {
	return makevec3<R>(a.x - b.x, a.y - b.y, a.z - b.z);
}

template <typename R>
inline vec3<R> operator*(const vec3<R>& a, const vec3<R>& b) {
	return makevec3<R>(a.x*b.x, a.y*b.y, a.z*b.z);
}

template <typename R>
inline vec3<R> operator*(const vec3<R>& a, const R& b) {
	return makevec3<R>(a.x*b, a.y*b, a.z*b);
}

template <typename R>
inline vec3<R> operator*(const R& b, const vec3<R>& a) {
	return makevec3<R>(a.x*b, a.y*b, a.z*b);
}


template <typename R>
inline vec3<R> operator/(const vec3<R>& a, const vec3<R>& b) {
	return makevec3<R>(a.x / b.x, a.y / b.y, a.z / b.z);
}

template <typename R>
inline vec3<R> operator/(const vec3<R>& a, const R& b) {
	return makevec3<R>(a.x / b, a.y / b, a.z / b);
}

template <typename R>
inline vec3<R> operator/(const R& a, const vec3<R>& b) {
	return makevec3<R>(a / b.x, a / b.y, a / b.z);
}

template <typename R>
inline vec3<R>& operator+=(vec3<R>& a, const vec3<R>& b) {
	a.x += b.x; a.y += b.y; a.z += b.z;
	return a;
}

template <typename R>
inline vec3<R> operator-=(vec3<R>& a, const vec3<R>& b) {
	a.x -= b.x; a.y -= b.y; a.z -= b.z;
	return a;
}

template <typename R>
inline vec3<R>& operator*=(vec3<R>& a, const vec3<R>& b) {
	a.x *= b.x; a.y *= b.y; a.z *= b.z;
	return a;
}

template <typename R>
inline vec3<R>& operator*=(vec3<R>& a, const R& b) {
	a.x *= b; a.y *= b; a.z *= b;
	return a;
}

template <typename R>
inline vec3<R>& operator/=(vec3<R>& a, const vec3<R>& b) {
	a.x /= b.x; a.y /= b.y; a.z /= b.z;
	return a;
}

template <typename R>
inline vec3<R>& operator/=(vec3<R>& a, const R& b) {
	a.x /= b; a.y /= b; a.z /= b;
	return a;
}

template <typename R>
inline R operator%(const vec3<R>& a, const vec3<R>& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

template <typename R>
inline R dot(const vec3<R>& a, const vec3<R>& b) {
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

template <typename R>
inline vec3<R> operator^(const vec3<R>& a, const vec3<R>& b) {
	return makevec3<R>(a.y*b.z - a.z*b.y, -a.x*b.z + a.z*b.x, a.x*b.y - a.y*b.x);
}

template <typename R>
inline vec3<R> cross(const vec3<R>& a, const vec3<R>& b) {
	return makevec3<R>(a.y*b.z - a.z*b.y, -a.x*b.z + a.z*b.x, a.x*b.y - a.y*b.x);
}

template <typename R>
inline R length(const vec3<R>& a) {
	return sqrt(a % a);
}

template <typename R>
inline vec3<R> abs_component(const vec3<R>& a) {
	return makevec3(abs(a.x), abs(a.y), abs(a.z));
}

template <typename R>
inline vec3<R> exp_component(const vec3<R>& a) {
	return makevec3(exp(a.x), exp(a.y), exp(a.z));
}

template <typename R>
inline vec3<R> pow_component(const vec3<R>& a, R e) {
	return makevec3(pow(a.x, e), pow(a.y, e), pow(a.z, e));
}

template <typename R>
inline vec3<R> sqrt_component(const vec3<R>& a) {
	return makevec3((R)sqrt(a.x), (R)sqrt(a.y), (R)sqrt(a.z));
}

template <typename R>
inline R lengthSqr(const vec3<R>& a) {
	return a % a;
}

template <typename R>
inline vec3<R> normalize(const vec3<R>& a) {
	R l = sqrt(a % a);
	if (l == 0) return a;
	else return a * (1 / l);
}

template <typename R>
inline R average(const vec3<R>& a) {
	return (a.x + a.y + a.z) / 3;
}

template <typename R>
inline R luminance(const vec3<R>& a) {
	return (a.x * 0.299 + a.y * 0.587 + a.z * 0.114);
}

template <typename R>
inline bool iszero(const vec3<R>& a) {
	return a.x == 0 && a.y == 0 && a.z == 0;
}

template <typename R>
inline bool hasnan(const vec3<R>& a) {
	return isnan(a.x) || isnan(a.y) || isnan(a.z);
}

template <typename R>
inline bool hasinf(const vec3<R>& a) {
	return isinf(a.x) || isinf(a.y) || isinf(a.z);
}

template <typename R>
inline bool hasneg(const vec3<R>& a) {
	return a.x < 0 || a.y < 0 || a.z < 0;
}

template <typename R>
inline bool hasinvalid(const vec3<R>& a) {
	return hasinf(a) || hasnan(a) || hasneg(a);
}

template <typename R>
inline vec2<R> min_component(const vec2<R>& a, const vec2<R>& b) {
	return makevec2<R>(std::min(a.x, b.x), std::min(a.y, b.y));
}

template <typename R>
inline vec2<R> max_component(const vec2<R>& a, const vec2<R>& b) {
	return makevec2<R>(std::max(a.x, b.x), std::max(a.y, b.y));
}

template <typename R>
inline vec3<R> min_component(const vec3<R>& a, const vec3<R>& b) {
	return makevec3<R>(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
}

template <typename R>
inline vec3<R> max_component(const vec3<R>& a, const vec3<R>& b) {
	return makevec3<R>(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
}

template <typename R>
inline int max_component_index(const vec3<R>& v) {
	if (v.x >= v.y && v.x >= v.z) return 0;
	else if (v.y >= v.x && v.y >= v.z) return 1;
	else return 2;
}

template <typename R>
inline R get_component(const vec3<R>& v, int index) {
	return ((R*)&v.x)[index];
}

template <typename R>
inline R& get_component(vec3<R>& v, int index) {
	return ((R*)&v.x)[index];
}

template<typename R>
inline vec4<R> makevec4(const R& x, const R& y, const R& z, const R& w) {
	vec4<R> v; v.x = x; v.y = y; v.z = z; v.w = w; return v;
}

typedef vec2<float> vec2f;
typedef vec2<double> vec2d;
typedef vec2<int> vec2i;

typedef vec3<float> vec3f;
typedef vec3<double> vec3d;
typedef vec3<int> vec3i;

typedef vec4<float> vec4f;
typedef vec4<double> vec4d;
typedef vec4<int> vec4i;

inline vec2f makevec2f(float x, float y) { return makevec2<float>(x, y); }
inline vec2i makevec2i(int x, int y) { return makevec2<int>(x, y); }
inline vec2d makevec2d(double x, double y) { return makevec2<double>(x, y); }

inline vec3f makevec3f(float x, float y, float z) { return makevec3<float>(x, y, z); }
inline vec3i makevec3i(int x, int y, int z) { return makevec3<int>(x, y, z); }
inline vec3d makevec3d(double x, double y, double z) { return makevec3<double>(x, y, z); }

inline vec4f makevec4f(float x, float y, float z, float w) { return makevec4<float>(x, y, z, w); }
inline vec4i makevec4i(int x, int y, int z, int w) { return makevec4<int>(x, y, z, w); }
inline vec4d makevec4d(double x, double y, double z, double w) { return makevec4<double>(x, y, z, w); }

const vec3f zero3f = makevec3f(0.0f, 0.0f, 0.0f);
const vec3d zero3d = makevec3d(0.0, 0.0, 0.0);
const vec3i zero3i = makevec3i(0, 0, 0);

const vec2f zero2f = makevec2f(0.0f, 0.0f);
const vec2d zero2d = makevec2d(0.0, 0.0);
const vec2i zero2i = makevec2i(0, 0);

const vec3f half3f = makevec3f(0.5f, 0.5f, 0.5f);
const vec3d half3d = makevec3d(0.5, 0.5, 0.5);

const vec2f half2f = makevec2f(0.5f, 0.5f);
const vec2d half2d = makevec2d(0.5, 0.5);

const vec3f one3f = makevec3f(1.0f, 1.0f, 1.0f);
const vec3d one3d = makevec3d(1.0, 1.0, 1.0);
const vec3i one3i = makevec3i(1, 1, 1);

const vec2f one2f = makevec2f(1.0f, 1.0f);
const vec2d one2d = makevec2d(1.0, 1.0);
const vec2i one2i = makevec2i(1, 1);

const vec3f x3f = makevec3f(1, 0, 0);
const vec3f y3f = makevec3f(0, 1, 0);
const vec3f z3f = makevec3f(0, 0, 1);

const vec3d x3d = makevec3d(1, 0, 0);
const vec3d y3d = makevec3d(0, 1, 0);
const vec3d z3d = makevec3d(0, 0, 1);

