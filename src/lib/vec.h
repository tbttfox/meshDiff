#pragma once

#include "vecext.h"

const vec4f zero4f = makevec4f(0, 0, 0, 0);
const vec2f x2f = makevec2f(1, 0);
const vec2f y2f = makevec2f(0, 1);

inline vec2i makevec2i(int *v) { return makevec2<int>(v[0], v[1]); }
inline vec2f makevec2f(float *v) { return makevec2<float>(v[0], v[1]); }

inline vec3f makevec3f(float *v) { return makevec3<float>(v[0], v[1], v[2]); }
inline vec3i makevec3i(int *v) { return makevec3<int>(v[0], v[1], v[2]); }

inline vec4f makevec4f(float *v) { return makevec4<float>(v[0], v[1], v[2], v[3]); }
inline vec4i makevec4i(int *v) { return makevec4<int>(v[0], v[1], v[2], v[3]); }

inline vec3i floor(vec3f v) { return makevec3i((int)floorf(v.x), (int)floorf(v.y), (int)floorf(v.z)); }

inline vec3f add(const vec3f &v0, const vec3i &v1) { return makevec3f(v0.x + (float)v1.x, v0.y + (float)v1.y, v0.z + (float)v1.z); }
inline vec3f add(const vec3i &v0, const vec3f &v1) { return makevec3f((float)v0.x + v1.x, (float)v0.y + v1.y, (float)v0.z + v1.z); }

inline bool operator==(const vec3i &lhs, const vec3i &rhs) { return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z; }

