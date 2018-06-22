#pragma once

#include <map>
#include "vec.h"
#include "geom.h"

struct HalfEdgeAdjacency {
	struct HEdge;
	struct Face {
		int idx = -1;
		HEdge* hedge = 0;
	};
	struct Edge {
		int idx = -1;
		HEdge* hedge = 0;
		bool is_boundary() { return !hedge->face || !hedge->opp->face; }
	};
	struct Vert {
		int idx = -1;
		HEdge* hedge = 0;
	};
	struct HEdge {
		int idx = -1;
		HEdge* next = nullptr;
		HEdge* prev = 0;
		HEdge* opp = 0;
		Face* face = 0;
		Edge* edge = 0;
		Vert* vert = 0; // start
		bool is_boundary() { return !face; }
	};

	std::vector<Face*>   faces;
	std::vector<Vert*>   verts;
	std::vector<Edge*>   edges;
	std::vector<HEdge*>  hedges;

	template<typename T>
	HalfEdgeAdjacency(const std::vector<T>& face, int nvert) { from_faces(face, nvert); }

	~HalfEdgeAdjacency() { _clear(); }

	HEdge* halfedge(int v0idx, int v1idx) {
		Vert* v0 = verts[v0idx];
		Vert* v1 = verts[v1idx];
		HEdge* he = v0->hedge;
		HEdge* start = v0->hedge;
		do {
			if (he->opp->vert == v1) return he;
			he = he->prev->opp;
		} while (he != start);
		return 0;
	}

	template<typename T>
	void from_faces(const std::vector<T>& face, int nvert) {
		// add face
		for (auto f : face) add_face();
		// add vert
		for (int i = 0; i < nvert; i++) add_vert();
		// add edge and hedge
		map<pair<int, int>, HEdge*> edgemap;
		for (auto f : faces) {
			auto ff = face[f->idx];
			for (int i = 0; i < ff.size(); i++) {
				auto v0 = ff[i]; auto v1 = ff[(i + 1) % ff.size()];
				auto v01 = make_pair(v0, v1);
				auto v10 = make_pair(v1, v0);
				if (edgemap.find(v01) == edgemap.end()) {
					auto he01 = add_hedge();
					auto he10 = add_hedge();
					auto e = add_edge();
					he01->vert = verts[v0];
					he10->vert = verts[v1];
					he01->opp = he10;
					he10->opp = he01;
					e->hedge = he01;
					he01->edge = e;
					he10->edge = e;
					edgemap[v01] = he01;
					edgemap[v10] = he10;
					verts[v0]->hedge = he01;
					verts[v1]->hedge = he10;
				}
				auto he01 = edgemap[v01];
				he01->face = f;
				if (i == 0) f->hedge = he01;
			}
		}
		// connect hedge on face
		for (auto f : faces) {
			auto ff = face[f->idx];
			for (int i = 0; i < ff.size(); i++) {
				auto v0 = ff[i]; auto v1 = ff[(i + 1) % ff.size()]; auto v2 = ff[(i + 2) % ff.size()];
				auto v01 = make_pair(v0, v1);
				auto v12 = make_pair(v1, v2);
				edgemap[v01]->next = edgemap[v12];
			}
		}
		// connect hedge on boundary
		for (auto he : hedges) {
			if (not he->is_boundary()) continue;
			auto v = he->opp->vert;
			HEdge* next = 0;
			for (auto hen : hedges) if (hen->is_boundary() and hen->vert == v) { next = hen; break; }
			he->next = next;
		}
		// connect prev hedge
		for (auto he : hedges) he->next->prev = he;
	}

	void _clear() {
		for (auto v : faces) delete v;
		for (auto v : edges) delete v;
		for (auto v : verts) delete v;
		for (auto v : hedges) delete v;
	}

	template<typename T>
	T* _add(std::vector<T*>& ptr) {
		auto value = new T();
		value->idx = ptr.size();
		ptr.push_back(value);
		return value;
	}

	Face* add_face() { return _add(faces); }
	Edge* add_edge() { return _add(edges); }
	Vert* add_vert() { return _add(verts); }
	HEdge* add_hedge() { return _add(hedges); }
};

struct EdgePartialAdjacency {
	std::vector<vec2i> edges;
	std::map<std::pair<int, int>, int> edgemap;

	template<typename T>
	EdgePartialAdjacency(const std::vector<T>& face) {
		for (auto f : face) {
			for (int i = 0; i < f.size(); i++) {
				int v0idx = f[i];
				int v1idx = f[(i + 1) % f.size()];
				if (edgemap.find(pair<int, int>(v0idx, v1idx)) != edgemap.end()) continue;
				edges.push_back(makevec2i(v0idx, v1idx));
				edgemap[pair<int, int>(v0idx, v1idx)] = edges.size() - 1;
				edgemap[pair<int, int>(v1idx, v0idx)] = edges.size() - 1;
			}
		}
	}

	int edge(int v0, int v1) {
		if (edgemap.find(std::pair<int, int>(v0, v1)) != edgemap.end())
			return edgemap[std::pair<int, int>(v0, v1)];
		else return -1;
	}
};

inline vec3f mesh_face_normal(const vec3i& f, const std::vector<vec3f>& pos) { return triangle_normal(pos[f.x], pos[f.y], pos[f.z]); }
inline vec3f mesh_face_normal(const vec4i& f, const std::vector<vec3f>& pos) { return triangle_normal(pos[f.x], pos[f.y], pos[f.z]); }

template<typename T>
inline std::vector<vec3f> mesh_smooth_normals(const std::vector<T>& face, const std::vector<vec3f>& pos) {
	std::vector<vec3f> norm(pos.size(), zero3f);
	std::vector<float> count(pos.size(), 0);
	for (auto f : face) for (auto vid : f) norm[vid] += mesh_face_normal(f, pos);
	for (vec3f& n : norm) n = normalize(n);
	return norm;
}

template<typename T>
inline std::vector<vec2i> mesh_edges(const std::vector<T>& face) {
	EdgePartialAdjacency adj(face);
	return adj.edges;
}

template<typename T>
inline void mesh_add_vert(T* mesh) {
	mesh->pos.push_back(zero3f);
	if (not mesh->norm.empty()) mesh->norm.push_back(zero3f);
	if (not mesh->uv.empty()) mesh->uv.push_back(zero2f);
}

template<typename T>
inline void mesh_zero_vert(T* mesh, int idx) {
	mesh->pos[idx] = zero3f;
	if (not mesh->norm.empty()) mesh->norm[idx] = zero3f;
	if (not mesh->uv.empty()) mesh->uv[idx] = zero2f;
}

template<typename T, typename VI>
inline void mesh_set_vert(T* mesh, int idx, const T* src, const VI& src_idx) {
	mesh_zero_vert(mesh, idx);
	float w = 1.0f / src_idx.size();
	for (auto vid : src_idx) {
		mesh->pos[idx] += src->pos[vid] * w;
		if (not mesh->norm.empty()) mesh->norm[idx] += src->norm[vid] * w;
		if (not mesh->uv.empty()) mesh->uv[idx] += src->uv[vid] * w;
	}
	if (not mesh->norm.empty()) mesh->norm[idx] = normalize(mesh->norm[idx]);
}

template<typename T, typename VI>
inline void mesh_add_vert(T* mesh, const T* src, const VI& src_idx) {
	mesh_add_vert(mesh);
	mesh_set_vert(mesh, mesh->pos.size() - 1, src, src_idx);
}

template<typename T>
void mesh_copy_verts(T* mesh, const T* src) {
	mesh->pos = src->pos;
	mesh->norm = src->norm;
	mesh->uv = src->uv;
}


// http://www.geometrictools.com/LibMathematics/Distance/Wm5DistPoint3Triangle3.cpp
/*inline vec3f closest_point_point_triangle( const vec3f &p, const vec3f &v0, const vec3f &v1, const vec3f &v2 ) {
	vec3f diff = v0 - p;
	vec3f edge0 = v1 - v0;
	vec3f edge1 = v2 - v0;
	float a00 = edge0 % edge0;
	float a01 = edge0 % edge1;
	float a11 = edge1 % edge1;
	float b0 = diff % edge0;
	float b1 = diff % edge1;
	float c = diff % diff;
	float det = max( 0.00001f, fabs(a00*a11 - a01*a01) );
	float s = a01*b1 - a11*b0;
	float t = a01*b0 - a00*b1;
	float sqrDistance;

	if (s + t <= det)
	{
		if (s < 0.0)
		{
			if (t < 0.0)  // region 4
			{
				if (b0 < 0.0)
				{
					t = 0.0;
					if (-b0 >= a00)
					{
						s = 1.0;
						sqrDistance = a00 + 2.0*b0 + c;
					}
					else
					{
						s = -b0/a00;
						sqrDistance = b0*s + c;
					}
				}
				else
				{
					s = 0.0;
					if (b1 >= 0.0)
					{
						t = 0.0;
						sqrDistance = c;
					}
					else if (-b1 >= a11)
					{
						t = 1.0;
						sqrDistance = a11 + 2.0*b1 + c;
					}
					else
					{
						t = -b1/a11;
						sqrDistance = b1*t + c;
					}
				}
			}
			else  // region 3
			{
				s = 0.0;
				if (b1 >= 0.0)
				{
					t = 0.0;
					sqrDistance = c;
				}
				else if (-b1 >= a11)
				{
					t = 1.0;
					sqrDistance = a11 + 2.0*b1 + c;
				}
				else
				{
					t = -b1/a11;
					sqrDistance = b1*t + c;
				}
			}
		}
		else if (t < 0.0)  // region 5
		{
			t = 0.0;
			if (b0 >= 0.0)
			{
				s = 0.0;
				sqrDistance = c;
			}
			else if (-b0 >= a00)
			{
				s = 1.0;
				sqrDistance = a00 + 2.0*b0 + c;
			}
			else
			{
				s = -b0/a00;
				sqrDistance = b0*s + c;
			}
		}
		else  // region 0
		{
			// minimum at interior point
			float invDet = 1.0f / det;
			s *= invDet;
			t *= invDet;
			sqrDistance = s*(a00*s + a01*t + 2.0*b0) + t*(a01*s + a11*t + 2.0*b1) + c;
		}
	}
	else
	{
		float tmp0, tmp1, numer, denom;

		if (s < 0.0)  // region 2
		{
			tmp0 = a01 + b0;
			tmp1 = a11 + b1;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a00 - 2.0*a01 + a11;
				if (numer >= denom)
				{
					s = 1.0;
					t = 0.0;
					sqrDistance = a00 + 2.0*b0 + c;
				}
				else
				{
					s = numer / denom;
					t = 1.0 - s;
					sqrDistance = s*(a00*s + a01*t + 2.0*b0) +
					t*(a01*s + a11*t + 2.0*b1) + c;
				}
			}
			else
			{
				s = 0.0;
				if (tmp1 <= 0.0)
				{
					t = 1.0;
					sqrDistance = a11 + 2.0*b1 + c;
				}
				else if (b1 >= 0.0)
				{
					t = 0.0;
					sqrDistance = c;
				}
				else
				{
					t = -b1/a11;
					sqrDistance = b1*t + c;
				}
			}
		}
		else if (t < 0.0)  // region 6
		{
			tmp0 = a01 + b1;
			tmp1 = a00 + b0;
			if (tmp1 > tmp0)
			{
				numer = tmp1 - tmp0;
				denom = a00 - 2.0*a01 + a11;
				if (numer >= denom)
				{
					t = 1.0;
					s = 0.0;
					sqrDistance = a11 + 2.0*b1 + c;
				}
				else
				{
					t = numer/denom;
					s = 1.0 - t;
					sqrDistance = s*(a00*s + a01*t + 2.0*b0) +
					t*(a01*s + a11*t + 2.0*b1) + c;
				}
			}
			else
			{
				t = 0.0;
				if (tmp1 <= 0.0)
				{
					s = 1.0;
					sqrDistance = a00 + 2.0*b0 + c;
				}
				else if (b0 >= 0.0)
				{
					s = 0.0;
					sqrDistance = c;
				}
				else
				{
					s = -b0/a00;
					sqrDistance = b0*s + c;
				}
			}
		}
		else  // region 1
		{
			numer = a11 + b1 - a01 - b0;
			if (numer <= 0.0)
			{
				s = 0.0;
				t = 1.0;
				sqrDistance = a11 + 2.0*b1 + c;
			}
			else
			{
				denom = a00 - 2.0*a01 + a11;
				if (numer >= denom)
				{
					s = 1.0;
					t = 0.0;
					sqrDistance = a00 + 2.0*b0 + c;
				}
				else
				{
					s = numer/denom;
					t = 1.0 - s;
					sqrDistance = s*(a00*s + a01*t + 2.0*b0) +
					t*(a01*s + a11*t + 2.0*b1) + c;
				}
			}
		}
	}

	//error_if_not_va( s>=0.0 && s <= 1.0001 && t>=0.0 && t <= 1.0001 && (s+t)<=1.0001, "s=%f, t=%f", s,t );
	return v0 + (s * edge0) + (t * edge1);
}*/

