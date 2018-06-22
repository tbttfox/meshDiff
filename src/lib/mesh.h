#pragma once

#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <assert.h>

#include "binheap.h"
#include "vec.h"

#include "nanoflann.hpp"
#include "kdtreeAdaptor.h"

typedef int iface;
typedef int ivert;
class Mesh;
typedef KDTreeVectorOfVectorsAdaptor<std::vector<vec3f>, float, 3, nanoflann::metric_L2, int> NKDTree;


inline vec3f triangle_normal(const vec3f& v0, const vec3f& v1, const vec3f& v2) {
	return normalize(cross(v1 - v0, v2 - v0));
}

struct Costs {
	float alpha = 1.0f; //3.0f;
	float dist = 1.0f, norm = 1.0f;  float dist_mult = 0.8f, norm_mult = 0.8f;
	float unm = 1.0f, mis = 1.0f; //4.0f;
	float valence = 0.0f, stretch = 1.0f; //0.1f;
	int knn = 10, srball = 2;
	bool expand_knn = true;
};


// Specialize the hash function for all those std::unordered_sets
// My initial test gave about a 7% speedup
class MGNode;
template <>
struct std::hash<MGNode *> {
	size_t operator()(const MGNode* val) const;
};

class MGNode {
public:
	enum Type { Vertex, Face };
	// back reference
	Mesh *m;
	int i;
	int i_comp;

	// node labels
	Type type;
	vec3f pos;
	vec3f norm;

	// edge (node-node) information
	std::unordered_set<MGNode*> avnodes;
	std::unordered_set<MGNode*> afnodes;
	std::unordered_set<MGNode*> anodes;

	MGNode(Mesh *m, int i, Type type, vec3f pos, vec3f norm) :
		m(m), i(i), type(type), pos(pos), norm(norm) {}
	~MGNode() {}

	void add_edge(MGNode *other) {
		if (this == other)
			return;

		if (other->type == Face)
			afnodes.insert(other);
		else
			avnodes.insert(other);

		if (type == Face)
			other->afnodes.insert(this);
		else
			other->avnodes.insert(this);

		anodes.insert(other);
		other->anodes.insert(this);
	}

	bool has_edge(MGNode *other) {
		return anodes.count(other) != 0;
	}
};

inline size_t std::hash<MGNode *>::operator()(const MGNode* val) const {
	static const size_t shift = (size_t)log2(1 + sizeof(MGNode));
	return (size_t)(val) >> shift;
}



class Mesh {
public:
	std::string                  filename;

	ivert						 n_v;
	std::vector<vec3f>           vpos;
	std::vector<vec3f>           vnorms;
	NKDTree                     *vkdt;
	std::vector<MGNode*>         vnodes;

	iface                             n_f;
	std::vector<std::vector<ivert>>   finds;
	std::vector<vec3f>                fpos;
	std::vector<vec3f>                fnorms;
	NKDTree                      *fkdt;

	std::vector<MGNode*>              fnodes;

	std::unordered_set<MGNode*>  nodes;

	std::vector<std::vector<ivert>>   vknn;
	std::vector<std::vector<iface>>   fknn;
	std::vector<std::vector<MGNode*>> vsrball;
	std::vector<std::vector<MGNode*>> fsrball;

	std::vector<std::unordered_set<MGNode*>> components;

	// exportable constructor
	Mesh(const float _vpos[][3], int vpSize, const int _finds[], int _fcounts[], int fcSize) {
		vpos.reserve(vpSize);
		for (int i = 0; i < vpSize; ++i) {
			vpos.push_back(makevec3f(_vpos[i][0], _vpos[i][1], _vpos[i][2]));
		}
		
		finds.reserve(fcSize);
		int nxt = 0;
		for (int i = 0; i < fcSize; ++i) {
			std::vector<ivert> face;
			for (int c = 0; c < _fcounts[i]; ++c) {
				face.push_back(_finds[nxt++]);
			}
			finds.push_back(face);
		}
	}

	// stl-converting constructor
	template <typename T>
	Mesh(const std::vector<T> &_vpos, const std::vector<std::vector<int>> &_finds):
	finds(_finds) {
		vpos.reserve(_vpos.size());
		for (auto vp : _vpos)
			vpos.push_back(makevec3f(vp[0], vp[1], vp[2]));
	}
	// local constructor
	Mesh(const std::vector<vec3f> &_vpos, const std::vector<std::vector<ivert>> &_finds):
	vpos(_vpos), finds(_finds) { }

	void startup() {
		n_v = vpos.size();
		n_f = finds.size();
		build_extra_data();
		vknn.resize(n_v); vsrball.resize(n_v);
		fknn.resize(n_f); fsrball.resize(n_f);
	}

	void init(float scale, vec3f &o);
	void build_extra_data();

	void buildKDTrees() {
		vkdt = new NKDTree(3, vpos, 16);
		vkdt->index->buildIndex();
		fkdt = new NKDTree(3, fpos, 16);
		fkdt->index->buildIndex();
	}

	void getknn(const NKDTree *nkdt, const vec3f &point, const size_t num_results,std::vector<int> &ret_indexes, std::vector<float> &out_dists_sqr) {
		ret_indexes.resize(num_results);
		out_dists_sqr.resize(num_results);
		nanoflann::KNNResultSet<float, int> resultSet(num_results);
		resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
		float pt[3] = { point[0], point[1], point[2] };
		nkdt->index->findNeighbors(resultSet, &pt[0], nanoflann::SearchParams(10));
	}

	void getknn(const NKDTree *nkdt, const vec3f &point, const size_t num_results, std::vector<int> &ret_indexes) {
		std::vector<float> out_dists_sqr;
		getknn(nkdt, point, num_results, ret_indexes, out_dists_sqr);
	}

	// k-nearest neighbors to self by euclidean distance
	void recompute_knn(int knn);

	// all neighbors within surface "radius" ball
	void recompute_srball(int r);
};


class MeshMatchChange {
public:
	MGNode *n_0;
	MGNode *n_1;
	bool addmatch;
	float cost_delta;

	MeshMatchChange(MGNode *n_0, MGNode *n_1, bool addmatch) :
		n_0(n_0), n_1(n_1), addmatch(addmatch) {}
};


class MeshMatch {
public:
	std::string name;

	Mesh *m0; int n_v0, n_f0;
	Mesh *m1; int n_v1, n_f1;

	// must assert consistency!
	std::vector<int> vm01, vm10;
	std::vector<int> fm01, fm10;

	Costs costs;

	// greedy algorithm variables
	BinaryHeap<MeshMatchChange*, float> *bh_changes = nullptr;
	std::vector<std::vector<int>> v0_knn1, v1_knn0, f0_knn1, f1_knn0;
	bool map_initialized = false; int l_knn = -1, l_srball = -1;
	std::unordered_map<MGNode*, std::unordered_set<BinaryHeapNode<MeshMatchChange*, float>*>> map_node_bhn;
	std::unordered_map<MGNode*, std::unordered_set<MGNode*>> map_node_nodes;

	MeshMatch(Mesh *m0, Mesh *m1) : m0(m0), m1(m1) {
		n_v0 = m0->n_v;
		n_f0 = m0->n_f;
		n_v1 = m1->n_v;
		n_f1 = m1->n_f;
		reset();
	}

	static void init_meshes(std::vector<Mesh*> &meshes);
	static void init_meshes(std::vector<Mesh*> &meshes, vec3f &offset, float &s);

	float compute_cost(MGNode *n_0, MGNode *n_1) const;
	float compute_cost_after_change(MGNode *n_0, MGNode *n_1) const;
	float compute_cost_delta(MeshMatchChange *mmc) const;

	void algorithm();
	void init_greedy();
	void assign_greedy();
	void assign_greedy_step();
	void greedy_backtrack(float per);
	
	void assign_no_geo_cost();

	void unassign_mismatched();
	void unassign_mismatched_faces();
	void unassign_twisted();
	void unassign_verts_with_any_unmatched_faces();
	void unassign_verts_with_all_unmatched_faces();
	void unassign_faces_with_unmatched_verts();
	void unassign_small_patches(int sz);
	void unassign_small_patches(float per);
	void add_potential_changes(MGNode *n, BinaryHeap<MeshMatchChange*, float> *bh);
	void build_change_heap(BinaryHeap<MeshMatchChange*, float>* bh, std::vector<MGNode *> nodes);
	void build_change_heap(BinaryHeap<MeshMatchChange*, float>* bh, std::vector<MGNode *> nodes, int threads);
	std::vector<MeshMatchChange *> MeshMatch::get_potential_changes(MGNode *n) const;
	void add_found_changes(MeshMatchChange *mmc, BinaryHeap<MeshMatchChange*, float> *bh);

	MGNode *get_match(MGNode *n) const;
	bool has_mismatch(MGNode *n);
	void get_mismatches(MGNode *n, std::unordered_set<MGNode*> &ln_o);
	std::vector<MGNode*> MeshMatch::get_mismatches(MGNode *n);
	bool is_twisted(MGNode *n);

	void cleanup() {
		unassign_twisted();
		unassign_mismatched();
		unassign_faces_with_unmatched_verts();
		unassign_verts_with_all_unmatched_faces();
	}
	void reset() {
		vm01.assign(n_v0, -1); fm01.assign(n_f0, -1);
		vm10.assign(n_v1, -1); fm10.assign(n_f1, -1);
	}
	void assert_consistency() {
		bool good = true;
		for (int i_v0 = 0; i_v0 < n_v0; i_v0++) good &= (vm01[i_v0] == -1 || vm10[vm01[i_v0]] == i_v0);
		for (int i_v1 = 0; i_v1 < n_v1; i_v1++) good &= (vm10[i_v1] == -1 || vm01[vm10[i_v1]] == i_v1);
		for (int i_f0 = 0; i_f0 < n_f0; i_f0++) good &= (fm01[i_f0] == -1 || fm10[fm01[i_f0]] == i_f0);
		for (int i_f1 = 0; i_f1 < n_f1; i_f1++) good &= (fm10[i_f1] == -1 || fm01[fm10[i_f1]] == i_f1);
		if (good) return;

		printf("\n\n\nConsistency assertion failed\n");
		printf("v0\n"); for (int i_v0 = 0; i_v0 < n_v0; i_v0++) printf("%d->%d  ", i_v0, vm01[i_v0]);
		printf("v1\n"); for (int i_v1 = 0; i_v1 < n_v1; i_v1++) printf("%d->%d  ", i_v1, vm10[i_v1]);
		printf("f0\n"); for (int i_f0 = 0; i_f0 < n_f0; i_f0++) printf("%d->%d  ", i_f0, fm01[i_f0]);
		printf("f1\n"); for (int i_f1 = 0; i_f1 < n_f1; i_f1++) printf("%d->%d  ", i_f1, fm10[i_f1]);
		assert(false);
	}
	void vassign(int i_v0, int i_v1) {
		assert(vm01[i_v0] == -1 && vm10[i_v1] == -1);
		vm01[i_v0] = i_v1; vm10[i_v1] = i_v0;
	}
	void fassign(int i_f0, int i_f1) {
		assert(fm01[i_f0] == -1 && fm10[i_f1] == -1);
		fm01[i_f0] = i_f1; fm10[i_f1] = i_f0;
	}
	void assign(MGNode *n_0, MGNode *n_1) {
		assert(n_0->type == n_1->type);
		if (n_0->type == MGNode::Vertex) {
			vassign(n_0->i, n_1->i);
		}
		else {
			fassign(n_0->i, n_1->i);
		}
	}
	void vclear0(int i_v0) {
		vclear(i_v0, vm01[i_v0]);
	}
	void vclear1(int i_v1) {
		vclear(vm10[i_v1], i_v1);
	}
	void vclear(int i_v0, int i_v1) {
		assert(i_v0 != -1 && i_v1 != -1 && vm01[i_v0] == i_v1 && vm10[i_v1] == i_v0);
		vm01[i_v0] = -1; vm10[i_v1] = -1;
	}
	void fclear0(int i_f0) {
		fclear(i_f0, fm01[i_f0]);
	}
	void fclear1(int i_f1) {
		fclear(fm10[i_f1], i_f1);
	}
	void fclear(int i_f0, int i_f1) {
		assert(i_f0 != -1 && i_f1 != -1 && fm01[i_f0] == i_f1 && fm10[i_f1] == i_f0);
		fm01[i_f0] = -1; fm10[i_f1] = -1;
	}
	void clear(MGNode *n_0, MGNode *n_1) {
		assert(n_0->type == n_1->type);
		if (n_0->type == MGNode::Vertex) {
			vclear(n_0->i, n_1->i);
		}
		else {
			fclear(n_0->i, n_1->i);
		}
	}
	void clear(MGNode *n) {
		assert(n->m == m0 || n->m == m1);
		if (n->type == MGNode::Vertex) {
			if (n->m == m0) vclear0(n->i);
			else vclear1(n->i);
		}
		else {
			if (n->m == m0) fclear0(n->i);
			else fclear1(n->i);
		}
	}
	void dochange(MeshMatchChange *mmc) {
		if (mmc->addmatch)
			assign(mmc->n_0, mmc->n_1);
		else
			clear(mmc->n_0, mmc->n_1);
	}
	void undochange(MeshMatchChange *mmc) {
		if (mmc->addmatch)
			clear(mmc->n_0, mmc->n_1);
		else
			assign(mmc->n_0, mmc->n_1);
	}
};


