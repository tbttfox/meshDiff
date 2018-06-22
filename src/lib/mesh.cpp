#include "mesh.h"
#include "binheap.h"
#include "vec.h"
#include "omp.h"


void Mesh::build_extra_data() {
	n_v = vpos.size();
	n_f = finds.size();

	// build extra data (face position and normal, vertex normal)
	vnorms.resize(n_v, zero3f);
	fpos.resize(n_f, zero3f);
	fnorms.resize(n_f, zero3f);

	for (int i_f = 0; i_f < n_f; i_f++) {
		std::vector<ivert> &vinds = finds[i_f];
		int n = vinds.size();
		vec3f &fp = fpos[i_f];
		vec3f norm = zero3f;
		for (int i = 0; i < n; i++) {
			vec3f &v0 = vpos[vinds[i]];
			vec3f &v1 = vpos[vinds[(i + 1) % n]];
			vec3f &v2 = vpos[vinds[(i + 2) % n]];
			fp += v0;
		}
		if (n) fp /= (float)n;
		for (int i = 0; i < n; i++) {
			vec3f &v0 = vpos[vinds[i]];
			vec3f &v1 = vpos[vinds[(i + 1) % n]];
			norm += triangle_normal(fp, v0, v1);
		}
		norm = normalize(norm);
		fnorms[i_f] = norm;
		for (int i_v : vinds) 
			vnorms[i_v] += norm;
	}

	for (int i_v = 0; i_v < n_v; i_v++) {
		vnorms[i_v] = normalize(vnorms[i_v]);
	}
}

void Mesh::init(float scale, vec3f &o) {
	for (int i_v = 0; i_v < n_v; i_v++)
		vpos[i_v] = (vpos[i_v] + o) * scale;
	for (int i_f = 0; i_f < n_f; i_f++)
		fpos[i_f] = (fpos[i_f] + o) * scale;

	buildKDTrees();

	vnodes.reserve(n_v);
	fnodes.reserve(n_f);
	nodes.reserve(n_v + n_f);

	for (ivert i_v = 0; i_v < n_v; i_v++) { // ~ 5 sec debug
		MGNode *n = new MGNode(this, i_v, MGNode::Vertex, vpos[i_v], vnorms[i_v]);
		vnodes.push_back(n); nodes.insert(n);
	}
	for (iface i_f = 0; i_f < n_f; i_f++) {// ~ 5 sec debug
		MGNode *n = new MGNode(this, i_f, MGNode::Face, fpos[i_f], fnorms[i_f]);
		fnodes.push_back(n); nodes.insert(n);
	}

	for (iface i_f = 0; i_f < n_f; i_f++) {
		MGNode *f = fnodes[i_f];
		std::vector<ivert> &vinds = finds[i_f];
		int n_vinds = vinds.size();

		// face-vert edges
		for (ivert i_v : finds[i_f])
			vnodes[i_v]->add_edge(f);

		// vert-vert edges
		for (int i0 = 0; i0 < n_vinds; i0++)
			vnodes[vinds[i0]]->add_edge(vnodes[vinds[(i0 + 1) % n_vinds]]);
	}
	// face-face edges
	for (MGNode *v : vnodes) {// ~ 55 sec debug
		std::unordered_set<MGNode*> &fnodes = v->afnodes;
		for (MGNode *av : v->avnodes) {
			std::unordered_set<MGNode*> &afnodes = av->afnodes;
			std::unordered_set<MGNode*> intersection;
			for (MGNode *f : fnodes)
				if (afnodes.count(f))
					intersection.insert(f);
			for (MGNode *f : intersection)
				for (MGNode *af : intersection)
					f->add_edge(af);
		}
	}

	std::unordered_set<MGNode*> left(nodes.begin(), nodes.end());
	while (!left.empty()) {
		std::unordered_set<MGNode*> comp_nodes;
		std::unordered_set<MGNode*> grow;
		MGNode *n = *left.begin(); left.erase(n); grow.insert(n);
		while (!grow.empty()) {
			MGNode *n = *grow.begin();
			grow.erase(n);
			left.erase(n);
			comp_nodes.insert(n);
			for (MGNode *an : n->anodes)
				if (left.count(an))
					grow.insert(an);
		}
		int i_comp = components.size();
		components.push_back(comp_nodes);
		for (MGNode *n : comp_nodes)
			n->i_comp = i_comp;
	}
}

void Mesh::recompute_knn(int knn) {
#   pragma omp parallel for default(none) shared(knn)
	for (ivert i_v = 0; i_v < n_v; i_v++) {
		getknn(vkdt, vpos[i_v], knn, vknn[i_v]);
	}
#   pragma omp parallel for default(none) shared(knn)
	for (iface i_f = 0; i_f < n_f; i_f++) {
		getknn(fkdt, fpos[i_f], knn, fknn[i_f]);
	}
}

void Mesh::recompute_srball(int r) {
#   pragma omp parallel for default(none) shared(r)
	for (int i_v = 0; i_v < n_v; i_v++) {
		MGNode *v = vnodes[i_v];
		std::unordered_set<MGNode*> srball; srball.reserve(150); srball.insert(v);
		std::unordered_set<MGNode*> search; search.reserve(150); search.insert(v);
		for (int d = 0; d < r; d++) {
			std::unordered_set<MGNode*> next; next.reserve(150);
			for (MGNode *n : search) next.insert(n->anodes.begin(), n->anodes.end());
			for (MGNode *n : srball) next.erase(n);
			for (MGNode *n : next) srball.insert(n);
			if (i_v < n_v - 1) { search.clear(); search.insert(next.begin(), next.end()); }
		}
		vsrball[v->i].assign(srball.begin(), srball.end());
	}
#   pragma omp parallel for default(none) shared(r)
	for (int i_f = 0; i_f < n_f; i_f++) {
		MGNode *f = fnodes[i_f];
		std::unordered_set<MGNode*> srball; srball.reserve(150); srball.insert(f);
		std::unordered_set<MGNode*> search; search.reserve(150); search.insert(f);
		for (int d = 0; d < r; d++) {
			std::unordered_set<MGNode*> next; next.reserve(150);
			for (MGNode *n : search) next.insert(n->anodes.begin(), n->anodes.end());
			for (MGNode *n : srball) next.erase(n);
			for (MGNode *n : next) srball.insert(n);
			if (i_f < n_f - 1) { search.clear(); search.insert(next.begin(), next.end()); }
		}
		fsrball[f->i].assign(srball.begin(), srball.end());
	}
}

void MeshMatch::init_meshes(std::vector<Mesh*> &meshes) {
	vec3f offset; float s;
	init_meshes(meshes, offset, s);
}

void MeshMatch::init_meshes(std::vector<Mesh*> &meshes, vec3f &offset, float &s) {
	float el = 0.0f; int n = 0;

	for (int i = 0; i < meshes.size(); ++i)
		meshes[i]->startup();

	vec3f pmin = meshes[0]->vpos[0], pmax = meshes[0]->vpos[0];

#   pragma omp parallel for
	for (int i = 0; i < meshes.size(); ++i) {
		Mesh *m = meshes[i];
		float xx;
		for (std::vector<ivert> &vinds : m->finds) {
			int c = vinds.size();
			for (int i = 0; i < c; i++) {
				xx = length(m->vpos[vinds[i]] - m->vpos[vinds[(i + 1) % c]]);
#				pragma omp atomic
				el += xx;
			}
#			pragma omp atomic
			n += c;
		}
	}

	for (Mesh *m : meshes) {
		for (vec3f &p : m->vpos) {
			pmin = min_component(pmin, p); pmax = max_component(pmax, p);
		}
	}

	s = (float)n / el;
	offset = -(pmin + pmax) / 2.0f;

#   pragma omp parallel for shared(s) shared(offset)
	for (int i = 0; i < meshes.size(); ++i)
		meshes[i]->init(s, offset);

}

float MeshMatch::compute_cost(MGNode *n_0, MGNode *n_1) const {
	float cost = 0.0f;

	assert(n_0->type == n_1->type);
	assert(n_0->m == m0 && n_1->m == m1);

	MGNode *n_01 = get_match(n_0), *n_10 = get_match(n_1);
	assert((n_01 == n_1 && n_10 == n_0) || (!n_01 && !n_10));

	float sz0 = n_0->anodes.size();
	float sz1 = n_1->anodes.size();
	vec3f &p_0 = n_0->pos;
	vec3f &p_1 = n_1->pos;

	if (!n_01 && !n_10) {                                           // unmatched
		cost += costs.alpha * 2.0f;
	}
	else {
		float dist = length(p_0 - p_1);
		float norm = dot(n_0->norm, n_1->norm);
		cost += costs.dist * dist / (dist + 1.0f);
		cost += costs.norm * (1.0f - norm);

		cost += costs.valence * fabs(sz0 - sz1) / (float)(sz0 + sz1);
	}

	for (MGNode *a_0 : n_0->anodes) {
		MGNode *a_1 = get_match(a_0);
		float sz = sz0 + a_0->anodes.size();
		if (!n_01 || !a_1) {
			cost += costs.unm / sz;
		}
		else {
			if (!n_1->has_edge(a_1)) {
				cost += costs.mis / sz;
			}

			float d0 = length(p_0 - a_0->pos), d1 = length(p_1 - a_1->pos);
			cost += costs.stretch * fabs(d0 - d1) / (d0 + d1 + 0.0001f) / sz;
		}
	}


	for (MGNode *a_1 : n_1->anodes) {
		MGNode *a_0 = get_match(a_1);
		float sz = sz1 + a_1->anodes.size();
		if (!n_10 || !a_0) {
			cost += costs.unm / sz;
		}
		else {
			if (!n_0->has_edge(a_0)) {
				cost += costs.mis / sz;
			}

			float d0 = length(p_0 - a_0->pos), d1 = length(p_1 - a_1->pos);
			cost += costs.stretch * fabs(d0 - d1) / (d0 + d1 + 0.0001f) / sz;
		}
	}

	return cost;
}

float MeshMatch::compute_cost_after_change(MGNode *n_0, MGNode *n_1) const {
	float cost = 0.0f;

	assert(n_0->type == n_1->type);
	assert(n_0->m == m0 && n_1->m == m1);

	MGNode *n_01 = get_match(n_0), *n_10 = get_match(n_1);
	assert((n_01 == n_1 && n_10 == n_0) || (!n_01 && !n_10));

	float sz0 = n_0->anodes.size();
	float sz1 = n_1->anodes.size();
	vec3f &p_0 = n_0->pos;
	vec3f &p_1 = n_1->pos;

	if (n_01 && n_10) {
		cost += costs.alpha * 2.0f;
	}
	else {
		float dist = length(p_0 - p_1);
		float norm = dot(n_0->norm, n_1->norm);
		cost += costs.dist * dist / (dist + 1.0f);
		cost += costs.norm * (1.0f - norm);

		cost += costs.valence * fabs(sz0 - sz1) / (float)(sz0 + sz1);
	}

	for (MGNode *a_0 : n_0->anodes) {
		MGNode *a_1 = get_match(a_0);
		float sz = sz0 + a_0->anodes.size();
		if (n_01 || !a_1) {
			cost += costs.unm / sz;
		}
		else {
			if (!n_1->has_edge(a_1)) {
				cost += costs.mis / sz;
			}

			float d0 = length(p_0 - a_0->pos), d1 = length(p_1 - a_1->pos);
			cost += costs.stretch * fabs(d0 - d1) / (d0 + d1 + 0.0001f) / sz;
		}
	}

	for (MGNode *a_1 : n_1->anodes) {
		MGNode *a_0 = get_match(a_1);
		float sz = sz1 + a_1->anodes.size();
		if (n_10 || !a_0) {
			cost += costs.unm / sz;
		}
		else {
			if (!n_0->has_edge(a_0)) {
				cost += costs.mis / sz;
			}

			float d0 = length(p_0 - a_0->pos), d1 = length(p_1 - a_1->pos);
			cost += costs.stretch * fabs(d0 - d1) / (d0 + d1 + 0.0001f) / sz;
		}
	}

	return cost;
}

float MeshMatch::compute_cost_delta(MeshMatchChange *mmc) const {
	mmc->cost_delta = compute_cost_after_change(mmc->n_0, mmc->n_1);
	mmc->cost_delta -= compute_cost(mmc->n_0, mmc->n_1);
	return mmc->cost_delta;
}

void MeshMatch::build_change_heap(BinaryHeap<MeshMatchChange*, float>* bh, std::vector<MGNode *> nodes) {
	build_change_heap(bh, nodes, omp_get_max_threads());
}

void MeshMatch::build_change_heap(BinaryHeap<MeshMatchChange*, float>* bh, std::vector<MGNode *> nodes, int threads){
	std::vector<std::vector<MeshMatchChange *>> changes;
	changes.resize(nodes.size());
	// 10 threads chosen by testing. The general idea is that nodes.size() 
	// is usually 18 if we're working on a quadded mesh
	// This way we get a relatively saturated and even split on the threads
	// without too much extra overhead. Saves about half a second on my tests

#	pragma omp parallel for num_threads(threads)
	for (int i = 0; i<nodes.size(); ++i)
		changes[i] = get_potential_changes(nodes[i]);

	for (std::vector<MeshMatchChange *> &mmcv : changes)
		for (MeshMatchChange *mmc : mmcv)
			add_found_changes(mmc, bh);

}




void MeshMatch::add_potential_changes(MGNode *n, BinaryHeap<MeshMatchChange*, float> *bh){
	for (auto && mmc : get_potential_changes(n))
		add_found_changes(mmc, bh);
}

std::vector<MeshMatchChange *> MeshMatch::get_potential_changes(MGNode *n) const {
	int i = n->i;
	Mesh *m = n->m;
	MGNode::Type t = n->type;

	std::vector<MeshMatchChange *> out;
	auto fn_add = [&](MGNode *n, MGNode *n_, bool add) {
		if (n->m == m1) std::swap(n, n_);

		auto f = map_node_nodes.find(n);
		if (f == map_node_nodes.end()) return;
		if (f->second.count(n_)) return;

		MeshMatchChange *mmc = new MeshMatchChange(n, n_, add);
		compute_cost_delta(mmc);
		if (mmc->cost_delta > 0.0f) { delete mmc; return; }
		out.push_back(mmc);
	};

	MGNode *n_ = get_match(n);

	if (n_) {
		fn_add(n, n_, false);
		return out;
	}

	// add knn as potential matches
	const std::vector<std::vector<int>> *knn_;
	const std::vector<MGNode*> *nodes_;
	const std::vector<int> *match;
	if (m == m0) {
		if (t == MGNode::Vertex) {
			knn_ = &v0_knn1;
			nodes_ = &m1->vnodes;
			match = &vm10;
		}
		else {
			knn_ = &f0_knn1;
			nodes_ = &m1->fnodes;
			match = &fm10;
		}
	}
	else {
		if (t == MGNode::Vertex) {
			knn_ = &v1_knn0;
			nodes_ = &m0->vnodes;
			match = &vm01;
		}
		else {
			knn_ = &f1_knn0;
			nodes_ = &m0->fnodes;
			match = &fm01;
		}
	}

	for (int i_ : (*knn_)[i]) {
		if ((*match)[i_] == -1) fn_add(n, (*nodes_)[i_], true);
	}

	// add srball and expanded knn as potential matches
	for (MGNode *an : n->anodes) {
		MGNode *an_ = get_match(an); if (!an_) continue;

		Mesh *m_ = an_->m;
		bool isvert = (an_->type == MGNode::Vertex);
		std::vector<MGNode*> &ln_pot = ((isvert ? m_->vsrball : m_->fsrball)[an_->i]);
		std::vector<int> &li_pot_knn = ((isvert ? m_->vknn : m_->fknn)[an_->i]);
		std::vector<MGNode*> &ln_pot_knn = (isvert ? m_->vnodes : m_->fnodes);

		for (MGNode *n_ : ln_pot) {
			if (n_->type != n->type) continue;         // wrong type
			if (get_match(n_)) continue;               // already matched
			fn_add(n, n_, true);
		}

		if (costs.expand_knn) {
			for (int i_ : li_pot_knn) {
				MGNode *n_ = ln_pot_knn[i_];
				if (n_->type != n->type) continue;         // wrong type
				if (get_match(n_)) continue;               // already matched
				fn_add(n, n_, true);
			}
		}
	}
	return out;
}

void MeshMatch::add_found_changes(MeshMatchChange *mmc, BinaryHeap<MeshMatchChange*, float> *bh) {
	auto chg = bh->insert(mmc, mmc->cost_delta);
	map_node_bhn[mmc->n_0].insert(chg);
	map_node_bhn[mmc->n_1].insert(chg);
	map_node_nodes[mmc->n_0].insert(mmc->n_1);
	map_node_nodes[mmc->n_1].insert(mmc->n_0);
};

void MeshMatch::assign_no_geo_cost() {
	int knn = 15;
	float epsilon = 0.01f;
	std::vector<int> li_knn;

	int count = 0;

	for (int i_v0 = 0; i_v0 < n_v0; i_v0++) {
		if (vm01[i_v0] != -1) continue;
		vec3f &v0 = m0->vpos[i_v0];
		vec3f &n0 = m0->vnorms[i_v0];

		m1->getknn(m1->vkdt, v0, knn, li_knn);
		
		for (int i_v1 : li_knn) {
			if (vm10[i_v1] != -1) continue;
			float cost = 0.0f;
			float dist = length(v0 - m1->vpos[i_v1]);
			float norm = dot(n0, m1->vnorms[i_v1]);
			cost += costs.dist * dist / (dist + 1.0f);
			cost += costs.norm * (1.0f - norm);
			if (cost > epsilon) continue;
			count++;
			vassign(i_v0, i_v1); break;
		}
	}
	for (int i_f0 = 0; i_f0 < n_f0; i_f0++) {
		if (fm01[i_f0] != -1) continue;
		vec3f &f0 = m0->fpos[i_f0];
		vec3f &n0 = m0->fnorms[i_f0];

		m1->getknn(m1->fkdt, f0, knn, li_knn);
		for (int i_f1 : li_knn) {
			if (fm10[i_f1] != -1) continue;
			float cost = 0.0f;
			float dist = length(f0 - m1->fpos[i_f1]);
			float norm = dot(n0, m1->fnorms[i_f1]);
			cost += costs.dist * dist / (dist + 1.0f);
			cost += costs.norm * (1.0f - norm);
			if (cost > epsilon) continue;
			count++;
			fassign(i_f0, i_f1); break;
		}
	}
}

void MeshMatch::unassign_mismatched() {
	std::unordered_set<MGNode*> nodes;
	auto gamm = [this, &nodes](std::vector<MGNode *> nns) {
		std::vector<std::vector<MGNode *> > nodevec;
		nodevec.resize(nns.size());
#		pragma omp parallel for
		for (int i = 0; i < nns.size(); ++i) {
			auto &n = nns[i];
			nodevec[i] = get_mismatches(n);
		}
		for (auto &nv : nodevec)
			nodes.insert(nv.begin(), nv.end());
	};

	gamm(m0->vnodes);
	gamm(m0->fnodes);
	gamm(m1->vnodes);
	gamm(m1->fnodes);

	for (MGNode *n : nodes)
		if (get_match(n))
			clear(n);

}

void MeshMatch::unassign_mismatched_faces() {
	/*std::unordered_set<MGNode*> nodes;
	for( MGNode *n : m0->fnodes ) {
		for( MGNode *nv : n->avnodes ) {

		}
		for( MGNode *an : n->anodes ) {
			if( n->type == an->type ) continue;
			MGNode *an_o = get_match(an); if( !an_o ) continue;
			if( !n_o->has_edge(an_o) ) {
				//ln_o.insert(n); ln_o.insert(n_o); ln_o.insert(an_o); ln_o.insert(an);
				if( n->type == MGNode::Face ) { ln_o.insert(n); ln_o.insert(n_o); }
				else { ln_o.insert(an_o); ln_o.insert(an); }
			}
		}

	}

	for( MGNode *n : m0->vnodes ) get_mismatches(n,nodes);
	for( MGNode *n : m0->fnodes ) get_mismatches(n,nodes);
	for( MGNode *n : m1->vnodes ) get_mismatches(n,nodes);
	for( MGNode *n : m1->fnodes ) get_mismatches(n,nodes);

	for( MGNode *n : nodes ) {
		if( get_match(n) ) clear(n);
	}*/
}

void MeshMatch::unassign_twisted() {
	std::unordered_set<MGNode*> nodes;
	for (MGNode *n : m0->fnodes) if (is_twisted(n)) {
		//nodes.insert(n);
		nodes.insert(n->avnodes.begin(), n->avnodes.end());
	}

	for (MGNode *n : nodes) {
		if (get_match(n)) clear(n);
	}
}

void MeshMatch::unassign_verts_with_any_unmatched_faces() {
	std::unordered_set<MGNode*> nodes;
	for (MGNode *n : m0->vnodes) {
		MGNode *n_ = get_match(n);
		if (!n_) continue;

		for (MGNode *n_f : n->afnodes) if (!get_match(n_f)) {
			nodes.insert(n);
			nodes.insert(n_);
			break;
		}
		for (MGNode *n_f : n_->afnodes) if (!get_match(n_f)) {
			nodes.insert(n);
			nodes.insert(n_);
			break;
		}
	}

	for (MGNode *n : nodes) {
		if (get_match(n)) clear(n);
	}
}

void MeshMatch::unassign_verts_with_all_unmatched_faces() {
	std::unordered_set<MGNode*> nodes;
	for (MGNode *n : m0->vnodes) {
		MGNode *n_ = get_match(n);
		if (!n_) continue;

		bool all = true;
		for (MGNode *n_f : n->afnodes) if (get_match(n_f)) { all = false; break; }
		for (MGNode *n_f : n_->afnodes) if (get_match(n_f)) { all = false; break; }
		if (all) { nodes.insert(n); nodes.insert(n_); }
	}

	for (MGNode *n : nodes) {
		if (get_match(n)) clear(n);
	}
}

void MeshMatch::unassign_faces_with_unmatched_verts() {
	std::unordered_set<MGNode*> nodes;
	for (MGNode *n : m0->fnodes) {
		if (!get_match(n)) continue;
		bool any = false;
		for (MGNode *v : n->avnodes) {
			if (!get_match(v)) { any = true; break; }
		}
		if (any) nodes.insert(n);
	}
	for (MGNode *n : m1->fnodes) {
		if (!get_match(n)) continue;
		bool any = false;
		for (MGNode *v : n->avnodes) {
			if (!get_match(v)) { any = true; break; }
		}
		if (any) nodes.insert(n);
	}

	for (MGNode *n : nodes) {
		if (get_match(n)) clear(n);
	}
}

void MeshMatch::unassign_small_patches(int sz) {
	std::unordered_set<MGNode*> nodes;

	for (MGNode *n : m0->fnodes) if (get_match(n)) nodes.insert(n);
	while (!nodes.empty()) {
		std::unordered_set<MGNode*> patch;
		std::unordered_set<MGNode*> grow;

		MGNode *n = *nodes.begin();
		grow.insert(n);
		while (!grow.empty()) {
			MGNode *n = *grow.begin(); grow.erase(n);
			if (patch.count(n)) continue;
			patch.insert(n);
			for (MGNode *an : n->afnodes) if (nodes.count(an)) grow.insert(an);
		}

		for (MGNode *n : patch) nodes.erase(n);

		if (patch.size() >= sz) continue;

		for (MGNode *f : patch) {
			for (MGNode *v : f->avnodes) if (get_match(v)) clear(v);
			clear(f);
		}
	}

	for (MGNode *n : m1->fnodes) if (get_match(n)) nodes.insert(n);
	while (!nodes.empty()) {
		std::unordered_set<MGNode*> patch;
		std::unordered_set<MGNode*> grow;

		MGNode *n = *nodes.begin();
		grow.insert(n);
		while (!grow.empty()) {
			MGNode *n = *grow.begin(); grow.erase(n);
			if (patch.count(n)) continue;
			patch.insert(n);
			for (MGNode *an : n->afnodes) if (nodes.count(an)) grow.insert(an);
		}

		for (MGNode *n : patch) nodes.erase(n);

		if (patch.size() >= sz) continue;

		for (MGNode *f : patch) {
			for (MGNode *v : f->avnodes) if (get_match(v)) clear(v);
			clear(f);
		}
	}
}

void MeshMatch::unassign_small_patches(float per) {
	std::unordered_set<MGNode*> nodes;

	for (MGNode *n : m0->fnodes) if (get_match(n)) nodes.insert(n);
	while (!nodes.empty()) {
		std::unordered_set<MGNode*> patch;
		std::unordered_set<MGNode*> grow;

		MGNode *n = *nodes.begin();
		grow.insert(n);
		while (!grow.empty()) {
			MGNode *n = *grow.begin(); grow.erase(n);
			if (patch.count(n)) continue;
			patch.insert(n);
			for (MGNode *an : n->afnodes) if (nodes.count(an)) grow.insert(an);
		}

		for (MGNode *n : patch) nodes.erase(n);

		int i_comp = (*patch.begin())->i_comp;
		int sz_comp = m0->components[i_comp].size();
		if (patch.size() >= (int)((float)sz_comp*per)) continue;

		for (MGNode *f : patch) {
			//for( MGNode *v : f->avnodes ) if(get_match(v)) clear(v);
			clear(f);
		}
	}

	for (MGNode *n : m1->fnodes) if (get_match(n)) nodes.insert(n);
	while (!nodes.empty()) {
		std::unordered_set<MGNode*> patch;
		std::unordered_set<MGNode*> grow;

		MGNode *n = *nodes.begin();
		grow.insert(n);
		while (!grow.empty()) {
			MGNode *n = *grow.begin(); grow.erase(n);
			if (patch.count(n)) continue;
			patch.insert(n);
			for (MGNode *an : n->afnodes) if (nodes.count(an)) grow.insert(an);
		}

		for (MGNode *n : patch) nodes.erase(n);

		int i_comp = (*patch.begin())->i_comp;
		int sz_comp = m1->components[i_comp].size();
		if (patch.size() >= (int)((float)sz_comp*per)) continue;

		for (MGNode *f : patch) {
			//for( MGNode *v : f->avnodes ) if(get_match(v)) clear(v);
			clear(f);
		}
	}

	unassign_verts_with_all_unmatched_faces();
}

void MeshMatch::algorithm() {
	if (true) {
		reset();
		assign_no_geo_cost();
	}
	init_greedy();

	costs.dist_mult = costs.norm_mult = 0.75f;

	assign_greedy();
}

void MeshMatch::assign_greedy() {
	float d = costs.dist;
	float n = costs.norm;

	costs.dist = 1.0f;
	costs.norm = 1.0f;
	
	for (float per = 0.08f; per >= 0.01f; per *= 0.5f) {
		assign_greedy_step();
		cleanup();
		unassign_small_patches(per);

		if (true) {
			costs.dist *= costs.dist_mult;
			costs.norm *= costs.norm_mult;
		}
	}

	costs.dist = d;
	costs.norm = n;
}

void MeshMatch::init_greedy() {
	if (l_srball != costs.srball) {
		m0->recompute_srball(costs.srball);
		m1->recompute_srball(costs.srball);

		l_srball = costs.srball;
	}

	if (l_knn != costs.knn) {

		v0_knn1.resize(n_v0); f0_knn1.resize(n_f0);
		v1_knn0.resize(n_v1); f1_knn0.resize(n_f1);

#       pragma omp parallel for default(none)
		for (int i_v0 = 0; i_v0 < n_v0; i_v0++)
			m1->getknn(m1->vkdt, m0->vpos[i_v0], costs.knn, v0_knn1[i_v0]);

#       pragma omp parallel for default(none)
		for (int i_f0 = 0; i_f0 < n_f0; i_f0++)
			m1->getknn(m1->fkdt, m0->fpos[i_f0], costs.knn, f0_knn1[i_f0]);

#       pragma omp parallel for default(none)
		for (int i_v1 = 0; i_v1 < n_v1; i_v1++)
			m0->getknn(m0->vkdt, m1->vpos[i_v1], costs.knn, v1_knn0[i_v1]);

#       pragma omp parallel for default(none)
		for (int i_f1 = 0; i_f1 < n_f1; i_f1++)
			m0->getknn(m0->fkdt, m1->fpos[i_f1], costs.knn, f1_knn0[i_f1]);

		m0->recompute_knn(costs.knn);
		m1->recompute_knn(costs.knn);

		l_knn = costs.knn;
	}

	if (!map_initialized) {
		map_node_bhn.reserve(m0->n_v + m0->n_f + m1->n_v + m1->n_f);
		map_node_nodes.reserve(m0->n_v + m0->n_f + m1->n_v + m1->n_f);
		for (auto n : m0->vnodes) { map_node_bhn[n].reserve(75); map_node_nodes[n].reserve(75); }
		for (auto n : m0->fnodes) { map_node_bhn[n].reserve(75); map_node_nodes[n].reserve(75); }
		for (auto n : m1->vnodes) { map_node_bhn[n].reserve(75); map_node_nodes[n].reserve(75); }
		for (auto n : m1->fnodes) { map_node_bhn[n].reserve(75); map_node_nodes[n].reserve(75); }
		map_initialized = true;
	}
}

void MeshMatch::assign_greedy_step() {
	init_greedy();

	float delta = 0.0f;

	bh_changes = new BinaryHeap<MeshMatchChange*, float>();

	build_change_heap(bh_changes, m0->vnodes);
	build_change_heap(bh_changes, m1->vnodes);
	build_change_heap(bh_changes, m0->fnodes);
	build_change_heap(bh_changes, m1->fnodes);

	int iter = 0;
	double t_rem = 0.0f;
	double t_add = 0.0f;
	double s;
	while (!bh_changes->empty()) {
		auto bn_min = bh_changes->minimum();

		float c = bn_min->key();
		if (c > 0.0f) break;
		delta += c;

		MeshMatchChange *mmc = bn_min->data();

		// find nodes that are affected by mmc
		std::unordered_set<MGNode*> nodes;
		nodes.insert(mmc->n_0); nodes.insert(mmc->n_0->anodes.begin(), mmc->n_0->anodes.end());
		nodes.insert(mmc->n_1); nodes.insert(mmc->n_1->anodes.begin(), mmc->n_1->anodes.end());
		std::vector<MGNode *> nodeVec(nodes.begin(), nodes.end());
		// perform change
		dochange(mmc);

		// remove old changes involving nodes affected by mmc
		for (auto n : nodeVec) {
			for (auto bhn : map_node_bhn[n]) {
				auto mmc = bhn->data();
				if (mmc->n_0 == n) {
					map_node_bhn[mmc->n_1].erase(bhn);
					map_node_nodes[mmc->n_1].erase(n);
				}
				else {
					map_node_bhn[mmc->n_0].erase(bhn);
					map_node_nodes[mmc->n_0].erase(n);
				}
				bh_changes->remove_delete(bhn);
			}
			map_node_bhn[n].clear();
			map_node_nodes[n].clear();
		}

		auto bh_local = new BinaryHeap<MeshMatchChange*, float>();

		// add new potential changes for nodes affected by mmc
		//for (auto n : nodes)add_potential_changes(n, bh_local);
		build_change_heap(bh_local, nodeVec, 10);

		bh_changes->combine(bh_local);
		delete bh_local;

		iter++;
	}
	

	auto clearNodes = [&](std::vector<MGNode *> nodes) {
#		pragma omp parallel for
		for (int i = 0; i<nodes.size(); ++i) {
			map_node_bhn[nodes[i]].clear();
			map_node_nodes[nodes[i]].clear();
		}
	};
	clearNodes(m0->vnodes);
	clearNodes(m0->fnodes);
	clearNodes(m1->vnodes);
	clearNodes(m1->fnodes);

	bh_changes->clear_delete();
	delete bh_changes;
	bh_changes = nullptr;
}

void MeshMatch::greedy_backtrack(float per) {
	//unassign_twisted();
	//unassign_mismatched();
	//unassign_faces_with_unmatched_verts();
	//unassign_verts_with_all_unmatched_faces();
	unassign_small_patches(per);
}

MGNode* MeshMatch::get_match(MGNode *n) const {
	assert(n->m == m0 || n->m == m1);

	if (n->type == MGNode::Vertex) {
		if (n->m == m0) {
			int i1 = vm01[n->i];
			if (i1 == -1) return nullptr;
			else return m1->vnodes[i1];
		}
		else {
			int i0 = vm10[n->i];
			if (i0 == -1) return nullptr;
			else return m0->vnodes[i0];
		}
	}
	else {
		if (n->m == m0) {
			int i1 = fm01[n->i];
			if (i1 == -1) return nullptr;
			else return m1->fnodes[i1];
		}
		else {
			int i0 = fm10[n->i];
			if (i0 == -1) return nullptr;
			else return m0->fnodes[i0];
		}

	}
}

bool MeshMatch::has_mismatch(MGNode *n) {
	assert(n->m == m0 || n->m == m1);
	MGNode *n_o = get_match(n);
	if (!n_o) return false;
	for (MGNode *an : n->anodes) {
		MGNode *an_o = get_match(an);
		if (!an_o) continue;
		if (!n_o->has_edge(an_o)) return true;
	}
	return false;
}



std::vector<MGNode*> MeshMatch::get_mismatches(MGNode *n) {
	assert(n->m == m0 || n->m == m1);
	std::vector<MGNode*> ln_o;
	MGNode *n_o = get_match(n);
	if (!n_o)
		return ln_o;
	for (MGNode *an : n->anodes) {
		if (n->type == an->type)
			continue;
		MGNode *an_o = get_match(an);
		if (!an_o)
			continue;
		if (!n_o->has_edge(an_o)) {
			if (n->type == MGNode::Face) {
				ln_o.push_back(n);
				ln_o.push_back(n_o);
			}
			else {
				ln_o.push_back(an_o);
				ln_o.push_back(an);
			}
		}
	}
	return ln_o;
}


void MeshMatch::get_mismatches(MGNode *n, std::unordered_set<MGNode*> &ln_o) {
	assert(n->m == m0 || n->m == m1);
	MGNode *n_o = get_match(n);
	if (!n_o)
		return;
	for (MGNode *an : n->anodes) {
		if (n->type == an->type)
			continue;
		MGNode *an_o = get_match(an);
		if (!an_o)
			continue;
		if (!n_o->has_edge(an_o)) {
			if (n->type == MGNode::Face) {
				ln_o.insert(n);
				ln_o.insert(n_o);
			}
			else {
				ln_o.insert(an_o);
				ln_o.insert(an);
			}
		}
	}
}

bool MeshMatch::is_twisted(MGNode *n) {
	assert(n->m == m0 || n->m == m1);
	assert(n->type == MGNode::Face);

	MGNode *n_ = get_match(n); if (!n_) return false;

	std::vector<ivert> &vm_0 = (n->m == m0 ? vm01 : vm10);

	// tranlate vertex inds to other mesh
	std::vector<ivert> &lvi0 = (n->m == m0 ? m0->finds[n->i] : m1->finds[n->i]);
	std::vector<ivert> lvi_0; for (ivert i : lvi0) lvi_0.push_back(vm_0[i]);
	std::vector<ivert> &lvi_1 = (n_->m == m0 ? m0->finds[n_->i] : m1->finds[n_->i]);

	if (lvi_0.size() != lvi_1.size()) return false;
	for (ivert i : lvi_0) {
		bool f = false;
		for (ivert i_ : lvi_1) if (i == i_) { f = true; break; }
		if (!f) return false;
	}
	for (ivert i : lvi_0) if (i == -1) return false;
	for (ivert i : lvi_1) if (i == -1) return false;

	// check if inds are out-of-order, but allow for rotation
	int nv = lvi_0.size();
	bool all;
	for (int i0 = 0; i0 < nv; i0++) {
		all = true;
		for (int io = 0; io < nv; io++) {
			int i1 = (i0 + io) % nv;
			if (lvi_0[io] != lvi_1[i1]) { all = false; break; }
		}
		if (all) return false;
		all = true;
		for (int io = 0; io < nv; io++) {
			int i1 = (nv + i0 - io) % nv;
			if (lvi_0[io] != lvi_1[i1]) { all = false; break; }
		}
		if (all) return false;
	}
	return true;
}


