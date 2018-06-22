#pragma once

#include <functional>

template<typename Data, typename Key>
class BinaryHeapNode {
	Key myKey;
	Data myData;
	int side = 0;

	BinaryHeapNode<Data, Key> *pNode = nullptr;
	BinaryHeapNode<Data, Key> *lNode = nullptr;
	BinaryHeapNode<Data, Key> *rNode = nullptr;

	BinaryHeapNode(Data d, Key k) : myKey(k), myData(d) {}

	int size() { return 1 + (lNode ? lNode->size() : 0) + (rNode ? rNode->size() : 0); }

	void insert(BinaryHeapNode<Data, Key> *oNode) {
		if (!oNode) return;
		if (!lNode) { lNode = oNode; oNode->pNode = this; return; }
		if (!rNode) { rNode = oNode; oNode->pNode = this; return; }

		if (side == 0) {
			if (lNode->myKey < oNode->myKey) lNode->insert(oNode);
			else {
				oNode->insert(lNode);
				lNode = oNode;
				oNode->pNode = this;
			}
		}
		else {
			if (rNode->myKey < oNode->myKey) rNode->insert(oNode);
			else {
				oNode->insert(rNode);
				rNode = oNode;
				oNode->pNode = this;
			}
		}

		side = (side + 1) % 2;
	}

	void clear_delete(std::function<void(Data)> fn_deldata) {
		if (lNode) { lNode->clear_delete(fn_deldata); delete lNode; lNode = nullptr; }
		if (rNode) { rNode->clear_delete(fn_deldata); delete rNode; rNode = nullptr; }
		pNode = nullptr;
		fn_deldata(myData);
	}
	void clear_delete() {
		if (lNode) { lNode->clear_delete(); delete lNode; lNode = nullptr; }
		if (rNode) { rNode->clear_delete(); delete rNode; rNode = nullptr; }
		pNode = nullptr;
		delete myData;
	}

	BinaryHeapNode<Data, Key> *remove() {
		BinaryHeapNode<Data, Key> *minNode = nullptr;
		BinaryHeapNode<Data, Key> *maxNode = nullptr;

		if (lNode && rNode) {
			if (lNode->myKey < rNode->myKey) {
				minNode = lNode;
				maxNode = rNode;
			}
			else {
				minNode = rNode;
				maxNode = lNode;
			}
		}
		else minNode = !rNode ? lNode : rNode;

		if (minNode) {
			minNode->insert(maxNode);
			minNode->pNode = pNode;
		}

		if (pNode) {
			if (pNode->lNode == this) pNode->lNode = minNode;
			else pNode->rNode = minNode;
		}

		pNode = lNode = rNode = nullptr;

		return minNode;
	}

public:
	bool good = true;
	Key key() const { return myKey; }
	Data data() const { return myData; }
	template <typename D, typename K> friend class BinaryHeap;
};

template<typename Data, typename Key>
class BinaryHeap {
	BinaryHeapNode<Data, Key> *root = nullptr;

public:

	BinaryHeap() {}
	~BinaryHeap() {}

	bool empty() const { return !root; }
	int size() const {
		if (!root) return 0;
		return root->size();
	}

	BinaryHeapNode<Data, Key> *minimum() const { return root; }

	BinaryHeapNode<Data, Key> *insert(Data d, Key k) {
		BinaryHeapNode<Data, Key> *nNode = new BinaryHeapNode<Data, Key>(d, k);
		if (!root) root = nNode;
		else if (root->myKey < k) root->insert(nNode);
		else {
			nNode->insert(root);
			root = nNode;
		}
		return nNode;
	}

	void combine(BinaryHeap *o) {
		if (!o->root) return;
		if (!this->root) root = o->root;
		else if (root->myKey < o->root->myKey) root->insert(o->root);
		else {
			o->root->insert(root);
			root = o->root;
		}
		o->root = nullptr;
	}

	// WARNING!!! DELETES ALL NODES!!!
	void clear_delete(std::function<void(Data)> fn_deldata) {
		if (!root) return;
		root->clear_delete(fn_deldata);
		delete root; root = nullptr;
	}
	void clear_delete() {
		if (!root) return;
		root->clear_delete();
		delete root; root = nullptr;
	}

	BinaryHeapNode<Data, Key> *removeMinimum() {
		if (!root) return nullptr;
		BinaryHeapNode<Data, Key> *rNode = root;
		root = root->remove();
		return rNode;
	}

	BinaryHeapNode<Data, Key> *remove(BinaryHeapNode<Data, Key> *rNode) {
		if (!rNode) return nullptr;
		if (rNode == root) root = rNode->remove();
		else rNode->remove();
		return rNode;
	}

	void remove_delete(BinaryHeapNode<Data, Key> *rNode) {
		if (!rNode) return;
		if (rNode == root) root = rNode->remove();
		else rNode->remove();
		rNode->clear_delete();
		delete rNode;
	}

	void remove_delete(BinaryHeapNode<Data, Key> *rNode, std::function<void(Data)> fn_deldata) {
		if (!rNode) return;
		if (rNode == root) root = rNode->remove();
		else rNode->remove();
		rNode->clear_delete(fn_deldata);
		delete rNode;
	}
};

