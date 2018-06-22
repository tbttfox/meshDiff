#include <Python.h>
#include <vector>

#include "mesh.h"


template <typename T>
T readItem(PyObject *_pt) {
	T out;
	PyObject* pt = PySequence_Fast(_pt, "Expected A Sequence");
	auto len = PySequence_Size(pt);
	out.reserve(len);
	for (Py_ssize_t i = 0; i < len; i++) {
		PyObject *value = PySequence_Fast_GET_ITEM(pt, i);
		out.push_back(PyFloat_AsDouble(value));
	}
	Py_DECREF(pt);
	return out;
}

template <>
std::vector<int> readItem<std::vector<int>>(PyObject *_pt) {
	std::vector<int> out;
	PyObject* pt = PySequence_Fast(_pt, "Expected A Sequence");
	auto len = PySequence_Size(pt);
	out.reserve(len);
	for (Py_ssize_t i = 0; i < len; i++) {
		PyObject *value = PySequence_Fast_GET_ITEM(pt, i);
		out.push_back((int) PyLong_AsLong(value));
	}
	Py_DECREF(pt);
	return out;
}

template <>
vec3f readItem<vec3f>(PyObject*_pt) {
	vec3f out;
	if (PySequence_Size(_pt) < 3)
		throw std::logic_error("Expected sequence of length 3, but didn't get it");

	PyObject *pt = PySequence_Fast(_pt, "Expected A Sequence");
	PyObject *x = PySequence_Fast_GET_ITEM(pt, 0);
	PyObject *y = PySequence_Fast_GET_ITEM(pt, 1);
	PyObject *z = PySequence_Fast_GET_ITEM(pt, 2);

	out = makevec3f((float)PyFloat_AsDouble(x), (float)PyFloat_AsDouble(y), (float)PyFloat_AsDouble(z));
	Py_DECREF(pt);
	return out;
}

template <typename T>
std::vector<T> readItems(PyObject *_pts){
	std::vector<T> out;
	PyObject* pts = PySequence_Fast(_pts, "Expected A Sequence");
	//if (!pts) throw std::logic_error("Expected A Sequence");
	Py_ssize_t len = PySequence_Size(pts);
    out.reserve((size_t) len);
    for(Py_ssize_t i = 0; i < len; i++) {
		PyObject *value = PySequence_Fast_GET_ITEM(pts, i);
		out.push_back(readItem<T>(value));
	}
	Py_DECREF(pts);
	return out;
}

PyObject* buildIndexList(const std::vector<int> &data) {
	PyObject* listObj = PyList_New(data.size());
	if (!listObj) throw std::logic_error("Unable to allocate memory for Python list");
	for (unsigned int i = 0; i < data.size(); i++) {
		PyObject *num = PyLong_FromLong((long)data[i]);
		if (!num) {
			Py_DECREF(listObj);
			throw std::logic_error("Unable to allocate memory for Python list");
		}
		PyList_SET_ITEM(listObj, i, num);
	}
	return listObj;
}









static PyObject * meshDiff(PyObject *self, PyObject *args){
    PyObject *meshAPyPts = NULL, *meshAPyFaces = NULL, *meshBPyPts = NULL, *meshBPyFaces = NULL;
    if (!PyArg_ParseTuple(args, "OOOO", &meshAPyPts, &meshAPyFaces, &meshBPyPts, &meshBPyFaces)) {
        return NULL;
    }
	
	std::vector<vec3f> meshAPoints, meshBPoints;
	std::vector<std::vector<int>> meshAFaces, meshBFaces;

	try {
		meshAPoints = readItems<vec3f>(meshAPyPts);
		meshBPoints = readItems<vec3f>(meshBPyPts);
		meshAFaces = readItems<std::vector<int>>(meshAPyFaces);
		meshBFaces = readItems<std::vector<int>>(meshBPyFaces);
	}
	catch (std::logic_error) { return PyErr_NoMemory(); }


	auto meshA = new Mesh(meshAPoints, meshAFaces);
	auto meshB = new Mesh(meshBPoints, meshBFaces);

	std::vector<Mesh*> meshes;
	meshes.push_back(meshA);
	meshes.push_back(meshB);
	vec3f glob_offset;
	float glob_scale;
	MeshMatch::init_meshes(meshes, glob_offset, glob_scale);

	auto mm = new MeshMatch(meshA, meshB);
	mm->algorithm();
	
	PyObject *vm01, *vm10, *fm01, *fm10;
	try {
		vm01 = buildIndexList(mm->vm01);
		vm10 = buildIndexList(mm->vm10);
		fm01 = buildIndexList(mm->fm01);
		fm10 = buildIndexList(mm->fm10);
	}
	catch (std::logic_error) { return PyErr_NoMemory(); }

	delete mm, meshA, meshB;
	return Py_BuildValue("(NN)(NN)", vm01, vm10, fm01, fm10);
}

static PyMethodDef MeshDiffMethods[] = {
    {"meshDiff",  meshDiff, METH_VARARGS,
     "Get the diff between two meshes"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

PyMODINIT_FUNC
initmeshDiff(void) {
    (void) Py_InitModule("meshDiff", MeshDiffMethods);
}
