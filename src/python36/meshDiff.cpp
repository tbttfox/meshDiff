#include <Python.h>
#include <vector>
#include "mesh.h"

/*
 * Implements an example function.
 */
PyDoc_STRVAR(meshDiff_example_doc, "example(obj, number) Example function");

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
std::vector<T> readItems(PyObject *_pts) {
	std::vector<T> out;
	PyObject* pts = PySequence_Fast(_pts, "Expected A Sequence");
	//if (!pts) throw std::logic_error("Expected A Sequence");
	Py_ssize_t len = PySequence_Size(pts);
	out.reserve((size_t)len);
	for (Py_ssize_t i = 0; i < len; i++) {
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

PyObject *meshDiff_meshDiff(PyObject *self, PyObject *args) {
    /* Shared references that do not need Py_DECREF before returning. */
    PyObject *meshAPyPts = NULL, *meshAPyFaces = NULL, *meshBPyPts = NULL, *meshBPyFaces = NULL;
	
    /* Parse positional and keyword arguments */
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

    Py_RETURN_NONE;
}

/*
 * List of functions to add to meshDiff in exec_meshDiff().
 */
static PyMethodDef meshDiff_functions[] = {
    { "meshDiff", (PyCFunction)meshDiff_meshDiff, METH_VARARGS, meshDiff_example_doc },
    { NULL, NULL, 0, NULL } /* marks end of array */
};

/*
 * Initialize meshDiff. May be called multiple times, so avoid
 * using static state.
 */
int exec_meshDiff(PyObject *module) {
    PyModule_AddFunctions(module, meshDiff_functions);

    PyModule_AddStringConstant(module, "__author__", "tyler");
    PyModule_AddStringConstant(module, "__version__", "1.0.0");
    PyModule_AddIntConstant(module, "year", 2018);

    return 0; /* success */
}

/*
 * Documentation for meshDiff.
 */
PyDoc_STRVAR(meshDiff_doc, "The meshDiff module");


static PyModuleDef_Slot meshDiff_slots[] = {
    { Py_mod_exec, exec_meshDiff },
    { 0, NULL }
};

static PyModuleDef meshDiff_def = {
    PyModuleDef_HEAD_INIT,
    "meshDiff",
    meshDiff_doc,
    0,              /* m_size */
    NULL,           /* m_methods */
    meshDiff_slots,
    NULL,           /* m_traverse */
    NULL,           /* m_clear */
    NULL,           /* m_free */
};

PyMODINIT_FUNC PyInit_meshDiff() {
    return PyModuleDef_Init(&meshDiff_def);
}
