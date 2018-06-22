#include <string>
#include <iostream>
#include <fstream>

#include "mesh.h"


inline bool exists_test(const std::string& name) {
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	}
	else {
		return false;
	}
}

// OBJ loading utility functions
std::vector<std::string> splitstring(const std::string &s, const char &delim) {
	std::vector<std::string> out;
	if (s.empty()) return out;

	auto start = 0U;
	auto end = s.find(delim);
	while (end != std::string::npos) {
		out.push_back(s.substr(start, end - start));
		start = end + 1;
		end = s.find(delim, start);
	}
	out.push_back(s.substr(start, end));
	return out;
}

void loadObjStrings(const std::string & filename, std::vector<std::vector<int>> &finds, std::vector<vec3f> &vpos) {
	std::string line;
	std::ifstream f(filename);
	if (!f.is_open()) {
		perror("error while opening file");
	}

	std::vector<std::string> sp, fsp;
	std::vector<ivert> ff;

	while (getline(f, line)) {
		if (line.empty()) { continue; }
		switch (line[0]) {
		case 'v':
			if (line[1] == ' ') {
				sp = splitstring(line, ' ');
				vpos.push_back(makevec3f(std::stof(sp[1]), std::stof(sp[2]), std::stof(sp[3])));
			}
			break;
		case 'f':
			sp = splitstring(line, ' ');
			ff.clear();
			for (size_t i = 1; i<sp.size(); ++i) {
				if (sp[i].empty()) continue;
				fsp = splitstring(sp[i], '/');
				ff.push_back(std::stoi(fsp[0]) - 1);
			}
			finds.push_back(ff);
			break;
		}
	}

	if (f.bad())
		perror("error while reading file");
}

void loadBinObj(const std::string & filename, std::vector<std::vector<int>> &finds, std::vector<vec3f> &vpos) {
	std::ifstream f(filename, std::ios::in | std::ios::binary);
	if (!f.is_open()) {
		perror("error while opening file");
	}
	int size, count, idx, i, j;
	float x, y, z;

	f.read((char *)&size, sizeof(size));
	for (i = 0; i < size; ++i) {
		f.read((char *)&x, sizeof(x));
		f.read((char *)&y, sizeof(y));
		f.read((char *)&z, sizeof(z));
		vpos.push_back(makevec3f(x, y, z));
	}

	f.read((char *)&size, sizeof(size));
	for (i = 0; i < size; ++i) {
		f.read((char *)&count, sizeof(count));
		std::vector<int> face;
		for (j = 0; j < count; ++j) {
			f.read((char *)&idx, sizeof(idx));
			face.push_back(idx);
		}
		finds.push_back(face);
	}

	if (f.bad())
		perror("error while reading file");
}

void saveBinObj(const std::string & filename, const std::vector<std::vector<int>> &finds, const std::vector<vec3f> &vpos) {
	std::ofstream f(filename, std::ios::out | std::ios::binary);
	if (!f.is_open()) {
		perror("error while opening file");
	}

	int vps = (int)vpos.size();
	f.write((char *)&vps, sizeof(vps));
	for (int i = 0; i < vps; ++i) {
		f.write((char*) &(vpos[i].x), sizeof(vpos[i].x));
		f.write((char*) &(vpos[i].y), sizeof(vpos[i].y));
		f.write((char*) &(vpos[i].z), sizeof(vpos[i].z));
	}

	int count;
	int fps = (int)finds.size();
	f.write((char *)&vps, sizeof(vps));

	for (size_t i = 0; i < finds.size(); ++i) {
		count = (int)finds[i].size();
		f.write((char *)&count, sizeof(count));
		for (int j = 0; j < count; ++j) {
			f.write((char *) &(finds[i][j]), sizeof(finds[i][j]));
		}
	}

	if (f.bad())
		perror("error while reading file");
}

void loadobj(const std::string &filename, std::vector<vec3f> &vpos, std::vector<std::vector<ivert>> &finds) {
	std::string binout = filename + "b";

	if (!exists_test(binout)) {
		loadObjStrings(filename, finds, vpos);
		saveBinObj(binout, finds, vpos);
	}
	else {
		loadBinObj(binout, finds, vpos);
	}
}



int main() {
	vec3f glob_offset;
	float glob_scale;

	std::vector<vec3f> vposA, vposB;
	std::vector<std::vector<ivert>> findsA, findsB;

	loadobj("..\\..\\sphereA.obj", vposA, findsA);
	loadobj("..\\..\\sphereB.obj", vposB, findsB);

	auto meshA = new Mesh(vposA, findsA);
	auto meshB = new Mesh(vposB, findsB);

	std::vector<Mesh*> meshes;
	meshes.push_back(meshA);
	meshes.push_back(meshB);
	MeshMatch::init_meshes(meshes, glob_offset, glob_scale);

	auto mm = new MeshMatch(meshA, meshB);
	mm->algorithm();
	size_t i;
	return 0;
	
	/*
	std::cout << "\n\nVerts: 0 match 1\n";
	for (i = 0; i < mm->vm01.size(); ++i) {
		int &x = mm->vm01[i];
		if (x == -1) {
			std::cout << i << " ";
		}
	}

	std::cout << "\n\nVerts: 1 match 0\n";
	for (i = 0; i < mm->vm10.size(); ++i) {
		int &x = mm->vm10[i];
		if (x == -1) {
			std::cout << i << " ";
		}
	}
	*/
	
	std::cout << "\n\nFaces: 0 match 1\n";
	for (i = 0; i < mm->fm01.size(); ++i) {
		int &x = mm->fm01[i];
		std::cout << x << ",";
	}

	std::cout << "\n\nFaces: 1 match 0\n";
	for (i = 0; i < mm->fm10.size(); ++i) {
		int &x = mm->fm10[i];
		std::cout << x << ",";
	}

	delete meshA, meshB, mm;
	return 0;
}
