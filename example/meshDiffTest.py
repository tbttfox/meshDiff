import os, sys
from itertools import count
from meshDiff import meshDiff

BASE = os.path.dirname(__file__)

def loadObj(path):
	''' A simple wavefront .obj file loader '''
	vertices = []
	faces = []
	with open(path, 'r') as inFile:
		lines = inFile.readlines()

	for line in lines:
		sp = line.split()
		if sp == []:
			pass
		elif sp[0] == "v":
			v = [float(i) for i in sp[1:4]]
			vertices.append(v)

		elif sp[0] == "f":
			face = []
			for s in sp[1:]:
				face.append(int(s.split('/', 1)[0]) - 1)
			faces.append(face)

	return vertices, faces

def main():
	origPath = os.path.join(BASE, "Useful", "sphereA.obj")
	twistedPath = os.path.join(BASE, "Useful", "sphereB.obj")

	meshAPts, meshAFaces = loadObj(origPath)
	print "Loaded Origin"
	meshBPts, meshBFaces = loadObj(twistedPath)
	print "Loaded Twisted"

	out = meshDiff(meshAPts, meshAFaces, meshBPts, meshBFaces)
	#print out
	print zip(count(0, 1), out[0][0])
	print "\n"
	print zip(count(0, 1), out[0][1])
	print "\n\n\n"
	print zip(count(0, 1), out[1][0])
	print "\n"
	print zip(count(0, 1), out[1][1])
	print "Done"



if __name__ == "__main__":
	main()



