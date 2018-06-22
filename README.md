# meshDiff

Using the paper and source code from Jon Denning's MeshGit (SIGGRAPH 2013, http://cse.taylor.edu/~jdenning/research.html), I did a little optimization, and built a python module for computing the the minimal change set to convert between two meshes.

The current build system is only tested for Python2.7 on Windows, but it shouldn't be difficult to make it work elsewhere. There is wrapper code for Python3.6 too, if anybody wants that.

# Usage
See example/meshDiffTest.py for a working example

```mAVertMatch, mBVertMatch, mAFaceMatch, mBFaceMatch = meshDiff(meshAPts, meshAFaces, meshBPts, meshBFaces)```

In each match, each item contains the index in other mesh that it matches to, or -1 if it is unmatched.
