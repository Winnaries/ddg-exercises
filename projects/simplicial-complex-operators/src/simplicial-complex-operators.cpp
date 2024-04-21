// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    // size_t idx = 0;
    // for (Vertex v : mesh->vertices()) {
    //     idx = geometry->vertexIndices[v];
    // }

    // for (Edge e : mesh->edges()) {
    //     idx = geometry->edgeIndices[e];
    // }

    // for (Face f : mesh->faces()) {
    //     idx = geometry->faceIndices[f];
    // }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    // for (Vertex v : mesh->vertices()) {
    //     idx = v.getIndex(); // == geometry->vertexIndices[v])
    // }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> tripletList;
    tripletList.reserve(2 * mesh->nEdges());

    for (Edge e : mesh->edges()) {
        Vertex firstVertex = e.firstVertex();
        Vertex secondVertex = e.secondVertex();

        size_t edgeIndex = e.getIndex();
        size_t firstIndex = firstVertex.getIndex();
        size_t secondIndex = secondVertex.getIndex();

        tripletList.push_back(T(edgeIndex, firstIndex, 1));
        tripletList.push_back(T(edgeIndex, secondIndex, 1));
    }

    SparseMatrix<size_t> mat(mesh->nEdges(), mesh->nVertices());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
    typedef Eigen::Triplet<size_t> T;
    std::vector<T> tripletList;
    tripletList.reserve(3 * mesh->nFaces());

    for (Face f : mesh->faces()) {
        Halfedge he = f.halfedge();
        for (size_t i = 0; i < 3; i++) {
            size_t faceIndex = f.getIndex();
            size_t edgeIndex = he.edge().getIndex();
            tripletList.push_back(T(faceIndex, edgeIndex, 1));
            he = he.next();
        }
    }

    SparseMatrix<size_t> mat(mesh->nFaces(), mesh->nEdges());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat; // placeholder
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {
    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nVertices());

    for (auto index : subset.vertices) {
        vec[index] = 1;
    }

    return vec;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {
    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nEdges());

    for (auto index : subset.edges) {
        vec[index] = 1;
    }

    return vec;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {
    Vector<size_t> vec = Vector<size_t>::Zero(mesh->nFaces());

    for (auto index : subset.faces) {
        vec[index] = 1;
    }

    return vec;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    Vector<size_t> vertexVector = buildVertexVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset) + A0 * vertexVector;
    Vector<size_t> faceVector = buildFaceVector(subset) + A1 * edgeVector;

    std::set<size_t> newVertices;
    std::set<size_t> newEdges;
    std::set<size_t> newFaces;

    for (size_t i = 0; i < mesh->nVertices(); i++) {
        if (vertexVector[i] >= 1) {
            newVertices.insert(i);
        }
    }

    for (size_t i = 0; i < mesh->nEdges(); i++) {
        if (edgeVector[i] >= 1) {
            newEdges.insert(i);
        }
    }

    for (size_t i = 0; i < mesh->nFaces(); i++) {
        if (faceVector[i] >= 1) {
            newFaces.insert(i);
        }
    }

    return MeshSubset(newVertices, newEdges, newFaces);
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    Vector<size_t> faceVector = buildFaceVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset) + A1.transpose() * faceVector;
    Vector<size_t> vertexVector = buildVertexVector(subset) + A0.transpose() * edgeVector;

    std::set<size_t> newVertices;
    std::set<size_t> newEdges;
    std::set<size_t> newFaces;

    for (size_t i = 0; i < mesh->nVertices(); i++) {
        if (vertexVector[i] >= 1) {
            newVertices.insert(i);
        }
    }

    for (size_t i = 0; i < mesh->nEdges(); i++) {
        if (edgeVector[i] >= 1) {
            newEdges.insert(i);
        }
    }

    for (size_t i = 0; i < mesh->nFaces(); i++) {
        if (faceVector[i] >= 1) {
            newFaces.insert(i);
        }
    }

    return MeshSubset(newVertices, newEdges, newFaces);
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    MeshSubset stcl = star(closure(subset));
    MeshSubset clst = closure(star(subset));

    clst.deleteSubset(stcl);

    return clst;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
    return subset.equals(closure(subset));
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    if (!isComplex(subset)) {
        return -1;
    }

    Vector<size_t> faceA = buildFaceVector(subset);
    Vector<size_t> edgeA = buildEdgeVector(subset);
    Vector<size_t> vertexA = buildVertexVector(subset);

    Vector<size_t> edgeB = (faceA.transpose() * A1).cwiseMin(1);
    Vector<size_t> vertexB = (edgeA.transpose() * A0).cwiseMin(1);


    if (faceA.sum() > 0) {
        return edgeA == edgeB && vertexA == vertexB ? 2 : -1;
    }

    if (edgeA.sum() > 0) {
        return vertexA == vertexB ? 1 : -1;
    }

    if (vertexA.sum() > 0) {
        return 0;
    }

    return -1;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    Vector<size_t> edgeVector = A1.transpose() * buildFaceVector(subset);
    Vector<size_t> vertexVector = A0.transpose() * buildEdgeVector(subset);

    MeshSubset bd = MeshSubset();

    for (size_t i = 0; i < mesh->nEdges(); i++) {
        if (edgeVector[i] == 1) {
            bd.addEdge(i);
        }
    }

    for (size_t i = 0; i < mesh->nVertices(); i++) {
        if (vertexVector[i] == 1) {
            bd.addVertex(i);
        }
    }

    return closure(bd); // placeholder
}