// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh.nVertices());

    for (Vertex v : mesh.vertices()) {
        tripletList.push_back(T(v.getIndex(), v.getIndex(), barycentricDualArea(v)));
    }

    SparseMatrix<double> mat(mesh.nVertices(), mesh.nVertices());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh.nEdges());

    for (Edge e : mesh.edges()) {
        double cotanA = cotan(e.halfedge());
        double cotanB = cotan(e.halfedge().twin());
        tripletList.push_back(T(e.getIndex(), e.getIndex(), (cotanA + cotanB) * .5));
    }

    SparseMatrix<double> mat(mesh.nEdges(), mesh.nEdges());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh.nFaces());

    for (Face f : mesh.faces()) {
        Halfedge he = f.halfedge();
        Vector3 pA = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pB = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pC = inputVertexPositions[he.vertex()];

        double area = 0.5 * norm(cross(pB - pA, pC - pA));
        tripletList.push_back(T(f.getIndex(), f.getIndex(), 1.0 / area));
    }

    SparseMatrix<double> mat(mesh.nFaces(), mesh.nFaces());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh.nEdges() * 2);

    for (Edge e : mesh.edges()) {
        Vertex vA = e.firstVertex();
        Vertex vB = e.secondVertex();
        tripletList.push_back(T(e.getIndex(), vA.getIndex(), -1));
        tripletList.push_back(T(e.getIndex(), vB.getIndex(), 1));
    }

    SparseMatrix<double> mat(mesh.nEdges(), mesh.nVertices());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(mesh.nFaces() * 3);

    for (Face f : mesh.faces()) {
        Halfedge he = f.halfedge();
        for (size_t i = 0; i < 3; i++) {
            tripletList.push_back(T(f.getIndex(), he.edge().getIndex(), he.orientation() ? 1 : -1));
            he = he.next();
        }
    }

    SparseMatrix<double> mat(mesh.nFaces(), mesh.nEdges());
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    return mat;
}

} // namespace surface
} // namespace geometrycentral