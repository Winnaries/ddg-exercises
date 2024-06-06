// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {
    if (he.isInterior()) {
        Vector3 pB = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pC = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pA = inputVertexPositions[he.vertex()];

        Vector3 e1 = pB - pA;
        Vector3 e2 = pC - pA;

        return dot(e1, e2) / norm(cross(e1, e2));
    } else {
        return 0.;
    }
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {
    double A = 0.;
    for (Face f : v.adjacentFaces()) {
        Halfedge he = f.halfedge();
        Vector3 pA = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pB = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pC = inputVertexPositions[he.vertex()];

        double area = 0.5 * norm(cross(pB - pA, pC - pA));
        A += area / 3.;
    }
    return A;
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {
    Halfedge he = c.halfedge();
    Vector3 pA = inputVertexPositions[he.vertex()];
    he = he.next();
    Vector3 pB = inputVertexPositions[he.vertex()];
    he = he.next();
    Vector3 pC = inputVertexPositions[he.vertex()];
    return std::acos(dot((pB - pA).normalize(), (pC - pA).normalize()));
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {
    Vector3 pI = inputVertexPositions[he.vertex()];
    Vector3 pJ = inputVertexPositions[he.twin().vertex()];
    Vector3 pK = inputVertexPositions[he.next().next().vertex()];
    Vector3 pL = inputVertexPositions[he.twin().next().next().vertex()];

    Vector3 eIJ = pJ - pI;
    Vector3 nIJK = cross(pJ - pI, pK - pI).normalize();
    Vector3 nIJL = cross(pL - pI, pJ - pI).normalize();

    return std::atan2(dot(eIJ.normalize(), cross(nIJK, nIJL)), dot(nIJK, nIJL));
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {
    Vector3 normal = {0., 0., 0.};

    for (Face f : v.adjacentFaces()) {
        Halfedge he = f.halfedge();
        Vector3 pA = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pB = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pC = inputVertexPositions[he.vertex()];
        normal += cross(pB - pA, pC - pA).normalize();
    }

    return normal.normalize();
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {
    Vector3 normal = {0., 0., 0.};

    for (Halfedge he : v.outgoingHalfedges()) {
        Corner corner = he.corner();
        Vector3 pA = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pB = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pC = inputVertexPositions[he.vertex()];
        Vector3 n = cross(pB - pA, pC - pA).normalize();
        normal += angle(corner) * n;
    }

    return normal.normalize();
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {
    Vector3 normal = {0., 0., 0.};

    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 pA = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pB = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pC = inputVertexPositions[he.vertex()];
        Vector3 n = cross(pB - pA, pC - pA);
        normal += n / norm2(pB - pA) / norm2(pC - pA);
    }

    return normal.normalize();
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {
    Vector3 normal = {0., 0., 0.};

    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 pA = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pB = inputVertexPositions[he.vertex()];
        he = he.next();
        Vector3 pC = inputVertexPositions[he.vertex()];
        Vector3 n = cross(pB - pA, pC - pA);
        normal += norm(n) * n.normalize();
    }

    return normal.normalize();
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {
    Vector3 normal = {0., 0., 0.};

    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 pI = inputVertexPositions[he.vertex()];
        Vector3 pJ = inputVertexPositions[he.next().vertex()];
        Vector3 eIJ = pJ - pI;
        double psiIJ = dihedralAngle(he);
        normal += psiIJ / norm(eIJ) * eIJ;
    }

    return normal.normalize();
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {
    Vector3 normal = {0., 0., 0.};

    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 pI = inputVertexPositions[he.vertex()];
        Vector3 pJ = inputVertexPositions[he.next().vertex()];
        Vector3 eIJ = pJ - pI;
        double cotanA = cotan(he);
        double cotanB = cotan(he.twin());
        normal += (cotanA + cotanB) * eIJ;
    }

    return normal.normalize();
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {
    double K = 2. * PI;

    for (Halfedge he : v.outgoingHalfedges()) {
        K -= angle(he.corner());
    }

    return K;
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {
    double K = 0.;

    for (Vertex v : mesh.vertices()) {
        K += angleDefect(v);
    }

    return K;
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {
    double H = 0.;

    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 pI = inputVertexPositions[he.vertex()];
        Vector3 pJ = inputVertexPositions[he.next().vertex()];
        Vector3 eIJ = pJ - pI;
        H += dihedralAngle(he) * norm(eIJ) / 2.;
    }

    return H;
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {
    double A = 0.;

    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 pI = inputVertexPositions[he.vertex()];
        Vector3 pJ = inputVertexPositions[he.next().vertex()];
        A += norm2(pJ - pI) * cotan(he);
    }

    for (Halfedge he : v.incomingHalfedges()) {
        Vector3 pI = inputVertexPositions[he.vertex()];
        Vector3 pJ = inputVertexPositions[he.next().vertex()];
        A += norm2(pJ - pI) * cotan(he);
    }

    return A / 8.;
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {
    double dualArea = circumcentricDualArea(v);
    double h = scalarMeanCurvature(v) / dualArea;
    double k = angleDefect(v) / dualArea;

    double k1 = h - std::sqrt(h * h - k);
    double k2 = h + std::sqrt(h * h - k);

    return std::make_pair(k1, k2);
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    // TODO
    return identityMatrix<double>(1); // placeholder
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    // TODO
    return identityMatrix<std::complex<double>>(1); // placeholder
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral