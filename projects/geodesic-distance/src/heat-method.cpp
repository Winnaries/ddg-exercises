// Implement member functions HeatMethod class.
#include "heat-method.h"
#include "geometrycentral/numerical/linear_solvers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HeatMethod::HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {

    this->mesh = surfaceMesh;
    this->geometry = geo;

    double t = geo->meanEdgeLength() * geo->meanEdgeLength();
    this->A = geo->laplaceMatrix(); 
    this->F = geo->massMatrix() + (SparseMatrix<double>)(t * this->A); 
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {
    FaceData<Vector3> faceData(*mesh);

    for (Face f : mesh->faces()) {
        Vector3 faceGradient = Vector3::constant(0.0); 
        Vector3 faceNormal = geometry->faceNormal(f); 
        double faceArea = geometry->faceArea(f);

        for (Halfedge he : f.adjacentHalfedges()) {
            Vertex v = he.next().next().vertex(); 
            Vector3 edgeVector = geometry->halfedgeVector(he); 
            faceGradient += u[v.getIndex()] * edgeVector; 
        }

        faceData[f] = -(cross(faceNormal, faceGradient) / faceArea / 2.).normalize();
    }

    return faceData; 
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {
    Vector<double> gradX = Vector<double>::Zero(mesh->nVertices());

    for (Vertex v : mesh->vertices()) {
        for (Halfedge heA : v.outgoingHalfedges()) {
            Face face = heA.face(); 
            Halfedge heB = heA.next().next(); 
            Vector3 eA = geometry->halfedgeVector(heA);
            Vector3 eB = geometry->halfedgeVector(heB.twin());
            Vector3 fGrad = X[face]; 
            gradX[v.getIndex()] += .5 * geometry->cotan(heA) * dot(fGrad, eA);
            gradX[v.getIndex()] += .5 * geometry->cotan(heB) * dot(fGrad, eB); 
        }
    }

    return gradX; 
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {
    SparseMatrix<double> F = this->F; 
    SparseMatrix<double> A = this->A; 
    Vector<double> ut = solvePositiveDefinite(F, delta); 
    Vector<double> div = computeDivergence(computeVectorField(ut)); 
    Vector<double> phi = -solvePositiveDefinite(A, div);

    // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    this->subtractMinimumDistance(phi);

    return phi;
}