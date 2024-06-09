// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    SparseMatrix<double> invM(M.rows(), M.cols());

    for (int i = 0; i < M.rows(); i++) {
        tripletList.push_back(T(i, i, 1.0 / M.coeff(i, i)));
    }

    invM.setFromTriplets(tripletList.begin(), tripletList.end());
    SparseMatrix<double> laplaceMat = geometry->laplaceMatrix();
    SparseMatrix<double> identity = identityMatrix<double>(laplaceMat.rows());
    return identity + h * invM * laplaceMat;
}

/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {
    SparseMatrix<double> flowOperator = buildFlowOperator(geometry->massMatrix(), h);

    Vector<double> vecX(mesh->nVertices());
    Vector<double> vecY(mesh->nVertices());
    Vector<double> vecZ(mesh->nVertices());

    for (Vertex v : mesh->vertices()) {
        vecX[v.getIndex()] = geometry->inputVertexPositions[v].x;
        vecY[v.getIndex()] = geometry->inputVertexPositions[v].y;
        vecZ[v.getIndex()] = geometry->inputVertexPositions[v].z;
    }

    Vector<double> dX = solveSquare(flowOperator, vecX);
    Vector<double> dY = solveSquare(flowOperator, vecY);
    Vector<double> dZ = solveSquare(flowOperator, vecZ);

    for (Vertex v : mesh->vertices()) {
        geometry->inputVertexPositions[v].x = dX[v.getIndex()];
        geometry->inputVertexPositions[v].y = dY[v.getIndex()];
        geometry->inputVertexPositions[v].z = dZ[v.getIndex()];
    }
}