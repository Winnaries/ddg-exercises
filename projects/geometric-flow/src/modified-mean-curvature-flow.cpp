// Implement member functions for ModifiedMeanCurvatureFlow class.
#include "modified-mean-curvature-flow.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ModifiedMeanCurvatureFlow::ModifiedMeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry, and A (Laplace matrix)
    mesh = inputMesh;
    geometry = inputGeo;

    this->A = geometry->laplaceMatrix();
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> ModifiedMeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    SparseMatrix<double> invM(M.rows(), M.cols());

    for (int i = 0; i < M.rows(); i++) {
        tripletList.push_back(T(i, i, 1.0 / M.coeff(i, i)));
    }

    invM.setFromTriplets(tripletList.begin(), tripletList.end());
    SparseMatrix<double> identity = identityMatrix<double>(A.rows());
    return identity + h * invM * A;
}