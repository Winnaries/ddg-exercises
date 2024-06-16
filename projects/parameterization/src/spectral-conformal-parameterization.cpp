// Implement member functions for SpectralConformalParameterization class.
#include "spectral-conformal-parameterization.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
SpectralConformalParameterization::SpectralConformalParameterization(ManifoldSurfaceMesh* inputMesh,
                                                                     VertexPositionGeometry* inputGeo) {

    this->mesh = inputMesh;
    this->geometry = inputGeo;
}

/*
 * Builds the complex conformal energy matrix EC = ED - A.
 *
 * Input:
 * Returns: A complex sparse matrix representing the conformal energy
 */
SparseMatrix<std::complex<double>> SpectralConformalParameterization::buildConformalEnergy() const {
    typedef Eigen::Triplet<std::complex<double>> T; 
    std::vector<T> tripletList; 
    tripletList.reserve(mesh->nExteriorHalfedges() * 2); 

    for (BoundaryLoop bl : mesh->boundaryLoops()) {
        for (Halfedge he : bl.adjacentHalfedges()) {
            size_t vertexI = he.tailVertex().getIndex(); 
            size_t vertexJ = he.tipVertex().getIndex(); 
            tripletList.push_back(T(vertexI, vertexJ, std::complex<double>(0, 0.25))); 
            tripletList.push_back(T(vertexJ, vertexI, std::complex<double>(0, -0.25))); 
        }
    }

    SparseMatrix<std::complex<double>> matL = geometry->complexLaplaceMatrix() * 0.5; 
    SparseMatrix<std::complex<double>> matA(mesh->nVertices(), mesh->nVertices()); 
    matA.setFromTriplets(tripletList.begin(), tripletList.end());  
    return matL - matA; 
}


/*
 * Flattens the input surface mesh with 1 or more boundaries conformally.
 *
 * Input:
 * Returns: A MeshData container mapping each vertex to a vector of planar coordinates.
 */
VertexData<Vector2> SpectralConformalParameterization::flatten() const {
    Vector<std::complex<double>> z = solveInversePowerMethod(buildConformalEnergy()); 

    VertexData<Vector2> result(*mesh);
    for (Vertex v : mesh->vertices()) {
        size_t index = v.getIndex(); 
        result[v] = Vector2::fromComplex(z[index]); 
    }

    return result; 
}