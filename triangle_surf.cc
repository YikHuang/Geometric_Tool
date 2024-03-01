#include "triangle_surf.hh"
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>

auto TriangleSurface::laplacian() const -> SpMat {
    SpMat L;
    igl::cotmatrix(V, F, L);
    return L;
}

Eigen::VectorXd TriangleSurface::voronoiAreas() const {
    SpMat M;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    return M.diagonal();
}

Eigen::MatrixXd TriangleSurface::normals() const {
    // TODO 3.1: construct the per-vertex area-weighted normal vector field.
    Eigen::MatrixX3d vNormals = Eigen::MatrixXd::Zero(V.rows(), 3);


    for (int i = 0; i < F.rows(); i++) {
        Eigen::Vector3d vertexA = V.row(F.row(i)(0));
        Eigen::Vector3d vertexB = V.row(F.row(i)(1));
        Eigen::Vector3d vertexC = V.row(F.row(i)(2));

        Eigen::Vector3d nT = (vertexB - vertexA).cross(vertexC - vertexA);
        double AT = nT.norm() / 2;

        Eigen::Vector3d areaWeightedNormal = AT * nT;

        vNormals.row(F.row(i)(0)) += areaWeightedNormal;
        vNormals.row(F.row(i)(1)) += areaWeightedNormal;
        vNormals.row(F.row(i)(2)) += areaWeightedNormal;
    }

    for (int i = 0; i < vNormals.rows(); i++) {
        vNormals.row(i) = vNormals.row(i) / vNormals.row(i).norm();
    }

    return vNormals;
}

Eigen::VectorXd TriangleSurface::gaussianCurvatures() const {
    // TODO 3.2: calculate the vector of per-vertex discrete Gaussian
    // curvatures using the angle defect formula.
    // Be sure to zero out the boundary values using `zeroOutBoundaryValues`
    VXd vAngles = VXd::Zero(V.rows());

    for (int i = 0; i < F.rows(); i++) {
        Eigen::Vector3d vertexA = V.row(F.row(i)(0));
        Eigen::Vector3d vertexB = V.row(F.row(i)(1));
        Eigen::Vector3d vertexC = V.row(F.row(i)(2));

        double angleA = atan2(((vertexB - vertexA).cross((vertexC - vertexA))).norm(), (vertexB - vertexA).dot((vertexC - vertexA)));
        double angleB = atan2(((vertexA - vertexB).cross((vertexC - vertexB))).norm(), (vertexA - vertexB).dot((vertexC - vertexB)));
        double angleC = atan2(((vertexA - vertexC).cross((vertexB - vertexC))).norm(), (vertexA - vertexC).dot((vertexB - vertexC)));
        
        vAngles(F.row(i)(0)) += angleA;
        vAngles(F.row(i)(1)) += angleB;
        vAngles(F.row(i)(2)) += angleC;
    }

    zeroOutBoundaryValues(vAngles);

    Eigen::VectorXd K = Eigen::VectorXd::Zero(V.rows());
    Eigen::VectorXd voronoi = TriangleSurface::voronoiAreas();
    Eigen::VectorXd voronoiReciprocal = Eigen::VectorXd::Zero(voronoi.rows());

    for (int i = 0; i < voronoi.rows(); i++) {
        voronoiReciprocal(i) = 1.0 / voronoi(i);
    }

    for (int i = 0; i < vAngles.rows(); i++) {
        K(i) = voronoiReciprocal(i) * (2 * M_PI - vAngles(i));
    }

    return K;
}

Eigen::VectorXd TriangleSurface::meanCurvatures() const {
    // TODO 3.3: calculate the vector of per-vertex discrete mean
    // curvatures using the discrete Laplace-Beltrami operator.
    // Be sure to zero out the boundary values using `zeroOutBoundaryValues`
    return VXd::Zero(V.rows());
}

Eigen::VectorXd TriangleSurface::kappa_1() const {
    // TODO 3.4: calculate the vector of per-vertex first principal curvatures.
    return VXd::Zero(V.rows());
}

Eigen::VectorXd TriangleSurface::kappa_2() const {
    // TODO 3.4: calculate the vector of per-vertex second principal curvatures.
    return VXd::Zero(V.rows());
}
