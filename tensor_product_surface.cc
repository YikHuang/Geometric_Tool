#include "tensor_product_surface.hh"
#include "spline.hh"
#include <array>

template<size_t Dimension>
Eigen::Matrix<float, Dimension + 1, 1> toHomogeneousCoordinates(const Eigen::Matrix<float, Dimension, 1> &x, float w) {
    Eigen::Matrix<float, Dimension + 1, 1> result;
    // TODO 1.1: convert rational B-spline control point data (x, w) into a single
    // vector [w x; w] in homogeneous coordinates.

    for (int i = 0; i < x.rows(); i++) {
        result(i) = x(i) * w;
    }
    result(Dimension) = w;

    return result;
}

template<size_t Dimension>
Eigen::Matrix<float, Dimension, 1> fromHomogeneousCoordinates(const Eigen::Matrix<float, Dimension + 1, 1> &x) {
    // TODO 1.1: convert point `x` expressed in homogeneous coordinates into
    // its Euclidean coordinates representation.

    Eigen::Matrix<float, Dimension, 1> euclideanPt;
    
    for (int i = 0; i < Dimension; i++) {
        euclideanPt(i) = x(i) / x(Dimension);
    }

    return euclideanPt;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::evalPtHomogCoord(double u, double v) const -> HVec {
    size_t I_u = u_spline.findSegmentContaining(u);
    size_t I_v = v_spline.findSegmentContaining(v);
    size_t L_u = u_spline.numControlPts() - 1;
    size_t L_v = v_spline.numControlPts() - 1;
    size_t n_u = u_spline.degree();
    size_t n_v = v_spline.degree();

    // TODO 1.1: evaluate the NURBS surface point's homogeneous coordinates xtilde(u, v)

    for (int i = I_u - (n_u - 1); i <= I_u + 1; i++) {
        for (int j = I_v - (n_v - 1); j <= I_v + 1; j++) {
            HVec vSplineHomoPt = toHomogeneousCoordinates<Dimension>(controlPts(i, j), weights(i, j));
            v_spline.setControlPt(j, vSplineHomoPt);
        }

        HVec uSplineHomoPt = v_spline.evalPt(SplineEvalMethod::DE_BOOR, v);
        u_spline.setControlPt(i, uSplineHomoPt);
    }

    HVec xtilde = u_spline.evalPt(SplineEvalMethod::DE_BOOR, u);
    return xtilde;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::evalPt(double u, double v) const -> Vec {

    // TODO 1.1: evaluate the NURBS surface point's Euclidean coordinates x(u, v)

    HVec xtilde = evalPtHomogCoord(u, v);
    Vec euclideanX = fromHomogeneousCoordinates<Dimension>(xtilde);

    return euclideanX;
}

template<size_t Dimension>
void TensorProductSurface_T<Dimension>::eval(size_t resolution, MX2f &U, MXDf &V, Eigen::MatrixX3i &F) const {
    const size_t numCellsU = resolution - 1;
    const size_t numCellsV = resolution - 1;
    const size_t   numVtxU = resolution;
    const size_t   numVtxV = resolution;

    const auto &cellIdx = [=](size_t a, size_t b) { return a * numCellsV + b; };
    const auto & vtxIdx = [=](size_t a, size_t b) { return a *   numVtxV + b; };

    float u_min = u_spline.domainStart();
    float u_max = u_spline.domainEnd();
    float v_min = v_spline.domainStart();
    float v_max = v_spline.domainEnd();

    // TODO 1.2: sample the NURBS surface on a triangulated regular grid in the
    // UV domain, obtaining triangle mesh (V, F). Also record the UV
    // coordinates (parameter domain) of each vertex in `U`.
    
    std::vector<Vec> points;
    std::vector<Eigen::Vector3i> faces;
    std::vector<Eigen::Vector2f> uvSets;

    // Generate vertices and uvSets
    for (size_t i = 0; i < numVtxU; i++) {
        for (size_t j = 0; j < numVtxV; j++) {

            double alpha = (double)i / numCellsU;
            double uNow = (1 - alpha) * u_min + alpha * u_max;

            double beta = (double)j / numCellsV;
            double vNow = (1 - beta) * v_min + beta * v_max;

            Vec x = evalPt(uNow, vNow);
            points.push_back(x);
            uvSets.emplace_back(uNow, vNow);
        }
    }

    // Faces index assigned
    for (int i = 0; i < numCellsU; i++) {
        for (int j = 0; j < numCellsV; j++) {
            faces.emplace_back(vtxIdx(i, j), vtxIdx(i + 1, j), vtxIdx(i + 1, j + 1));
            faces.emplace_back(vtxIdx(i, j), vtxIdx(i + 1, j + 1), vtxIdx(i, j + 1));
        }
    }


    U.resize(uvSets.size(), 2);
    V.resize(points.size(), 3);
    F.resize(faces.size(), 3);


    for (size_t i = 0; i < points.size(); ++i) V.row(i) = points[i];
    for (size_t i = 0; i < faces.size(); ++i) F.row(i) = faces[i];
    for (size_t i = 0; i < uvSets.size(); ++i) U.row(i) = uvSets[i];
}

template<size_t Dimension>
void TensorProductSurface_T<Dimension>::insertVKnot(double vbar) {
    // TODO 1.3: insert the knot `vbar` into each of the "v curves", creating a
    // new column of data points.

    size_t L_u = u_spline.numControlPts() - 1;
    size_t L_v = v_spline.numControlPts() - 1;

    ControlPtGrid newControlPts = controlPts;
    WeightGrid newWeights = weights;

    HSpline vSplineCopy = HSpline(v_spline.degree(), v_spline.knots, v_spline.controlPts);

    newControlPts.conservativeResize(newControlPts.rows(), newControlPts.cols() + 1);
    newWeights.conservativeResize(newWeights.rows(), newWeights.cols() + 1);


    for (int i = 0; i <= L_u; i++) {
        for (int j = 0; j <= L_v; j++) {
            HVec vSplineHomoPt = toHomogeneousCoordinates<Dimension>(controlPts(i, j), weights(i, j));
            v_spline.setControlPt(j, vSplineHomoPt);
        }

        vSplineCopy = HSpline(v_spline.degree(), v_spline.knots, v_spline.controlPts);
        vSplineCopy.insertKnot(vbar);

        // Update control points and weights to newControlPts and newWeights        
        for (int j = 0; j < vSplineCopy.getControlPts().rows(); j++) {
            Vec newEuclideanControlPt = fromHomogeneousCoordinates<Dimension>(vSplineCopy.controlPt(j));
            newControlPts(i, j) = newEuclideanControlPt;

            float w = vSplineCopy.controlPt(j)(Dimension);
            newWeights(i, j) = w;
        }
    }

    controlPts = newControlPts;
    weights = newWeights;
    v_spline = vSplineCopy;

}

template<size_t Dimension>
void TensorProductSurface_T<Dimension>::insertUKnot(double ubar) {
    // TODO 1.3: insert the knot `ubar` into each of the "u curves", creating a
    // new row of data points.
    size_t L_u = u_spline.numControlPts() - 1;
    size_t L_v = v_spline.numControlPts() - 1;

    ControlPtGrid newControlPts = controlPts;
    WeightGrid newWeights = weights;

    HSpline uSplineCopy = HSpline(u_spline.degree(), u_spline.knots, u_spline.controlPts);

    newControlPts.conservativeResize(newControlPts.rows() + 1, newControlPts.cols());
    newWeights.conservativeResize(newWeights.rows() + 1, newWeights.cols());


    for (int i = 0; i <= L_v; i++) {
        for (int j = 0; j <= L_u; j++) {
            HVec uSplineHomoPt = toHomogeneousCoordinates<Dimension>(controlPts(j, i), weights(j, i));
            u_spline.setControlPt(j, uSplineHomoPt);
        }

        uSplineCopy = HSpline(u_spline.degree(), u_spline.knots, u_spline.controlPts);
        uSplineCopy.insertKnot(ubar);

        // Update control points and weights to newControlPts and newWeights        
        for (int j = 0; j < uSplineCopy.getControlPts().rows(); j++) {
            Vec newEuclideanControlPt = fromHomogeneousCoordinates<Dimension>(uSplineCopy.controlPt(j));
            newControlPts(j, i) = newEuclideanControlPt;

            float w = uSplineCopy.controlPt(j)(Dimension);
            newWeights(j, i) = w;
        }
    }

    controlPts = newControlPts;
    weights = newWeights;
    u_spline = uSplineCopy;
}

template<size_t Dimension>
std::shared_ptr<TensorProductSurface_T<Dimension>> ruled_surface(const NURBS<Dimension> &spline_1, const NURBS<Dimension> &spline_2) {
    auto result_ptr = std::make_shared<TensorProductSurface_T<Dimension>>();
    auto &result = *result_ptr;

    const size_t n = spline_1.degree();
    if (n != spline_2.degree()) throw std::runtime_error("Splines must be of the same degree");

    result.controlPts.resize(spline_1.numControlPts(), 2);
    result.weights   .resize(spline_1.numControlPts(), 2);
    result.getControlPtColFlattened(0) = spline_1.getControlPts();
    result.getControlPtColFlattened(1) = spline_2.getControlPts();
    result.weights.col(0) = spline_1.weights();
    result.weights.col(1) = spline_2.weights();
    result.u_spline = spline_1.homogeneousSpline();

    // Convert v_spline into a simple spline/Bézier curve with just two control
    // points and knots.
    result.v_spline.degree() = 1;
    result.v_spline.controlPts.conservativeResize(2, Dimension + 1);
    result.v_spline.knots.setLinSpaced(2, 0, 1);

    return result_ptr;
}

template<size_t Dimension>
std::shared_ptr<TensorProductSurface_T<Dimension>> surface_of_revolution(const NURBS<Dimension> &spline, size_t numSegments) {
    NURBS<Dimension> circ = NURBSCircle<Dimension>(1.0, numSegments);
    auto result_ptr = std::make_shared<TensorProductSurface_T<Dimension>>();
    auto &result = *result_ptr;

    result.controlPts.resize(circ.numControlPts(), spline.numControlPts());
    result.weights   .resize(circ.numControlPts(), spline.numControlPts());

    for (size_t i = 0; i < circ.numControlPts(); ++i) {
        for (size_t j = 0; j < spline.numControlPts(); ++j) {
            float r = spline.controlPt(j)[0];
            // Use the circ control point to rotate (and scale) around the y axis.
            result.controlPts(i, j) << r * circ  .controlPt(i)[1],
                                           spline.controlPt(j)[1],
                                       r * circ  .controlPt(i)[0];
            result.weights(i, j) = circ.weight(i) * spline.weight(j);
        }
    }

    result.u_spline = circ.homogeneousSpline();
    result.v_spline = spline.homogeneousSpline();

    return result_ptr;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::partial_xtilde(double u, double v, size_t diff_u, size_t diff_v) const -> HVec {
    // TODO 2.2: compute the partial derivative of the homogeneous coordinate
    // vector "xtilde" `diff_u` times with respect to `u` and `diff_v` times
    // with respect to `v` (applying `diff_u + diff_v` derivatives in total).

    size_t I_u = u_spline.findSegmentContaining(u);
    size_t I_v = v_spline.findSegmentContaining(v);
    size_t L_u = u_spline.numControlPts() - 1;
    size_t L_v = v_spline.numControlPts() - 1;
    size_t n_u = u_spline.degree();
    size_t n_v = v_spline.degree();

    for (int i = I_u - (n_u - 1); i <= I_u + 1; i++) {
        for (int j = I_v - (n_v - 1); j <= I_v + 1; j++) {
            HVec vSplineHomoPt = toHomogeneousCoordinates<Dimension>(controlPts(i, j), weights(i, j));
            v_spline.setControlPt(j, vSplineHomoPt);
        }

        // Get Derivatives v spline
        HSpline dvSpline = v_spline;
        for (int j = 0; j < diff_v; j++) {
            dvSpline = derivative(dvSpline);
        }

        HVec uSplineHomoPt = dvSpline.evalPt(SplineEvalMethod::BASIS, v);
        u_spline.setControlPt(i, uSplineHomoPt);
    }

    // Get Derivatives u spline
    HSpline duSpline = u_spline;
    for (int i = 0; i < diff_u; i++) {
        duSpline = derivative(duSpline);
    }

    HVec xtilde = duSpline.evalPt(SplineEvalMethod::BASIS, u);
    return xtilde;

    //return HVec::Zero();
}

template<size_t Dimension>
Eigen::Matrix<float, Dimension, 2> TensorProductSurface_T<Dimension>::jacobian(double u, double v) const {
    HVec xtilde = evalPtHomogCoord(u, v);

    Eigen::Matrix<float, Dimension + 1, 2> J_xtilde;
    // TODO 2.3: Compute the Jacobian of `xtilde` in `J_xtilde`.
    J_xtilde.col(0) = partial_xtilde(u, v, 1, 0);
    J_xtilde.col(1) = partial_xtilde(u, v, 0, 1);


    //HVec xtilde = evalPtHomogCoord(u, v);

    return J_xtilde.template topRows<Dimension>() / xtilde[Dimension]
            - xtilde.template head<Dimension>() * (J_xtilde.template bottomRows<1>() / std::pow(xtilde[Dimension], 2));
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::normal(double u, double v) const -> Vec {
    // TODO 2.3: Compute the unit surface normal at (u, v)
    
    Eigen::Matrix<float, Dimension, 2> J = jacobian(u, v);

    Vec xCrossProduct;
    xCrossProduct(0) = J(1, 0) * J(2, 1) - J(2, 0) * J(1, 1);
    xCrossProduct(1) = J(2, 0) * J(0, 1) - J(0, 0) * J(2, 1);
    xCrossProduct(2) = J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1);

    double xCrossProductNorm = xCrossProduct.norm();
    return xCrossProduct / xCrossProductNorm;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::firstFundamentalForm(double u, double v) const -> M2d {
    // TODO 2.3: Calculate the first fundamental form at (u, v)
    Eigen::Matrix<float, Dimension, 2> J = jacobian(u, v);
    Eigen::Matrix<float, 2, Dimension> JT = J.transpose();

    return JT * J;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::secondFundamentalForm(double u, double v) const -> M2d {
    // Quantities needed in the chain rule terms for the perspective divide.
    HVec xtilde = evalPtHomogCoord(u, v);

    Eigen::Matrix<float, Dimension + 1, 2> J_xtilde;
    // TODO 2.4: copy your J_xtilde computation from 2.3 here.
    J_xtilde.col(0) = partial_xtilde(u, v, 1, 0);
    J_xtilde.col(1) = partial_xtilde(u, v, 0, 1);


    std::array<std::array<HVec, 2>, 2> d2xtilde;
    // TODO 2.4: calculate the "vector-valued matrix" of second partial derivatives of `xtilde`
    d2xtilde[0][0] = partial_xtilde(u, v, 2, 0);
    d2xtilde[0][1] = partial_xtilde(u, v, 1, 1);
    d2xtilde[1][0] = partial_xtilde(u, v, 1, 1);
    d2xtilde[1][1] = partial_xtilde(u, v, 0, 2);


    // The following applies the chain rule to obtain the matrix of second partial derivatives
    // of Euclidean surface point `x`.
    std::array<std::array<Vec, 2>, 2> d2x;
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            // Convert derivatives of the surface point in homogeneous
            // coordinates into derivatives of the Euclidean coordinates.
            d2x[i][j] = d2xtilde[i][j].template head<Dimension>() / xtilde[Dimension]
                      -  (J_xtilde(Dimension, i) * J_xtilde.col(j).template head<Dimension>()
                        + J_xtilde(Dimension, j) * J_xtilde.col(i).template head<Dimension>()
                        + xtilde.template head<Dimension>() * d2xtilde[i][j][Dimension]) / std::pow(xtilde[Dimension], 2)
                      + xtilde.template head<Dimension>() * (2 * (J_xtilde(Dimension, i) * J_xtilde(Dimension, j)) / std::pow(xtilde[Dimension], 3))
                    ;
        }
    }
    d2x[0][1] = d2x[1][0]; // The loop above only computes the lower triangle (exploiting symmetry)

    // TODO 2.4: evaluate the second fundamental form.
    
    Vec xNormal = normal(u, v);
    M2d II = M2d::Zero();

    for (int i = 0; i < d2x.size(); i++) {
        for (int j = 0; j < d2x[i].size(); j++) {
            II(i, j) = xNormal.dot(d2x[i][j]);
        }
    }
    
    return II;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::principalCurvatures(double u, double v) const -> V2d {
    // TODO: 2.4: solve a generalized eigenvalue problem for the principal curvatures.
    M2d II = secondFundamentalForm(u, v);
    M2d I = firstFundamentalForm(u, v);

    Eigen::GeneralizedSelfAdjointEigenSolver<M2d> es(II, I, Eigen::DecompositionOptions::Ax_lBx);

    V2d curvatures;
    curvatures(0, 0) = es.eigenvalues()(1);
    curvatures(1, 0) = es.eigenvalues()(0);

    return curvatures;
}

template<size_t Dimension>
auto TensorProductSurface_T<Dimension>::gaussianAndMeanCurvature(double u, double v) const -> V2d {
    // TODO: 2.4: evaluate (K, H) using the principal curvatures.
    V2d curvatures = principalCurvatures(u, v);

    double gaussian = curvatures(0, 0) * curvatures(1, 0);
    double mean = (curvatures(0, 0) + curvatures(1, 0)) / 2;

    V2d gaussianAndMeanMatrix;
    gaussianAndMeanMatrix(0, 0) = gaussian;
    gaussianAndMeanMatrix(1, 0) = mean;

    return gaussianAndMeanMatrix;
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation of all class and function templates and functions.
////////////////////////////////////////////////////////////////////////////////
template struct TensorProductSurface_T<3>;

template std::shared_ptr<TensorProductSurface> ruled_surface<3>(const NURBS<3> &spline_1, const NURBS<3> &spline_2);
template std::shared_ptr<TensorProductSurface> surface_of_revolution<3>(const NURBS<3> &spline, size_t numSegments);
