#include "bezier.hh"
#include <iostream>

Eigen::Vector3f BezierCurve::evalPt(EvalMethod method, double t) const {
    if (method == EvalMethod::DE_CASTELJAU) {
        // TODO: Copy your implementation of HW2 Problem 1.1
        return Eigen::Vector3f::Zero();
    }
    if (method == EvalMethod::BERNSTEIN) {
        // TODO: Copy your implementation of HW2 Problem 2.1
        return Eigen::Vector3f::Zero();
    }
    if (method == EvalMethod::HORNER) {
        // TODO: Copy your implementation of HW2 Problem 2.2
        return Eigen::Vector3f::Zero();
    }
    throw std::runtime_error("Unknown evaluation method");
}

void BezierCurve::eval(size_t resolution, EvalMethod method, std::vector<Eigen::Vector3f> &result) const {
    // TODO: Copy your implementation of HW2 Problem 1.2
    result.assign(resolution, Eigen::Vector3f::Zero());
}

void BezierCurve::visualizeDeCasteljau(double t, Eigen::MatrixX3f &V, Eigen::MatrixX2i &E) const {
    // TODO: Copy your implementation of HW2 Problem 1.3
    size_t n = degree();
    V.resize(0, 3);
    E.resize(0, 2);
}

// Get the tangent from the last evaluation
Eigen::Vector3f BezierCurve::probeCurve(double t, Eigen::Vector3f &tangent, Eigen::Vector3f &normal) const {
    // TODO: Copy your implementation of HW2 Problem 4.1
    tangent << 1, 0, 0;
    normal  << 0, 1, 0;
    Eigen::Vector3f position(0, 0, 0);
    return position;
}

void BezierCurve::setDegree(size_t new_degree) {
    // TODO: Copy your implementation of HW2 Problem 3.1 to be used for degree elevation.

    // TODO: HW3 Problem 2.1
    // Implement degree reduction by repeatedly calling `reduceDegree()`
    // while `degree() < new_degree`.

    if (new_degree < 1) throw std::runtime_error("Invalid degree");
    size_t npts = new_degree + 1;

    const auto &oldPts = getControlPts();

    auto newPts = oldPts;

    setControlPts(newPts);
}

void BezierCurve::subdivide(double t, BezierCurve &c1, BezierCurve &c2) const {
    // TODO: Copy your implementation of HW2 Problem 3.2
    c1 = *this;
    c2 = *this;
}

BoundingBox BezierCurve::boundingBox() const {
    BoundingBox bb;
    // TODO: Copy your implementation of HW2 Problem 4.2
    return bb;
}

bool BoundingBox::overlaps(const BoundingBox &b) const {
    // TODO: Copy your implementation of HW2 Problem 4.2
    // Return whether this bounding box overlaps the other bounding box `b`.
    return false;
}

void getIntersections(const BezierCurve &c1, const BezierCurve &c2, std::vector<Eigen::Vector3f> &result, float tol) {
    // TODO: Copy your implementation of HW2 Problem 4.2
    // *Append* approximations of *all* intersections between curves `c1` and
    // `c2` into the result collection `result`.
}

void mergeDuplicateIntersections(std::vector<Eigen::Vector3f> &ipoints, float tol) {
    // TODO: Copy your implementation of HW2 Bonus (if attempted).
}

void BezierCurve::reduceDegree() {
    // TODO: HW3 Problem 2.1
    // Implement one step of degree reduction, computing the control points of
    // the one-degree-lower Bézier curve whose *elevation* best approximates
    // the current control points.
    Eigen::MatrixXf X_tilde = pointsToMatrixRows(b[0]);
    setControlPts(pointsFromMatrixRows(X_tilde));
}

void BezierCurve::drag(float t, const Eigen::Vector3f &p) {
    // TODO: HW3 Problem 2.2
    // Calculate and apply the minimum-norm displacement to the current control points
    // such that the curve passes through `p` at parameter `t` (i.e., `c(t) = p`).
}

void BezierCurve::approximate(const Eigen::VectorXd &t_values, const std::vector<Eigen::Vector3f> &pts, double smoothingWeight) {
    size_t m = t_values.size(), n = degree();
    if (pts.size() != m) throw std::runtime_error("Time and point size mismatch");

    // TODO: HW3 Problem 3.1
    // Solve for the Bézier control point positions that represent the best-fit
    // polynomial curve for the data {(t_values[i], pts[i])}.
    Eigen::MatrixXf B = pointsToMatrixRows(b[0]);
    setControlPts(pointsFromMatrixRows(B));

    // TODO: HW3 Problem 3.2
    // Incorporate a smoothing regularization term into the least-squares approximation
    // using the weight `smoothingWeight`.
}
