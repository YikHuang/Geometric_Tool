#include "spline.hh"
#include <iostream>

template<size_t Dimension>
double SplineBase<Dimension>::bsplineBasisFunction(size_t n, size_t i, double u) const {
    // HW4 TODO 1.3: Evaluate the B-Spline basis function N_i^n(u) using recursion.
    size_t head = 0;
    size_t tail = knots.size() - 1;

    if (n == 0) {
        if (u >= knots[clamp(i - 1, head, tail)] && u < knots[clamp(i, head, tail)]) {
            return 1.0f;
        }
        else {
            return 0.0f;
        }
    }


    double p1 = 0.0;
    double n1 = bsplineBasisFunction(n - 1, i, u);
    if (n1 != 0) {
        double w1 = (u - knots[clamp(i - 1, head, tail)]) / (knots[clamp(n - 1 + i, head, tail)] - knots[clamp(i - 1, head, tail)]);
        p1 = w1 * n1;
    }

    double p2 = 0.0;
    double n2 = bsplineBasisFunction(n - 1, i + 1, u);
    if (n2 != 0) {
        double w2 = (knots[clamp(n + i, head, tail)] - u) / (knots[clamp(n + i, head, tail)] - knots[clamp(i, head, tail)]);
        p2 = w2 * n2;
    }

    return p1 + p2;
}

template<size_t Dimension>
size_t Spline_T<Dimension>::findSegmentContaining(double u) const {
    size_t I = std::distance(knots.data(),
                     std::upper_bound(knots.data(), knots.data() + knots.size(),  u)) - 1;
    const size_t L = numControlPts() - 1;
    return std::min(I, L - 1); // Clamp to the last curve segment [U_(L - 1), u_L]; this is needed when evaluating exactly at the upper knot bound.
}

template<size_t Dimension>
auto Spline_T<Dimension>::getControlPtsInfluencingSegment(size_t I) const ->std::vector<Vec> {
    const size_t n = degree();
    // HW4 TODO 1.1: extract the de Boor control points influencing the segment [u_I, u_{I + 1}].
    std::vector<Vec> d0;
    for (int i = I - (n - 1); i <= I + 1; i++) {
        d0.push_back(controlPts.row(i));
    }

    return d0;
}

template<size_t Dimension>
typename Spline_T<Dimension>::Vec Spline_T<Dimension>::evalPt(EvalMethod method, double u) const {
    const size_t n = degree();
    int I;

    if (n == 0 && controlPts.rows() == 1) {
        return controlPts.row(0);
    }
    
    if (u < knots[0]) {
        I = -1;
    }
    else {
        I = findSegmentContaining(u);
    }

    if (n == 0) {
        return controlPts.row(I + 1);
    }


    if (method == EvalMethod::DE_BOOR) {
        // HW4 TODO 1.1: implement the triangular scheme of the de Boor algorithm.
        d.clear();

        std::vector<Vec> d0 = Spline_T<Dimension>::getControlPtsInfluencingSegment(I);
        d.push_back(d0);

        for (int j = 1; j <= n; j++) {
            std::vector<Vec> emptyVec;
            d.push_back(emptyVec);

            for (int i = 0; i <= n - j; i++) {
                double alpha = (u - knots[I - n + j + i]) / (knots[I + 1 + i] - knots[I - n + j + i]);
                Vec interpolationPt = (1 - alpha) * d[j - 1][i] + alpha * d[j - 1][i + 1];

                d[j].push_back(interpolationPt);
            }
        }

        return d[n][0];
    }
    if (method == EvalMethod::BASIS) {
        // HW4 TODO 1.3: evaluate the B-Spline using the summation formula and your
        // `bsplineBasisFunction` implementation.
        Vec su = Vec::Zero();

        for (size_t j = I - (n - 1); j <= I + 1; j++) {
            double basisVal = bsplineBasisFunction(n, j, u);
            Vec point = controlPts.row(j) * basisVal;

            su = su + point;
        }

        return su;
    }
    throw std::runtime_error("Unimplemented method");
}

template<size_t Dimension>
void SplineBase<Dimension>::eval(size_t resolution, EvalMethod method, std::vector<Vec> &result) const {
    result.clear();

    // HW4 TODO 1.1: evaluate the curve at `resolution` equispaced points between
    // `domainStart()` and `domainEnd()` by calling `evalPt`.
    double stride = (domainEnd() - domainStart()) / (resolution - 1);
    double u = domainStart();

    for (int r = 0; r < resolution; r++) {
        Vec dn = evalPt(method, u);
        result.push_back(dn);

        u += stride;
    }
}


template<size_t Dimension>
void Spline_T<Dimension>::visualizeDeBoor(double u, MXDf &V, Eigen::MatrixX2i &E) const {
    // HW4 TODO 1.2: Visualization of the generations of the de Boor algorithm.
    evalPt(EvalMethod::DE_BOOR, u);

    std::vector<Vec> points;
    std::vector<Eigen::Vector2i> edges;

    for (int j = 1; j < d.size() - 1; j++) {
        for (int i = 0; i < d[j].size() - 1; i++) {
            points.push_back(d[j][i]);
            points.push_back(d[j][i + 1]);

            edges.emplace_back(points.size() - 2, points.size() - 1);
        }
    }

    V.resize(points.size(), 3);
    E.resize(edges.size(), 2);

    for (size_t i = 0; i < points.size(); ++i) V.row(i) = points[i];
    for (size_t i = 0; i < edges.size(); ++i) E.row(i) = edges[i];
}

template<size_t Dimension>
void Spline_T<Dimension>::inferKnots() {
    const int n = degree();
    const int L = numControlPts() - 1;
    const int K = L + n - 1;
    knots.setLinSpaced(K + 1, 0, 1);


    // HW4 TODO 1.4.1 - 1.4.3: fill in `knots` by inferring the parameter values
    // using various approaches depending on the `paramType` and `repeatedEndKnots`
    // member variables.

    if (paramType == ParametrizationType::UNIFORM) {
        if (repeatedEndKnots == true) {
            knots = VXf::Zero(K + 1);

            double uPos = 0.0;
            double spacing = 1.0 / (K - 2 * (n - 1));

            for (int i = 0; i < K + 1; i++) {
                if (i < n) {
                    knots(i) = 0.0;
                }
                else if (i > K - n) {
                    knots(i) = 1.0;
                }
                else {
                    uPos += spacing;
                    knots(i) = uPos;
                }
            }
        }
    }

    else if (paramType == ParametrizationType::CHORDLEN || paramType == ParametrizationType::CENTRIPETAL) {
        knots = VXf::Zero(K + 1);
        std::vector<double> knotsInterval;
        double intervalTotal = 0.0;

        if (repeatedEndKnots == true) {
            for (int i = 0; i < K; i++) {
                if (i >= 0 && i < n - 1) {
                    knotsInterval.push_back(0.0);
                }
                else if (i >= L && i < K) {
                    knotsInterval.push_back(0.0);
                }
                else {
                    double edgeSum = 0.0;
                    for (int c = i - (n - 1); c < i + 1; c++) {
                        double edgeLength = (controlPts.row(c + 1) - controlPts.row(c)).norm();
                        edgeSum += edgeLength;
                    }

                    if (paramType == ParametrizationType::CENTRIPETAL) {
                        edgeSum = std::sqrt(edgeSum);
                    }

                    knotsInterval.push_back(edgeSum);
                    intervalTotal += edgeSum;
                }
            }
        }
        else if (repeatedEndKnots == false) {
            for (int i = 0; i < K; i++) {
                double edgeSum = 0.0;

                for (int c = i - (n - 1); c < i + 1; c++) {
                    if (c < 0 || c >= L) {
                        continue;
                    }

                    double edgeLength = (controlPts.row(c + 1) - controlPts.row(c)).norm();
                    edgeSum += edgeLength;
                }

                if (paramType == ParametrizationType::CENTRIPETAL) {
                    edgeSum = std::sqrt(edgeSum);
                }

                knotsInterval.push_back(edgeSum);
                intervalTotal += edgeSum;
            }
        }

        for (int i = 1; i < knots.size(); i++) {
            knotsInterval[i - 1] = knotsInterval[i - 1] / intervalTotal;
            knots(i) = knots(i - 1) + knotsInterval[i - 1];
        }
    }


    Base::validateKnots(knots);
}


template<size_t Dimension>
void Spline_T<Dimension>::insertKnot(double ubar) {
    // HW4 TODO 2.1: Calculate the new, expanded `knots` and `controlPts`
    // produced by inserting the knot `ubar`.

    size_t I = findSegmentContaining(ubar);
    size_t n = degree();

    // Replace Control Points
    MXDf controlPts_new = MXDf::Zero(controlPts.rows() + 1, Dimension);
    evalPt(EvalMethod::DE_BOOR, ubar);

    int replacePos = I - n + 2;
    for (int j = 0; j < d[1].size(); j++) {
        controlPts_new.row(replacePos) = d[1][j];
        replacePos += 1;
    }

    for (int j = 0; j < I - n + 2; j++) {
        controlPts_new.row(j) = controlPts.row(j);
    }

    int newInd = controlPts_new.rows() - 1;
    for (int j = controlPts.rows() - 1; j > I; j--) {
        controlPts_new.row(newInd) = controlPts.row(j);
        newInd -= 1;
    }

    // Add a knot
    VXf knots_new = VXf::Zero(knots.size() + 1);

    for (int j = 0; j < knots.size(); j++) {
        if (j <= I) {
            knots_new(j) = knots(j);
        }
        else {
            knots_new(j + 1) = knots(j);
        }

        if (I == j) {
            knots_new(j + 1) = ubar;
        }
    }


    // Apply the new curve data.
    knots = knots_new;
    controlPts = controlPts_new;

    Base::validateSizes();
}

template<size_t Dimension>
std::vector<BezierCurve> splineToBezierSegments(Spline_T<Dimension> spline /* copy is intentional */) {
    auto& knots = spline.knots;
    const size_t n = spline.degree();

    // HW4 TODO 2.2: extract the Bézier curve representation of each polynomial piece
    // of the B-Spline `spline`. This should be done by inserting each knot up
    // to multiplicity `m` and then extracting the de Boor points influencing
    // each nonempty curve segment.
    std::vector<BezierCurve> result;

    size_t L = spline.numControlPts() - 1;

    // Multiplicity to n
    for (int i = n - 1; i <= L; i++) {
        int m = 0;

        for (int j = 0; j < knots.size(); j++) {
            if (knots[i] == knots[j]) {
                m += 1;
            }
        }

        int insertTimes = n - m;
        for (int j = 0; j < insertTimes; j++) {
            spline.insertKnot(knots[i]);
        }

        L = spline.numControlPts() - 1;
    }


    // Convert to Bezier Curve
    L = spline.numControlPts() - 1;
    for (int i = n - 1; i < L; i++) {
        if (knots[i + 1] - knots[i] != 0) {
            std::vector<Eigen::Vector3f> segPts = spline.getControlPtsInfluencingSegment(i);
            BezierCurve bCurve = BezierCurve(segPts);

            result.push_back(bCurve);
        }
    }

    return result;
}

////////////////////////////////////////////////////////////////////////////////
// NURBS functionality (using homogeneous coordinates)
////////////////////////////////////////////////////////////////////////////////
template<size_t Dimension>
typename NURBS<Dimension>::Vec NURBS<Dimension>::controlPt(size_t i) const {
    // HW4 TODO 4: get the Euclidean coordinates corresponding to the homogenized
    // coordinate vector `m_homogeneousSpline.controlPt(i)`.
    Vec euclideanPt = Vec::Zero(Dimension);

    for (int j = 0; j < euclideanPt.size(); j++) {
        euclideanPt(j) = m_homogeneousSpline.controlPt(i)(j) / m_homogeneousSpline.controlPt(i)(Dimension);
    }

    return euclideanPt;
}

template<size_t Dimension>
typename NURBS<Dimension>::MXDf NURBS<Dimension>::getControlPts() const {
    // HW4 TODO 4: get the Euclidean coordinates corresponding to all of the
    // homogenized coordinate vectors stored in the rows of
    // `m_homogeneousSpline.getControlPts()`.
    MXDf euclideanPts = MXDf::Zero(m_homogeneousSpline.numControlPts(), Dimension);

    for (int i = 0; i < euclideanPts.rows(); i++) {
        euclideanPts.row(i) = controlPt(i);
    }

    return euclideanPts;
}

template<size_t Dimension>
void NURBS<Dimension>::setControlPt(size_t i, const Vec &v) {
    // HW4 TODO 4: set the Euclidean coordinates of control point `i`
    // while leaving its weight unchanged.

    HSpline::Vec homogeneousPt = HSpline::Vec::Zero(Dimension + 1);

    for (int j = 0; j < Dimension; j++) {
        homogeneousPt(j) = v(j) * weight(i);
    }
    homogeneousPt(Dimension) = weight(i);

    m_homogeneousSpline.setControlPt(i, homogeneousPt);
}

template<size_t Dimension>
void NURBS<Dimension>::setWeight(size_t i, float w) {
    // HW4 TODO 4: set the weight of control point `i` while leaving its Euclidean
    // coordinates unchanged.

    HSpline::Vec homogeneousPt = HSpline::Vec::Zero(Dimension + 1);
    Vec euclideanPt = controlPt(i);

    for (int j = 0; j < Dimension; j++) {
        homogeneousPt(j) = euclideanPt(j) * w;
    }
    homogeneousPt(Dimension) = w;

    m_homogeneousSpline.setControlPt(i, homogeneousPt);
}

template<size_t Dimension>
typename NURBS<Dimension>::Vec NURBS<Dimension>::evalPt(EvalMethod method, double u) const {
    // HW4 TODO 4: evaluate s(u) by first evaluating a point in homogenized
    // coordinates using `m_homogeneousSpline.evalPt` and then doing a
    // perspective divide.
    HSpline::Vec su = m_homogeneousSpline.evalPt(method, u);
    Vec euclideanSu = Vec::Zero(Dimension);

    for (int i = 0; i < Dimension; i++) {
        euclideanSu(i) = su(i) / su(Dimension);
    }

    return euclideanSu;
}

// Calculate the derivative of a spline.
template<size_t Dimension>
Spline_T<Dimension> derivative(const Spline_T<Dimension> &s) {
    using VXf  = typename Spline_T<Dimension>::VXf;
    using MXDf = typename Spline_T<Dimension>::MXDf;

    // HW5 TODO 2.1: calculate the one-lower degree spline that is the derivative of `s`.
    // Note: if you've copied this method into your old `spline.cc`, be sure also to copy
    // the two explicit template instantiations of this method at the bottom of the file.
    size_t n = s.degree();
    VXf  u = s.knots;
    MXDf d = s.getControlPts();


    // Special Case if attempting to differentiate degree 0
    if (n == 0) {
        for (int i = 0; i < d.rows(); i++) {
            d.row(i) = Eigen::Matrix<float, Dimension, 1>::Zero();
        }
        return Spline_T<Dimension>(n, u, d);
    }


    // Derivatives Control Points
    std::vector<Eigen::Matrix<float, Dimension, 1>> points;

    for (int i = 0; i < d.rows() - 1; i++) {
        Eigen::Matrix<float, Dimension, 1> dDelta = d.row(i + 1) - d.row(i);
        double uDelta = u[i + n] - u[i];

        Eigen::Matrix<float, Dimension, 1> newPt = n * (dDelta / uDelta);
        points.push_back(newPt);
    }

    MXDf newPts;
    newPts.resize(points.size(), Dimension);

    for (size_t i = 0; i < points.size(); ++i) newPts.row(i) = points[i];

    // Derivatives Knots
    VXf newKnots = VXf::Zero(u.size() - 2);
    for (int i = 0; i < newKnots.size(); i++) {
        newKnots[i] = u[i + 1];
    }



    return Spline_T<Dimension>(n - 1, newKnots, newPts);
}

////////////////////////////////////////////////////////////////////////////////
// Explicit instantiation of all class and function templates and functions.
////////////////////////////////////////////////////////////////////////////////
template struct SplineBase<3>;
template struct Spline_T<3>;
template struct SplineBase<4>;
template struct Spline_T<4>;
template struct NURBS<3>;
template std::vector<BezierCurve> splineToBezierSegments<3>(Spline_T<3> spline);

template Spline_T<3> derivative(const Spline_T<3> &s);
template Spline_T<4> derivative(const Spline_T<4> &s);
