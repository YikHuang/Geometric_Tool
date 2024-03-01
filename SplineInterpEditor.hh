////////////////////////////////////////////////////////////////////////////////
// SplineinterpEditor.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Spline interpolation user interface for ECS 278.
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
//  Company:  University of California, Davis
//  Created:  10/28/2022 10:34:11
*///////////////////////////////////////////////////////////////////////////////
#ifndef SPLINEINTERPEDITOR_HH
#define SPLINEINTERPEDITOR_HH

#include "EditorBase.hh"
#include "BezierEditor.hh"
#include "PolyinterpEditor.hh"
#include "SplineEditor.hh"
#include "spline.hh"
#include "spline_interpolation.hh"

struct SplineInterpEditor : public EditorBase {
    SplineInterpEditor(AssignmentGUI *gui, DataSource &ds) : EditorBase(gui), dataSource(ds) {
        interpolant = std::make_shared<Spline>();
        interpPts = ds.generate();
        recomputeInterpolant();
    }
    DataSource &dataSource;

    Spline::ParametrizationType paramType = Spline::ParametrizationType::UNIFORM;

    std::shared_ptr<Spline> interpolant;
    std::vector<Eigen::Vector3f> interpPts;

    struct Selection {
        int controlPt = 0;
        bool dragging = false;
    };

    int   resolution = 100;
    Selection selection;

    virtual void enter(IGLViewer &v) override {
        SplineEditor::knotViewer(v).attachEditor(this);
        SplineEditor::knotViewer(v).enable(interpolant);
        EditorBase::enter(v);
    }

    virtual void leave(IGLViewer &v) override {
        SplineEditor::knotViewer(v).disable();
        EditorBase::leave(v);
    }

    virtual void updateView(IGLViewer &v) override {
        auto &plotter        = *dynamic_cast<LinePlotter         *>(v.plugins[0]);
        auto &vectorRenderer = *dynamic_cast<VectorFieldRenderer *>(v.plugins[1]);
        vectorRenderer.clear();
        plotter.clear();
        auto &vdata = v.data();
        vdata.clear();
        float dpiScale = static_cast<Viewer &>(v).getDPIScale();

        Eigen::Vector4f cageColor = Eigen::Vector4f(0.0f, 0.0f, 0.0f, 1.0f);
        Eigen::Vector4f curveColor = Eigen::Vector4f(0.3, 0.6f, 1.0f, 1.0f);

        vdata.add_points(pointsToMatrixRows(interpPts).cast<double>(), /* C = */ BezierEditor::cageColor(true).head<3>().transpose().cast<double>());
        SplineEditor::drawSpline(v, resolution, Spline::EvalMethod::DE_BOOR, *interpolant, /* selected */ false, /* showCage */ true);

        SplineEditor::knotViewer(v).enable(interpolant);

        vdata.point_size = 10.0f * dpiScale;
    }

    virtual bool callback_key_pressed(IGLViewer &viewer, unsigned int key, int modifiers) override {
        key = std::toupper(key);

        bool handled = false;
        switch (key) {
            case '=':
            case '+':
            case '-':
            case '_':
                resolution = std::min(std::max(resolution + ((key == '=' || key == '+') ? 1 : -1), 2), 1000);
                handled = true;
        }
        if (handled) updateView(viewer);
        return handled;
    }

    virtual bool callback_mouse_down(IGLViewer &viewer, int button, int modifier) override {
        if (button == int(IGLViewer::MouseButton::Right))
            return false;

        MouseRay ray(viewer, viewer.down_mouse_x, viewer.down_mouse_y);

        float selectionRadius = 0.05f;
        float minDist = selectionRadius;

        float dist;
        size_t closestPt;
        ray.closestCollectionPtAndDistance(interpPts, dist, closestPt);
        if (dist < minDist) {
            minDist = dist;
            selection.controlPt = closestPt;
            selection.dragging  = true;
        }

        return selection.dragging;
    }

    void recomputeInterpolant() {
        *interpolant = naturalCubicSplineInterpolant(interpPts, paramType);
    }

    virtual bool callback_mouse_move(IGLViewer &viewer, int mouse_x, int mouse_y) override {
        if (selection.dragging) {
            Eigen::Vector3f &cp = interpPts[selection.controlPt];

            const auto &c = viewer.core();
            // Determine screen-space depth of control point (to be preserved by drag)
            Eigen::Vector3f win = igl::project(cp, c.view, c.proj, c.viewport);

            win.head<2>() << mouse_x, c.viewport(3) - mouse_y; // GLFW's y coordinate is flipped vs unproject_ray's
            cp = igl::unproject(win, c.view, c.proj, c.viewport);

            if (editMode == EditingMode::TWO_D)
                cp[2] = 0.0f; // Enforce a precisely zero component in 2D mode;
                              // if roundoff error introduces a small nonzero z component the
                              // intersection test will break...

            recomputeInterpolant();

            updateView(viewer);
            return true;
        }
        return false;
    }

    virtual bool callback_mouse_up(IGLViewer &viewer, int button, int modifier) override {
        if (selection.dragging) {
            selection.dragging = false;
            return true;
        }

        return false;
    }

    virtual bool draw_menu(IGLViewer &viewer) override {
        bool handled = false, changed = false;
        if (ImGui::InputInt("Curve Eval Resolution", &resolution, 1)) {
            resolution = std::min(std::max(resolution, 2), 1000);
            changed = true;
        }
        int editModeInt = int(editMode);
        if (ImGui::Combo("Dimension", &editModeInt, "2D\0003D\0")) {
            editMode = EditingMode(editModeInt);
            applyEditingMode(viewer);
            if (editMode == EditingMode::TWO_D) {
                // When switching to 2D editing mode, we need to project the design onto the Y = 0 axis.
                for (auto &p : interpPts)
                    p[2] = 0.0;
                recomputeInterpolant();
                
                changed = true;
            }
            handled = true;
        }

        if (dataSource.draw_menu(viewer)) {
            interpPts = dataSource.generate();
            recomputeInterpolant();
            changed = true;
        }

        int ptypeInt = int(paramType);
        if (ImGui::Combo("Parametrization Type", &ptypeInt, "Uniform\0Chord Length\0Centripetal\0")) {
            paramType = Spline::ParametrizationType(ptypeInt);
            recomputeInterpolant();
            changed = true;
        }

        if (ImGui::Button("Convert to Bézier Curves")) {
            auto result = splineToBezierSegments(*interpolant);
            BezierEditor *bezierEditor = dynamic_cast<BezierEditor *>(gui->getEditor("bezier"));
            if (bezierEditor == nullptr) throw std::logic_error("Failed to access Bézier editor");
            bezierEditor->bezierCurves = result;
            gui->switchEditor("bezier");
        }

        if (changed) updateView(viewer);
        return changed || handled;
    }
};

#endif /* end of include guard: SPLINEINTERPEDITOR_HH */
