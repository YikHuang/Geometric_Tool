////////////////////////////////////////////////////////////////////////////////
// Viewer.hh
////////////////////////////////////////////////////////////////////////////////
/*! @file
//  Simple wrapper around libigl's viewer that sets our preferred defaults.
*/
//  Author:  Julian Panetta (jpanetta), jpanetta@ucdavis.edu
////////////////////////////////////////////////////////////////////////////////
#ifndef VIEWER_HH
#define VIEWER_HH

#include <imgui.h>
#include "igl/unproject_ray.h"
#include <igl/opengl/glfw/Viewer.h>
#include <GLFW/glfw3.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

using IGLViewer       = igl::opengl::glfw::Viewer;
using IGLViewerPlugin = igl::opengl::glfw::ViewerPlugin;

struct MouseRay {
    // Construct the ray corresponding to a click location (mouse_x, mouse_y) in
    // screen space.
    MouseRay(IGLViewer &viewer, int mouse_x, int mouse_y) {
        const auto &c = viewer.core();
        mouse_y = c.viewport(3) - mouse_y; // GLFW's y coordinate is flipped vs unproject_ray's
        igl::unproject_ray(Eigen::Vector3f(mouse_x, mouse_y, -1.0f), c.view, c.proj, c.viewport, p, v);
        v.normalize();
    }

    // p x----->
    //      v  (unit direction)
    Eigen::Vector3f p, v;

    float closestRayParameter(const Eigen::Vector3f &pt) const {
        Eigen::Vector3f d = pt - p;
        float t = d.dot(v);
        // For query points in the opposite direction from the ray, closest
        // point is the start point.
        return std::max(t, 0.0f);
    }

    Eigen::Vector3f closestRayPt(const Eigen::Vector3f &pt) const {
        return p + v * closestRayParameter(pt);
    }

    float distanceToPt(const Eigen::Vector3f &pt) const {
        float t = closestRayParameter(pt);
        return std::sqrt((pt - p).squaredNorm() - t * t);
    }

    template<class Pt>
    void closestCollectionPtAndDistance(const std::vector<Pt> &pts, float &minDist, size_t &closest) const {
        minDist = std::numeric_limits<float>::max();
        closest = std::numeric_limits<size_t>::max();
        for (size_t i = 0; i < pts.size(); ++i) {
            float dist = distanceToPt(pts[i]);
            if (dist < minDist) {
                minDist = dist;
                closest = i;
            }
        }
    }

    template<class Derived>
    void closestCollectionPtAndDistance(const Eigen::MatrixBase<Derived> &pts, float &minDist, size_t &closest) const {
        minDist = std::numeric_limits<float>::max();
        closest = std::numeric_limits<size_t>::max();
        for (int i = 0; i < pts.rows(); ++i) {
            float dist = distanceToPt(pts.row(i).transpose().template cast<float>());
            if (dist < minDist) {
                minDist = dist;
                closest = i;
            }
        }
    }
};

struct Viewer : public IGLViewer {
    Viewer(const std::string &title = "ECS278 Viewer", const std::vector<IGLViewerPlugin*> &customPlugins = std::vector<IGLViewerPlugin*>())
        : windowTitle(title)
    {
        data().set_face_based(true);
        core().set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL);

        plugins = customPlugins;
        plugins.push_back(&imguiMenuPlugin);
        imguiMenuPlugin.widgets.push_back(&menu);

        launch_init(/* resizeable = */ true,
                    /* fullscreen = */ false,
                    /*       name = */ windowTitle);

        // Redraw when the window is resizing
        // This must be set after `launch_init` otherwise we'll try to redraw
        // before the menu plugin is initialized and crash.
        bool reentered_from_draw = false;
        IGLViewer::callback_post_resize = [this, reentered_from_draw](IGLViewer &v, int w, int h) mutable {
            if (Viewer::callback_post_resize)
                return Viewer::callback_post_resize(v, w, h);
            if (!reentered_from_draw) {
                reentered_from_draw = true;
                v.draw();
                glfwSwapBuffers(v.window);
                reentered_from_draw = false;
            }
            return true;
        };

        float menu_width = 275.f * menu.menu_scaling();
        // Customize the ImGui menu window
        menu.callback_draw_viewer_window = [&, menu_width]() {
          ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiCond_FirstUseEver);
          ImGui::SetNextWindowSize(ImVec2(menu_width, 0.0f), ImGuiCond_FirstUseEver);
          ImGui::SetNextWindowSizeConstraints(ImVec2(menu_width, 100.f), ImVec2(2 * menu_width, core().viewport[3]));
          bool _viewer_menu_visible = true;
          ImGui::Begin(
              "Settings",
              &_viewer_menu_visible,
              ImGuiWindowFlags_NoSavedSettings
          );
          ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.5f);
          if (menu.callback_draw_viewer_menu) { menu.callback_draw_viewer_menu(); }
          else { menu.draw_viewer_menu(); }
          ImGui::PopItemWidth();
          ImGui::End();
        };
    }

    // The libigl viewer specifies its current camera as a transformation of a
    // "base camera," which is chosen to approximately fit the mesh
    // data in view.
    // Normally they configure the base camera at `launch_init`, but this won't happen
    // for us since we call this before any mesh data is added. We expose this
    // convenience method for the user to reinitialize the base camera after
    // setting data.
    void recalculate_base_camera() { core().align_camera_center(data().V, data().F); }

    // We also provide a convenience method that allows automatically updating
    // the base camera when setting new mesh data.
    void set_mesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, bool update_base_camera = false) {
        data().set_mesh(V, F);
        if (update_base_camera) recalculate_base_camera();
    }

    bool run() { return launch_rendering(/* loop = */ true); }

    // Still allow the user to specify their own custom post-resize callback
    // (in addition to our default one that enables redraw during resize).
    std::function<bool(IGLViewer &viewer, int w, int h)> callback_post_resize;

    std::string windowTitle;

    igl::opengl::glfw::imgui::ImGuiPlugin imguiMenuPlugin;
    igl::opengl::glfw::imgui::ImGuiMenu menu;

    ~Viewer() {
        launch_shut();
    }

    float getDPIScale() const {
        float xscale, yscale;
        glfwGetWindowContentScale(window, &xscale, &yscale);
        return 0.5 * (xscale + yscale);
    }
};

#endif /* end of include guard: VIEWER_HH */
