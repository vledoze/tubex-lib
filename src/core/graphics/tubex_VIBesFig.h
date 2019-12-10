/** 
 *  \file
 *  VIBesFig class
 * ----------------------------------------------------------------------------
 *  \date       2016
 *  \author     Simon Rohou
 *  \copyright  Copyright 2019 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#ifndef __TUBEX_VIBESFIG_H__
#define __TUBEX_VIBESFIG_H__

#include "tubex_Figure.h"
#include "tubex_Polygon.h"
#include "vibes.h"

namespace tubex
{
  /**
   * \class VIBesFig
   * \brief Two-dimensional graphical item based on the VIBes viewer
   */
  class VIBesFig : public Figure
  {
    public:

      /// \name Definition and properties
      /// @{

      /**
       * \brief Creates a VIBesFig
       *
       * \param fig_name a reference to the figure that will be displayed in the window's title
       */
      VIBesFig(const std::string& fig_name);

      /**
       * \brief VIBesFig destructor
       */
      ~VIBesFig();

      /**
       * \brief Closes this figure
       */
      void close();
      
      /**
       * \brief Sets the properties (coordinates and dimensions) of this figure
       *
       * \param x horizontal coordinate (in pixels)
       * \param y vertical coordinate (in pixels)
       * \param width width value (in pixels)
       * \param height height value (in pixels)
       */
      void set_properties(int x, int y, int width, int height);

      /**
       * \brief Sets the axis limits of this figure
       *
       * The function updates the viewbox and applies the changes.
       *
       * \param x_min lower horizontal value to be displayed
       * \param x_max upper horizontal value to be displayed
       * \param y_min lower vertical value to be displayed
       * \param y_max upper vertical value to be displayed
       * \param same_ratio if `true`, will compute the min/max values so
       *        that the previous ratio will be preserved (false by default)
       * \param margin adds a custom margin to the view box (none by default)
       * \return the updated view box of this figure
       */
      const ibex::IntervalVector& axis_limits(double x_min, double x_max, double y_min, double y_max, bool same_ratio = false, float margin = 0.);

      /**
       * \brief Sets the axis limits of this figure
       *
       * The function updates the viewbox and applies the changes.
       *
       * \param viewbox the 2d box defining lower/upper horizontal/vertical values
       * \param same_ratio if `true`, will compute the min/max values so
       *        that the previous ratio will be preserved (false by default)
       * \param margin adds a custom margin to the view box (none by default)
       * \return the updated view box of this figure
       */
      const ibex::IntervalVector& axis_limits(const ibex::IntervalVector& viewbox, bool same_ratio = false, float margin = 0.);

      /// @}
      /// \name Saving this figure
      /// @{

      /**
       * \brief Saves the figure in SVG/PNG/... format
       *
       * A file named {path}/{figure_name}{suffix}.{extension} will be created in the current directory.
       *
       * \param suffix optional part name that can be added to the figure name (none by default)
       * \param extension optional part to specify the type of the image ("svg" by default)
       * \param path optional path to a different directory ("." by default)
       */
      void save_image(const std::string& suffix = "", const std::string& extension = "svg", const std::string& path = ".") const;

      /// @}
      /// \name Figure's content
      /// @{

      /**
       * \brief Displays this figure
       */
      void show() {};

      /**
       * \brief Clears this figure by removing displayed items
       */
      void clear();
      
      /// @}
      /// \name Displaying objects
      /// @{

      /**
       * \brief Draws a box
       *
       * \param box the 2d IntervalVector to be displayed
       * \param params VIBes parameters related to the box
       */
      void draw_box(const ibex::IntervalVector& box, const vibes::Params& params);

      /**
       * \brief Draws a box
       *
       * \param box the 2d IntervalVector to be displayed
       * \param color the optional color of the box (black by default) 
       * \param params VIBes parameters related to the box (none by default)
       */
      void draw_box(const ibex::IntervalVector& box, const std::string& color = "", const vibes::Params& params = vibes::Params());

      /**
       * \brief Draws a line of points
       *
       * \param v_x vector of horizontal coordinates
       * \param v_y vector of vertical coordinates
       * \param params VIBes parameters related to the line
       */
      void draw_line(const std::vector<double>& v_x, const std::vector<double>& v_y, const vibes::Params& params);

      /**
       * \brief Draws a line of points
       *
       * \param v_x vector of horizontal coordinates
       * \param v_y vector of vertical coordinates
       * \param color the optional color of the line (black by default) 
       * \param params VIBes parameters related to the line (none by default)
       */
      void draw_line(const std::vector<double>& v_x, const std::vector<double>& v_y, const std::string& color = "", const vibes::Params& params = vibes::Params());

      /**
       * \brief Draws a circle
       *
       * \param x horizontal coordinate
       * \param y vertical coordinate
       * \param r radius
       * \param params VIBes parameters related to the circle
       */
      void draw_circle(double x, double y, double r, const vibes::Params& params);

      /**
       * \brief Draws a circle
       *
       * \param x horizontal coordinate
       * \param y vertical coordinate
       * \param r radius
       * \param color the optional color of the circle (black by default) 
       * \param params VIBes parameters related to the circle (none by default)
       */
      void draw_circle(double x, double y, double r, const std::string& color = "", const vibes::Params& params = vibes::Params());

      /**
       * \brief Draws a polygon
       *
       * \param p polygon
       * \param params VIBes parameters related to the polygon (none by default)
       */
      void draw_polygon(const Polygon& p, const vibes::Params& params);

      /**
       * \brief Draws a polygon
       *
       * \param p polygon
       * \param color the optional color of the polygon (black by default) 
       * \param params VIBes parameters related to the polygon (none by default)
       */
      void draw_polygon(const Polygon& p, const std::string& color = "", const vibes::Params& params = vibes::Params());

      /**
       * \brief Draws a point
       *
       * \param p the 2d Point to be displayed
       * \param params VIBes parameters related to the point
       */
      void draw_point(const Point& p, const vibes::Params& params);

      /**
       * \brief Draws a point
       *
       * \param p the 2d Point to be displayed
       * \param color the optional color of the point (black by default) 
       * \param params VIBes parameters related to the point (none by default)
       */
      void draw_point(const Point& p, const std::string& color = "", const vibes::Params& params = vibes::Params());

      /**
       * \brief Draws a point
       *
       * \param p the 2d Point to be displayed
       * \param size display size of the points
       * \param params VIBes parameters related to the point
       */
      void draw_point(const Point& p, float size, const vibes::Params& params);

      /**
       * \brief Draws a point
       *
       * \param p the 2d Point to be displayed
       * \param size display size of the points
       * \param color the optional color of the point (black by default) 
       * \param params VIBes parameters related to the point (none by default)
       */
      void draw_point(const Point& p, float size, const std::string& color = "", const vibes::Params& params = vibes::Params());

      /**
       * \brief Draws a set of points
       *
       * \param v_pts vector of Point objects
       * \param size display size of the points
       * \param params VIBes parameters related to the set (none by default)
       */
      void draw_points(const std::vector<Point>& v_pts, float size, const vibes::Params& params);

      /**
       * \brief Draws a set of points
       *
       * \param v_pts vector of Point objects
       * \param size display size of the points
       * \param color the optional color of the set (black by default) 
       * \param params VIBes parameters related to the set (none by default)
       */
      void draw_points(const std::vector<Point>& v_pts, float size, const std::string& color = "", const vibes::Params& params = vibes::Params());

      /// @}
  };
}

#endif