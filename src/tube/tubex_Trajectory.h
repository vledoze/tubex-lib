/* ============================================================================
 *  tubex-lib - Trajectory class
 * ============================================================================
 *  Copyright : Copyright 2017 Simon Rohou
 *  License   : This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 *
 *  Author(s) : Simon Rohou
 *  Bug fixes : -
 *  Created   : 2018
 * ---------------------------------------------------------------------------- */

#ifndef TRAJECTORY_HEADER
#define TRAJECTORY_HEADER

#include <map>
#include "ibex_Interval.h"

namespace tubex
{
  class Trajectory
  {
    public:

      /**
       * \brief Default constructor
       */
      Trajectory(const std::string& name = "", const std::string& color = "blue");
      Trajectory(const std::map<double,double>& m_map_values, const std::string& name = "", const std::string& color = "blue");
      const std::string& color() const;
      const std::string& name() const;
      void setName(const std::string& name);
      const std::map<double,double> getMap() const;
      const ibex::Interval domain() const;
      double& set(double t, double y);
      const double operator[](double t) const;
      void truncateDomain(const ibex::Interval& domain);
      void shiftDomain(double shift_ref);


    protected:

      /** Class variables **/

        std::map<double,double> m_map_values;

        // Graphics attributes
        std::string m_name, m_color;
        static int nb_traj;
        static std::vector<std::string> v_traj_names;
  };
}

#endif