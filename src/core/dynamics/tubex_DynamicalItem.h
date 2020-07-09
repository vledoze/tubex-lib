/** 
 *  \file
 *  DynamicalItem class
 * ----------------------------------------------------------------------------
 *  \date       2018
 *  \author     Simon Rohou
 *  \copyright  Copyright 2020 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#ifndef __TUBEX_DYNAMICALITEM_H__
#define __TUBEX_DYNAMICALITEM_H__

#include "ibex_Interval.h"
#include "ibex_IntervalVector.h"

namespace tubex
{
  /**
   * \class DynamicalItem
   * \brief Abstract class for common properties of Tube, TubeVector,
   *        Slice, Trajectory, TrajectoryVector objects
   */
  class DynamicalItem
  {
    public:

      /**
       * \brief DynamicalItem destructor
       */
      virtual ~DynamicalItem();

      /**
       * \brief Returns the dimension of the object
       *
       * \return n
       */
      virtual int size() const = 0;

      /**
       * \brief Returns the temporal definition domain of this object
       *
       * \return an Interval object \f$[t_0,t_f]\f$
       */
      virtual const ibex::Interval tdomain() const = 0;

      /**
       * \brief Returns the box of the feasible values along \f$[t_0,t_f]\f$
       *
       * \note Used for genericity purposes
       *
       * \return the envelope of codomain values
       */
      virtual const ibex::IntervalVector codomain_box() const = 0;

      /**
       * \brief Returns the name of this class
       *
       * \note Only used for some generic display method
       *
       * \return the predefined name
       */
      virtual const std::string class_name() const = 0;

      /**
       * \brief Verifies that this interval is a feasible tdomain
       *
       * \note The tdomain must be non-empty, bounded and not degenerated
       * \todo Allow unbounded tdomains such as \f$[t_0,\infty]\f$?
       *
       * \param tdomain temporal domain \f$[t_0,t_f]\f$ to be tested
       * \return true in case of valid temporal tdomain
       */
      static bool valid_tdomain(const ibex::Interval& tdomain);
  };
}

#endif