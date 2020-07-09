/** 
 *  \file
 *  VIBesFigPaving class
 * ----------------------------------------------------------------------------
 *  \date       2016
 *  \author     Simon Rohou
 *  \copyright  Copyright 2020 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#ifndef __TUBEX_VIBESFIGPAVING_H__
#define __TUBEX_VIBESFIGPAVING_H__

#include "tubex_VIBesFig.h"
#include "tubex_Paving.h"

namespace tubex
{
  /**
   * \class VIBesFigPaving
   * \brief Two-dimensional graphical item to display a Paving object
   *
   * One figure is linked to one paving, so that any update on this object
   * can be easily displayed on the figure.
   */
  class VIBesFigPaving : public VIBesFig
  {
    public:

      /**
       * \brief Creates a VIBesFigPaving
       *
       * \param fig_name a reference to the figure that will be displayed in the window's title
       * \param paving a const pointer to the paving to be displayed
       */
      VIBesFigPaving(const std::string& fig_name, const Paving *paving);

      /**
       * \brief Sets a custom color map
       *
       * \param color_map the color map `<paving_value,color_code>` to be set
       */
      void set_color_map(const std::map<SetValue,std::string>& color_map);

      /**
       * \brief Updates the display of the Paving object
       */
      void show();

    protected:
      
      /**
       * \brief Draws a paving object
       *
       * This method is recursive, because a Paving is implemented as a binary tree.
       *
       * \param paving a const pointer to the Paving object to be displayed
       */
      void draw_paving(const Paving *paving);

    protected:

      const Paving *m_paving; //!< const pointer to the object to be displayed
      std::map<SetValue,std::string> m_color_map; //!< custom color map
  };
}

#endif