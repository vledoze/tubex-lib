/* ================================================================
 *  tubex-lib - pyibex wrapper
 * ================================================================
 *  Copyright : ENSTA Bretagne (France)
 *  License   : This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 *
 *  Author(s) : Simon Rohou, Benoît Desrochers
 *  Bug fixes : -
 *  Created   : 2019
 * ------------------------------------------------------------- */

#include "tubex.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>

namespace py = pybind11;
using namespace pybind11::literals;
using py::class_;
using py::init;

using namespace tubex;

using ibex::Interval;
using ibex::IntervalVector;

PYBIND11_MODULE(tube, m)
{
  py::class_<Tube>(m, "Tube");
}