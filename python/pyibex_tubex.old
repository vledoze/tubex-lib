//============================================================================
//                               P Y I B E X
// File        : pyIbex_Tube.cpp
// Author      : Benoit Desrochers
// Copyright   : Benoit Desrochers
// License     : See the LICENSE file
// Created     : Feb 10, 2016
//============================================================================

#include "tubex.h"
//#include "VibesFigure_Tube.h"
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
/*
std::vector<Interval> invert_wapper(Tube& o, const ibex::Interval& itv) {
  std::vector<Interval> res;
  o.invert(itv, res);
  return res;
}

Tube subTube(const Tube& o, const Interval& t){
  Tube tube = Tube(t, o.dt());
  int idx0  = o.input2index(t.lb());
  // std::cerr << t << " " << o.input2index(t.lb()) << " " << o.input2index(t.ub()) << "\n";
  for (int i = idx0; i < o.input2index(t.ub()); i++){
    tube.set(o[i],i-idx0);
  }
  return tube;
}

std::vector<Interval> tolist(const Tube& o){
  std::vector<Interval> res(o.size());
  for (int i = 0; i < res.size(); i++){
    res[i] = o[i];
  }
  return res;
}

struct TubeIterator{
  TubeIterator(Tube& tube) : tube(tube), cpt(0) {}
  // TubeIterator(Tube& tube) : tube(tube), cpt(0) {};
  Interval next() {
    // std::cout << "nExt " << cpt<< " " << tube.size() <<  "\n";
    if (cpt < tube.size()) {
      // Interval itv = tube[cpt];
      // cpt++;
      return tube[cpt++];
    } else {
      throw py::stop_iteration();
    }
  };
  Tube& tube;
  int cpt;
};

PYBIND11_MODULE(tube, m)
{
  // py::module m("tube", "python binding of ibex-robotics Tube");
  // py::class_<Interval>(m, "Interval");

  py::class_<TubeIterator>(m, "TubeIterator")
    .def("__next__", &TubeIterator::next)
    ;


  py::class_<Tube>(m, "Tube")
    .def(init<const Interval&, double, const Interval&>(), py::arg("intv_t"), py::arg("time_step"),
            py::arg_t<ibex::Interval>("default_value", ibex::Interval::EMPTY_SET, "EMPTY_SET") )
            // see http://pybind11.readthedocs.io/en/latest/advanced.html#default-arguments-revisited
    .def(init<const Tube&>())
    .def(init<const Tube&, const Interval&>())
    .def(init<const std::string&> ())
    .def(init<const ibex::Interval&, double, const ibex::Function&, const ibex::Interval&>(),
          py::arg("domain"), py::arg("timestep"), py::arg("f"),
          py::arg_t<ibex::Interval>("default_value", ibex::Interval::EMPTY_SET, "EMPTY_SET")
          )
    .def("volume", &Tube::volume )
    .def("dist", &Tube::dist )
    .def("size", &Tube::size )
    .def("dt", &Tube::dt )
    .def("isSlice", &Tube::isSlice )
    .def("isEmpty", &Tube::isEmpty )
    .def("input2index", &Tube::input2index )
    .def("index2input", &Tube::index2input )
    .def("subTube", &subTube)
    .def("domain", (const ibex::Interval& (Tube::*) () const ) &Tube::domain )
    .def("domain", (const ibex::Interval& (Tube::*) (int) const ) &Tube::domain )
    .def("domain", (ibex::Interval (Tube::*) (double) const ) &Tube::domain )
    .def("image", &Tube::image)
    .def("tolist", &tolist)

    .def("__getitem__", [](const Tube & t, int idx) -> Interval { return t[idx];})
    .def("__getitem__", [](const Tube & t, double time) -> Interval { return t[time];})
    .def("__setitem__", [](Tube & t, int idx, Interval& itv) { t.set(itv, idx);})
    .def("__setitem__", [](Tube & t, double time, Interval& itv) { t.set(itv, time);})


    .def( "set", (void (Tube::*) (const ibex::Interval&, int ) )  &Tube::set)
    .def( "set", (void (Tube::*) (const ibex::Interval&, double) )  &Tube::set)
    .def( "set", (void (Tube::*) (const ibex::Interval&, const ibex::Interval&) )  &Tube::set,
                  py::arg("intv_y"), py::arg("intv_t")=ibex::Interval::ALL_REALS )

    .def("feed", (const ibex::Interval (Tube::*)(const ibex::Interval& , int ) ) &Tube::feed)
    .def("feed", (const ibex::Interval (Tube::*)(const ibex::Interval& , double ) ) &Tube::feed)
    .def("feed", (void (Tube::*)(const std::map<double,ibex::Interval>& ) ) &Tube::feed)
    .def("feed", (void (Tube::*)(const std::map<double,double>& , const ibex::Interval& ) ) &Tube::feed)
    .def("feed", (void (Tube::*)(const std::map<double,double>& , const std::map<double,double>& ) ) &Tube::feed)

    .def("primitive", &Tube::primitive, py::arg_t<ibex::Interval>("default_value", ibex::Interval(0), "[0,0]") )

    .def("integral", (ibex::Interval (Tube::*) (double ) const ) &Tube::integral, py::call_guard<py::gil_scoped_release>())
    .def("integral", (ibex::Interval (Tube::*) (const ibex::Interval& ) const) &Tube::integral, py::call_guard<py::gil_scoped_release>())
    .def("integral", (ibex::Interval (Tube::*) (const ibex::Interval&, const ibex::Interval& ) const) &Tube::integral, py::call_guard<py::gil_scoped_release>())

    //.def("timeIntegration", (ibex::Interval (Tube::*) (const ibex::Interval& , const ibex::Interval& ) const) &Tube::timeIntegration)
    .def("ctcFwdBwd", &Tube::ctcFwdBwd, py::arg("derivative_tube"), py::arg("initial_value")=ibex::Interval::ALL_REALS, py::call_guard<py::gil_scoped_release>())

    .def("eval", &Tube::eval)

    .def("ctcObs", (bool (Tube::*) (const Tube& , ibex::Interval&, const ibex::Interval&, bool) ) &Tube::ctcObs,
                  "derivative_tube"_a, "t"_a, "y"_a, "fwd_bwd"_a=true, py::call_guard<py::gil_scoped_release>())

    .def("invert", invert_wapper)
    .def("serialize", (bool (Tube::*)(const std::string&) const) &Tube::serialize)
    .def("serialize", (bool (Tube::*)(const std::string&, const std::map<double,double>&) const) &Tube::serialize)
    .def("serialize", (bool (Tube::*)(const std::string&, const std::vector<std::map<double,double> >&) const) &Tube::serialize)

    .def("__iter__", [](Tube& o) { return TubeIterator(o);})
    .def("__len__", [](Tube& o) { return o.size();})
    // __next__ methode is in TubeIterator object

    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self * py::self)
    .def(py::self / py::self)
    .def(float() * py::self)
    .def(py::self &= py::self)
    .def(py::self * double())
    // .def()
    .def("__add__", [](const Tube& x1,  const ibex::Interval& x2) { return  x1+x2;})
    .def("__radd__", [](const Tube& x1,  const ibex::Interval& x2) { return  x2+x1;})

    .def("__sub__", [](const Tube& x1,  const ibex::Interval& x2) { return  x1-x2;})
    .def("__rsub__", [](const Tube& x1,  const ibex::Interval& x2) { return  x2-x1;})

    .def("__mult__", [](const Tube& x1,  const ibex::Interval& x2) { return  x1*x2;})
    .def("__rmult__", [](const Tube& x1,  const ibex::Interval& x2) { return  x2*x1;})
    .def("__rsub__", [](const Tube& x1,  double x2) { return  x2 - x1;})
    .def("__eq__", [](const Tube& x1, Tube &x2) {return x1 == x2;})

    .def("__abs__", [](const Tube& x1){return abs(x1);})

    // std::pair<ibex::Interval,ibex::Interval> partialTimeIntegration(const ibex::Interval& t) const;
  ;

  py::class_<VibesFigure_Tube>(m, "VibesFigure_Tube")
    .def(init<const std::string&, Tube *>())
    .def("setColors", &VibesFigure_Tube::setColors,
          "slices_color"_a, "slices_contracted_color"_a="",
          "background_color"_a="lightGray[lightGray]",
          "truth_color"_a = "red")
    .def("show",  ( void (VibesFigure_Tube::*) ()  ) &VibesFigure_Tube::show)
    .def("show", ( void (VibesFigure_Tube::*) (bool, int, bool)  )&VibesFigure_Tube::show,
          "detail_slices"_a, "slices_limit"_a, "update_background"_a=true)
    .def_static("showTube", [](Tube *tube, const std::string& name = "", int x = 0, int y = 0){
                          VibesFigure_Tube::show(tube,name,x,y);}, "tube"_a, "name"_a = "", "x"_a=0, "y"_a=0)
    .def("showScalarValues", &VibesFigure_Tube::showScalarValues, "map_scalar_values"_a, "color"_a="red", "points_size"_a=0.)

    // }( void (VibesFigure_Tube::*) (Tube*, const std::string&, int, int)  )&VibesFigure_Tube::show)
          // "tube"_a, "name"_a = "", "x"_a=0, "y"_a=0):
  ;


  m.def("abs", ( Tube (*) (const Tube&) ) &abs);;
  m.def("abs", [] (const Tube& x){ return abs(x);});
  m.def("sqr", [] (const Tube& x){ return sqr(x);});
  m.def("sqrt", [] (const Tube& x){ return sqrt(x);});
  m.def("pow", [] (const Tube& x, int p){ return pow(x,p);});
  m.def("pow", [] (const Tube& x, double p){ return pow(x,p);});
  m.def("pow", [] (const Tube &x, const ibex::Interval& p){ return pow(x, p);});
  m.def("root", [] (const Tube& x, int p){ return root(x, p);});
  m.def("exp", [] (const Tube& x){ return exp(x);});
  m.def("log", [] (const Tube& x){ return log(x);});
  m.def("cos", [] (const Tube& x){ return cos(x);});
  m.def("sin", [] (const Tube& x){ return sin(x);});
  m.def("tan", [] (const Tube& x){ return tan(x);});
  m.def("acos", [] (const Tube& x){ return acos(x);});
  m.def("asin", [] (const Tube& x){ return asin(x);});
  m.def("atan", [] (const Tube& x){ return atan(x);});
  m.def("cosh", [] (const Tube& x){ return cosh(x);});
  m.def("sinh", [] (const Tube& x){ return sinh(x);});
  m.def("tanh", [] (const Tube& x){ return tanh(x);});
  m.def("acosh", [] (const Tube& x){ return acosh(x);});
  m.def("asinh", [] (const Tube& x){ return asinh(x);});
  m.def("atanh", [] (const Tube& x){ return atanh(x);});
  m.def("atan2", [] (const Tube& y, const Tube& x){ return atan2(y, x);});


// const ibex::Interval feed(const ibex::Interval& intv_y, int index);
// const ibex::Interval feed(const ibex::Interval& intv_y, double t);
// const std::pair<ibex::Interval,ibex::Interval> getEnclosedBounds(const ibex::Interval& intv_t = ibex::Interval::ALL_REALS) const;
// ibex::Interval setInversion(const ibex::Interval& intv_y) const;
// void setInversion(const ibex::Interval& intv_y, std::vector<ibex::Interval> &v_intv_t) const;
// Tube& operator|=(const Tube& x);
// Tube& operator&=(const Tube& x);
// void print(int precision = 0) const;
// friend std::ostream& operator<<(std::ostream& str, const Tube& x);
// Tube primitive(const ibex::Interval& initial_value = ibex::Interval(0.)) const;
// ibex::Interval timeIntegration(double t) const;
// ibex::Interval timeIntegration(const ibex::Interval& t) const;
// std::pair<ibex::Interval,ibex::Interval> partialTimeIntegration(const ibex::Interval& t) const;
// ibex::Interval timeIntegration(const ibex::Interval& t1, const ibex::Interval& t2) const;
// std::pair<ibex::Interval,ibex::Interval> partialTimeIntegration(const ibex::Interval& t1, const ibex::Interval& t2) const;

// bool ctcFwd(const Tube& derivative_tube);
// bool ctcBwd(const Tube& derivative_tube);
// bool ctcFwdBwd(const Tube& derivative_tube);
// bool ctcIn(const Tube& derivative_tube, ibex::Interval& y, ibex::Interval& t);
// bool ctcOut(const ibex::Interval& y, const ibex::Interval& t);

    // return m.ptr();

}
*/