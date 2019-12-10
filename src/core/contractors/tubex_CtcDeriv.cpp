/** 
 *  CtcDeriv class
 * ----------------------------------------------------------------------------
 *  \date       2015
 *  \author     Simon Rohou
 *  \copyright  Copyright 2019 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include "tubex_CtcDeriv.h"
#include "tubex_ConvexPolygon.h"

using namespace std;
using namespace ibex;

namespace tubex
{
  CtcDeriv::CtcDeriv() : Ctc()
  {

  }

  void CtcDeriv::contract(Tube& x, const Tube& v, TPropagation t_propa)
  {
    assert(x.domain() == v.domain());
    assert(Tube::same_slicing(x, v));

    Slice *x_slice = x.first_slice();
    const Slice *v_slice = v.first_slice();

    if(t_propa & FORWARD)
    {

      while(x_slice != NULL)
      {
        assert(v_slice != NULL);

        if(x_slice->domain().intersects(m_restricted_domain))
        {
          contract(*x_slice, *v_slice, t_propa);
          // If current slide is in restricted domain
          // Next slide must be in the restricted domain
          // Or we stop
          Slice *next_x_slice = x_slice->next_slice();
          if (next_x_slice == NULL)
            break;
          if (!next_x_slice->domain().intersects(m_restricted_domain))
            break;
        }
        x_slice = x_slice->next_slice();
        v_slice = v_slice->next_slice();
      }
    }
    
    if(t_propa & BACKWARD)
    {

      while(x_slice != NULL)
      {
        assert(v_slice != NULL);

        if(x_slice->domain().intersects(m_restricted_domain))
        {
          contract(*x_slice, *v_slice, t_propa);
          // If current slide is in restricted domain
          // Next slide must be in the restricted domain
          // Or we stop
          Slice *prev_x_slice = x_slice->prev_slice();
          if (prev_x_slice == NULL)
            break;
          if (!prev_x_slice->domain().intersects(m_restricted_domain))
            break;
        }
        x_slice = x_slice->prev_slice();
        v_slice = v_slice->prev_slice();
      }
    }
  }

  void CtcDeriv::contract(Slice& x, const Slice& v, TPropagation t_propa)
  {
    assert(x.domain() == v.domain());

    if(!x.domain().intersects(m_restricted_domain))
      return;

    Interval envelope = x.codomain(), ingate = x.input_gate(), outgate = x.output_gate();

    if(m_fast_mode) // Faster contraction without polygons
    {
      if(t_propa & FORWARD)
      {
        x.set_envelope(envelope & (ingate + Interval(0.,x.domain().diam()) * v.codomain()));
        x.set_output_gate(outgate & (ingate + x.domain().diam() * v.codomain()));
      }

      if(t_propa & BACKWARD)
      {
        x.set_envelope(envelope & (outgate - Interval(0.,x.domain().diam()) * v.codomain()));
        x.set_input_gate(ingate & (outgate - x.domain().diam() * v.codomain()));
      }
    }

    else // Optimal contraction
    {
      if(outgate == Interval::ALL_REALS) // Direct evaluation, polygons not needed
      {
        envelope &= ingate + Interval(0.,x.domain().diam()) * v.codomain();
        outgate &= ingate + x.domain().diam() * v.codomain();
      }

      else if(ingate == Interval::ALL_REALS) // Direct evaluation, polygons not needed
      {
        envelope &= outgate - Interval(0.,x.domain().diam()) * v.codomain();
        ingate &= outgate - x.domain().diam() * v.codomain();
      }

      else // Using polygons to compute the envelope
      {
        // todo: remove this: (or use Polygons with truncature)
        envelope &= Interval(-BOUNDED_INFINITY,BOUNDED_INFINITY);

        x.set_envelope(envelope);

        // Gates contraction
        outgate &= ingate + x.domain().diam() * v.codomain();
        ingate &= outgate - x.domain().diam() * v.codomain();

        // Gates needed for polygon computation
        x.set_input_gate(ingate);
        x.set_output_gate(outgate);

        // Optimal envelope
        envelope &= x.polygon(v).box()[1];

        // todo: remove this: (or use Polygons with truncature)
        if(envelope.ub() == BOUNDED_INFINITY) envelope = Interval(envelope.lb(),POS_INFINITY);
        if(envelope.lb() == -BOUNDED_INFINITY) envelope |= Interval(NEG_INFINITY,envelope.ub());
        if(ingate.ub() == BOUNDED_INFINITY) ingate = Interval(ingate.lb(),POS_INFINITY);
        if(ingate.lb() == -BOUNDED_INFINITY) ingate = Interval(NEG_INFINITY,ingate.ub());
        if(outgate.ub() == BOUNDED_INFINITY) outgate = Interval(outgate.lb(),POS_INFINITY);
        if(outgate.lb() == -BOUNDED_INFINITY) outgate = Interval(NEG_INFINITY,outgate.ub());
      }

      x.set_envelope(envelope);
      x.set_input_gate(ingate);
      x.set_output_gate(outgate);
    }
  }

  void CtcDeriv::contract_gates(Slice& x, const Slice& v)
  {
    assert(x.domain() == v.domain());
    
    Interval in_gate = x.input_gate(), out_gate = x.output_gate();

    Interval in_gate_proj = in_gate + x.domain().diam() * v.codomain();
    out_gate &= in_gate_proj;
    x.set_output_gate(out_gate);

    Interval out_gate_proj = out_gate - x.domain().diam() * v.codomain();
    in_gate &= out_gate_proj;
    x.set_input_gate(in_gate);
  }

  void CtcDeriv::contract(TubeVector& x, const TubeVector& v, TPropagation t_propa)
  {
    assert(x.size() == v.size());
    assert(x.domain() == v.domain());
    assert(TubeVector::same_slicing(x, v));

    for(int i = 0 ; i < x.size() ; i++)
      contract(x[i], v[i], t_propa);
  }
}