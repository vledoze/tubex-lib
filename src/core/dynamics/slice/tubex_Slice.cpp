/** 
 *  Slice class
 * ----------------------------------------------------------------------------
 *  \date       2018
 *  \author     Simon Rohou
 *  \copyright  Copyright 2019 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include <iomanip>
#include "tubex_Slice.h"
#include "tubex_CtcDeriv.h"
#include "tubex_arithmetic.h"

using namespace std;
using namespace ibex;

namespace tubex
{
  // Public methods

    // Definition

    Slice::Slice(const Interval& domain, const Interval& codomain)
      : m_domain(domain), m_codomain(codomain)
    {
      assert(valid_domain(domain));
      m_input_gate = new Interval(codomain);
      m_output_gate = new Interval(codomain);
    }

    Slice::Slice(const Slice& x)
      : Slice(x.domain(), x.codomain()) // in order to instantiate gates
    {
      *this = x;
    }

    Slice::~Slice()
    {
      // Links to other slices are destroyed
      if(m_prev_slice != NULL) m_prev_slice->m_next_slice = NULL;
      if(m_next_slice != NULL) m_next_slice->m_prev_slice = NULL;

      // Gates are deleted if not shared with other slices
      if(m_prev_slice == NULL) delete m_input_gate;
      if(m_next_slice == NULL) delete m_output_gate;
    }

    int Slice::size() const
    {
      return 1;
    }

    const Slice& Slice::operator=(const Slice& x)
    {
      m_domain = x.m_domain;
      m_codomain = x.m_codomain;
      *m_input_gate = *x.m_input_gate;
      *m_output_gate = *x.m_output_gate;
      
      if(m_synthesis_reference != NULL)
      {
        m_synthesis_reference->request_values_update();
        m_synthesis_reference->request_integrals_update();
      }
      
      return *this;
    }
    
    const Interval Slice::domain() const
    {
      return m_domain;
    }

    // Slices structure

    Slice* Slice::prev_slice()
    {
      return const_cast<Slice*>(static_cast<const Slice&>(*this).prev_slice());
    }

    const Slice* Slice::prev_slice() const
    {
      return m_prev_slice;
    }

    Slice* Slice::next_slice()
    {
      return const_cast<Slice*>(static_cast<const Slice&>(*this).next_slice());
    }

    const Slice* Slice::next_slice() const
    {
      return m_next_slice;
    }

    const Interval Slice::input_gate() const
    {
      return *m_input_gate;
    }
    
    const Interval Slice::output_gate() const
    {
      return *m_output_gate;
    }
    
    // Accessing values

    const Interval Slice::codomain() const
    {
      return m_codomain;
    }
    
    const IntervalVector Slice::box() const
    {
      IntervalVector box(2);
      box[0] = domain();
      box[1] = m_codomain;
      return box;
    }
    
    double Slice::diam(const Interval& interval)
    {
      if(interval.is_empty())
        return 0.;
      
      else if(interval.is_unbounded())
        return POS_INFINITY;
      
      else
        return interval.diam();
    }
    
    double Slice::diam() const
    {
      return diam(m_codomain);
    }
    
    double Slice::volume() const
    {
      return m_domain.diam() * diam(m_codomain);
    }

    const Interval Slice::operator()(double t) const
    {
      // Write access is not allowed for this operator
      // For write access: use set()
      assert(domain().contains(t));

      if(t == m_domain.lb())
        return input_gate();

      else if(t == m_domain.ub())
        return output_gate();

      return m_codomain;
    }

    const Interval Slice::operator()(const Interval& t) const
    {
      // Write access is not allowed for this operator
      // For write access: use set()
      assert(domain().is_superset(t));
      if(t.is_degenerated())
        return (*this)(t.lb());
      return m_codomain;
    }

    const pair<Interval,Interval> Slice::eval(const Interval& t) const
    {
      Interval intersection = t & m_domain;
      pair<Interval,Interval> p_eval = make_pair(Interval::EMPTY_SET, Interval::EMPTY_SET);

      if(!intersection.is_empty())
      {
        if(intersection.contains(m_domain.lb()))
        {
          p_eval.first |= input_gate().lb();
          p_eval.second |= input_gate().ub();
        }

        if(intersection.contains(m_domain.ub()))
        {
          p_eval.first |= output_gate().lb();
          p_eval.second |= output_gate().ub();
        }

        if(!intersection.is_degenerated()
          || (intersection != m_domain.lb() && intersection != m_domain.ub()))
        {
          p_eval.first |= m_codomain.lb();
          p_eval.second |= m_codomain.ub();
        }
      }
      
      return p_eval;
    }

    const Interval Slice::interpol(double t, const Slice& v) const
    {
      assert(domain() == v.domain());
      assert(domain().contains(t));
      return interpol(Interval(t), v);
    }

    const Interval Slice::interpol(const Interval& t, const Slice& v) const
    {
      assert(domain() == v.domain());
      assert(domain().is_superset(t));

      if(domain().is_subset(t))
        return codomain();

      else if(t.is_degenerated())
        return (output_gate() - (m_domain.ub() - t) * v.codomain())
             & (input_gate() + (t - m_domain.lb()) * v.codomain());

      else
      {
        IntervalVector input_box = box();
        input_box[0] &= t;
        return (polygon(v) & input_box)[1];
      }
    }

    const Interval Slice::invert(const Interval& y, const Interval& search_domain) const
    {
      Slice v(domain(), Interval::ALL_REALS); // todo: optimize this
      return invert(y, v, search_domain);
    }

    const Interval Slice::invert(const Interval& y, const Slice& v, const Interval& search_domain) const
    {
      assert(domain() == v.domain());
      // todo: use enclosed bounds also? in order to speed up computations

      if(!m_domain.intersects(search_domain))
        return Interval::EMPTY_SET;

      else if((m_domain & search_domain) == m_domain && m_codomain.is_subset(y))
        return m_domain;

      else if(search_domain == m_domain.lb())
      {
        if(y.intersects(input_gate()))
          return m_domain.lb();
        else
          return Interval::EMPTY_SET;
      }

      else if(search_domain == m_domain.ub())
      {
        if(y.intersects(output_gate()))
          return m_domain.ub();
        else
          return Interval::EMPTY_SET;
      }

      else
      {
        ConvexPolygon p = polygon(v);
        IntervalVector box(2);
        box[0] = search_domain; box[1] = y;
        box = p & box;
        return box[0];
      }
    }

    // Tests
    
    bool Slice::operator==(const Slice& x) const
    {
      return domain() == x.domain() &&
             codomain() == x.codomain() &&
             input_gate() == x.input_gate() &&
             output_gate() == x.output_gate();
    }
    
    bool Slice::operator!=(const Slice& x) const
    {
      return domain() != x.domain() ||
             codomain() != x.codomain() ||
             input_gate() != x.input_gate() ||
             output_gate() != x.output_gate();
    }

    #define macro_sets_comparison(f) \
      \
      assert(domain() == x.domain()); \
      return codomain().f(x.codomain()) \
          && input_gate().f(x.input_gate()) \
          && output_gate().f(x.output_gate()); \

    bool Slice::is_subset(const Slice& x) const
    {
      macro_sets_comparison(is_subset);
    }

    bool Slice::is_strict_subset(const Slice& x) const
    {
      return is_subset(x) && *this != x;
    }

    bool Slice::is_interior_subset(const Slice& x) const
    {
      macro_sets_comparison(is_interior_subset);
    }

    bool Slice::is_strict_interior_subset(const Slice& x) const
    {
      return is_interior_subset(x) && *this != x;
    }

    bool Slice::is_superset(const Slice& x) const
    {
      macro_sets_comparison(is_superset);
    }

    bool Slice::is_strict_superset(const Slice& x) const
    {
      return is_superset(x) && *this != x;
    }
    
    bool Slice::is_empty() const
    {
      return m_codomain.is_empty() || input_gate().is_empty() || output_gate().is_empty();
    }

    const BoolInterval Slice::contains(const Trajectory& x) const
    {
      assert(domain().is_subset(x.domain()));

      Interval traj_domain = x(m_domain);
      // We use x(Interval(double)) for reliable evaluation:
      Interval traj_input = x(Interval(m_domain.lb()));
      Interval traj_output = x(Interval(m_domain.ub()));

      if(!traj_domain.intersects(m_codomain)
      || !traj_input.intersects(input_gate())
      || !traj_output.intersects(output_gate()))
        return NO;

      else
      {
        if(!traj_input.is_subset(input_gate()) || !traj_output.is_subset(output_gate()))
          return MAYBE;

        else if(traj_domain.is_subset(m_codomain))
          return YES;

        else // too much pessimism for the trajectory evaluation on m_domain
        {
          // Bissections are performed to reach an accurate evaluation

          list<Interval> s_subdomains;
          s_subdomains.push_front(m_domain);

          while(!s_subdomains.empty())
          {
            Interval t = s_subdomains.front();
            s_subdomains.pop_front();

            Interval thinner_eval = x(t);

            if(!thinner_eval.intersects(m_codomain))
            {
              return NO;
            }

            else if(!thinner_eval.is_subset(m_codomain))
            {
              if(t.diam() < EPSILON_CONTAINS)
                return MAYBE;

              s_subdomains.push_front(Interval(t.lb(), t.lb() + t.diam() / 2.));
              s_subdomains.push_front(Interval(t.lb() + t.diam() / 2., t.ub()));
            }
          }

          return YES;
        }
      }
    }

    // Setting values

    void Slice::set(const Interval& y)
    {
      m_codomain = y;

      *m_input_gate = y;
      if(prev_slice() != NULL)
        *m_input_gate &= prev_slice()->codomain();

      *m_output_gate = y;
      if(next_slice() != NULL)
        *m_output_gate &= next_slice()->codomain();

      if(m_synthesis_reference != NULL)
      {
        m_synthesis_reference->request_values_update();
        m_synthesis_reference->request_integrals_update();
      }
    }
    
    void Slice::set_empty()
    {
      set(Interval::EMPTY_SET);
    }

    void Slice::set_envelope(const Interval& envelope)
    {
      m_codomain = envelope;
      *m_input_gate &= m_codomain;
      *m_output_gate &= m_codomain;

      if(m_synthesis_reference != NULL)
      {
        m_synthesis_reference->request_values_update();
        m_synthesis_reference->request_integrals_update();
      }
    }

    void Slice::set_inter(const Interval& y)
    {
      assert(m_codomain.intersects(y));
      m_codomain &= y;

      *m_input_gate &= m_codomain;
      if(prev_slice() != NULL)
        *m_input_gate &= prev_slice()->codomain();

      *m_output_gate &= m_codomain;
      if(next_slice() != NULL)
        *m_output_gate &= next_slice()->codomain();

      if(m_synthesis_reference != NULL)
      {
        m_synthesis_reference->request_values_update();
        m_synthesis_reference->request_integrals_update();
      }
    }

    void Slice::set_input_gate(const Interval& input_gate)
    {
      *m_input_gate = input_gate;
      *m_input_gate &= m_codomain;

      if(prev_slice() != NULL)
        *m_input_gate &= prev_slice()->codomain();

      if(m_synthesis_reference != NULL)
      {
        m_synthesis_reference->request_values_update();
        // Note: integrals are not impacted by gates
      }
    }

    void Slice::set_output_gate(const Interval& output_gate)
    {
      *m_output_gate = output_gate;
      *m_output_gate &= m_codomain;

      if(next_slice() != NULL)
        *m_output_gate &= next_slice()->codomain();

      if(m_synthesis_reference != NULL)
      {
        m_synthesis_reference->request_values_update();
        // Note: integrals are not impacted by gates
      }
    }
    
    const Slice& Slice::inflate(double rad)
    {
      assert(rad >= 0.);
      
      Interval e(-rad,rad);
      set_envelope(m_codomain + e);
      set_input_gate(*m_input_gate + e);
      set_output_gate(*m_output_gate + e);
    }

    void Slice::shift_domain(double shift_ref)
    {
      set_domain(m_domain - shift_ref);

      if(m_synthesis_reference != NULL)
      {
        // todo: update domain structure
      }
    }
    
    // String

    ostream& operator<<(ostream& str, const Slice& x)
    {
      str << "Slice " << x.domain() << "↦(" << x.input_gate() << ")" << x.codomain() << "(" << x.output_gate() << ")" << flush;
      return str;
    }


  // Protected methods
    
    void Slice::set_domain(const Interval& domain)
    {
      assert(valid_domain(domain));
      m_domain = domain;
    }

    // Slices structure

    void Slice::chain_slices(Slice *first_slice, Slice *second_slice)
    {
      if(first_slice != NULL)
        first_slice->m_next_slice = second_slice;

      if(second_slice != NULL)
        second_slice->m_prev_slice = first_slice;

      if(first_slice != NULL && second_slice != NULL)
      {
        if(second_slice->m_input_gate != NULL)
          *first_slice->m_output_gate &= *second_slice->m_input_gate;
        second_slice->m_input_gate = first_slice->m_output_gate;
      }
    }

    // Access values

    const IntervalVector Slice::codomain_box() const
    {
      return IntervalVector(m_codomain);
    }

    // Setting values
    
}