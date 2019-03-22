/** 
 *  tubex-lib - Examples
 *  Testing micro contractors
 * ----------------------------------------------------------------------------
 *
 *  \date       2019
 *  \author     Simon Rohou
 *  \copyright  Copyright 2019 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include "tubex.h"
#include "ibex_Function.h"
#include "ibex_CtcFwdBwd.h"

using namespace std;
using namespace ibex;
using namespace tubex;

// This method applies a CtcFwdBwd on slices
void ctc_fwdbwd_contract(Slice &x, Slice &v, ibex::CtcFwdBwd& ctc_fwdbwd);

int main(int argc, char *argv[])
{
  clock_t t_start;

  // Tubes
  double timestep = 0.01;
  Interval domain(0.,10.); // temporal definition domain
  TubeVector x(domain, timestep, 2), xold(x);

  // Pointers to slices
  Slice *sx, *sv; // slices are one-D: one for x, one for xdot (v)
  Slice sx_old(domain), sv_old(domain); // used for fixed point

  // Defining contractors
  CtcDeriv ctc_deriv;

  // Method 1: fixed point on slices, then propagation
  {
    // FwdBwd contractor for slices
    ibex::Function *f = new ibex::Function("x", "v", "v+sin(x)");
    ibex::CtcFwdBwd uctc_fwdbwd(*f);

    // Restoring tubes
    x = TubeVector(domain, timestep, 2); // 2d TubeVector: x,xdot
    x[0].set(Interval(1.), 0.); // initial condition

    t_start = clock();

    do
    {
      xold = x;
      sx = x[0].first_slice(); sv = x[1].first_slice();

      while(sx != NULL)
      {
        do
        {
          sx_old = *sx; sv_old = *sv;
          ctc_fwdbwd_contract(*sx, *sv, uctc_fwdbwd); // Cfwdbwd
          ctc_deriv.contract(*sx, *sv);               // Cd/dt
        } while(sx_old != *sx || sv_old != *sv);      // fixed point

        sx = sx->next_slice(); sv = sv->next_slice(); // slices iteration
      }
    } while(xold != x);                               // fixed point
    
    cout << "[method 1] Volume of tube x = " << x.volume() << endl;
    printf("           Time taken: %.2fs\n", (double)(clock() - t_start)/CLOCKS_PER_SEC);
    delete f;
  }

  // Method 2: global contractions only, up to a fixed point
  {
    // FwdBwd contractor for tubes
    tubex::Function *f = new tubex::Function("x", "v", "v+sin(x)");
    tubex::CtcFwdBwd ctc_fwdbwd(*f);

    // Restoring tubes
    x = TubeVector(domain, timestep, 2); // 2d TubeVector: x,xdot
    x[0].set(Interval(1.), 0.); // initial condition

    t_start = clock();

    do
    {
      xold = x;
      ctc_fwdbwd.contract(x);                         // Cfwdbwd
      ctc_deriv.contract(x[0], x[1]);                 // Cd/dt

    } while(xold != x);                               // fixed point

    cout << "[method 2] Volume of tube x = " << x.volume() << endl;
    printf("           Time taken: %.2fs\n", (double)(clock() - t_start)/CLOCKS_PER_SEC);
    delete f;
  }

  // Display
  vibes::beginDrawing();
  Trajectory truth(domain, tubex::Function("2.*atan(exp(-t)*tan(0.5))"));
  VIBesFigTube fig_x("x", &x[0], &truth);
  fig_x.set_properties(100, 100, 600, 600);
  fig_x.show();
  vibes::endDrawing();

  return EXIT_SUCCESS;
}

void ctc_fwdbwd_contract(Slice &x, Slice &v, ibex::CtcFwdBwd& ctc_fwdbwd)
{
  IntervalVector envelope(2);
  envelope[0] = x.codomain(); envelope[1] = v.codomain();
  ctc_fwdbwd.contract(envelope);
  x.set_envelope(envelope[0]); v.set(envelope[1]);

  IntervalVector ingate(2);
  ingate[0] = x.input_gate(); ingate[1] = v.input_gate();
  ctc_fwdbwd.contract(ingate);
  x.set_input_gate(ingate[0]); v.set_input_gate(ingate[1]);

  IntervalVector outgate(2);
  outgate[0] = x.output_gate(); outgate[1] = v.output_gate();
  ctc_fwdbwd.contract(outgate);
  x.set_output_gate(outgate[0]); v.set_output_gate(outgate[1]);
}