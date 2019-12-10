/**
 *  tubex-lib - Examples
 *  Robot evolution equations: causal chain
 * ----------------------------------------------------------------------------
 *
 *  \brief      Example from the paper "Guaranteed Computation of Robot Trajectories"
 *              Simon Rohou, Luc Jaulin, Lyudmila Mihaylova, Fabrice Le Bars, Sandor M. Veres
 *
 *  \date       2016
 *  \author     Simon Rohou
 *  \copyright  Copyright 2019 Simon Rohou
 *  \license    This program is distributed under the terms of
 *              the GNU Lesser General Public License (LGPL).
 */

#include "ibex.h"
#include "tubex.h"
#include "tubex-robotics.h"

using namespace std;
using namespace ibex;
using namespace tubex;

#define FINAL_CONDITION true

int main()
{
  /* =========== INITIALIZATION =========== */

    clock_t t_start = clock();
    double tf = 5.0;
    Interval domain(0.0, tf);
    double timestep = 0.1;

    // Domaine etat
    TubeVector x(domain, timestep, 6);
    x[0] = Tube(x[0], Interval(-999.9, 999.9)); // position x initial
    x[1] = Tube(x[1], Interval(-999.9, 999.9)); // position y initial
    x[2] = Tube(x[2], Interval(-1.0, 1.0));     // cos(th) initial
    x[3] = Tube(x[3], Interval(-1.0, 1.0));     // sin(th) initial
    x[4] = Tube(x[4], Interval(-1.0, 1.0));       // vx initial
    x[5] = Tube(x[5], Interval(-1.0, 1.0));       // vy initial

    // Conditions initiales etat
    IntervalVector x0(x.size());
    x0[0] = Interval(1.0); // position x initial
    x0[1] = Interval(0.0); // position y initial
    x0[2] = Interval(0.0); // cos(th) initial
    x0[3] = Interval(1.0); // sin(th) initial
    x0[4] = Interval(0.0); // vx initial
    x0[5] = Interval(1.0); // vy initial
    x.set(x0, 0.0);

    // Domaine de mesures
    TubeVector y(domain, timestep, 3);
    y[0] = Tube(y[0], Interval(-20, 20));  // acceleration axe i
    y[1] = Tube(y[1], Interval(-20, 20));  // acceleration axe j
    y[2] = Tube(y[2], Interval(-20, 20));  // vitesse de rotation

    // Trajectoire = cercle de rayon 1
    Trajectory tx(domain, tubex::Function("cos(t)"));
    Trajectory ty(domain, tubex::Function("sin(t)"));
    Trajectory tct(domain, tubex::Function("cos(t+3.14159265358979323846/2.0)"));
    Trajectory tst(domain, tubex::Function("sin(t+3.14159265358979323846/2.0)"));


  /* =========== CREATING CONTRACTORS =========== */

    tubex::CtcDeriv* ctc_dyn = new tubex::CtcDeriv();
    tubex::CtcFwdBwd* ctc_fwdbwd = new tubex::CtcFwdBwd(
      tubex::Function("x", "y", "cos_th", "sin_th", "vx", "vy",
                      "(cos_th^2 + sin_th^2 - 1.0; \
                        vx^2 + vy^2 - 1.0)"));

  /* =========== MESURES =========== */

    double dt = 0.02;
    for (double t=dt; t<tf; t=t+dt) {
      // Mesures
      double th = t+M_PI/2.0;
      double ax = -cos(t);
      double ay = -sin(t);
      // Mesures IMU à 50Hz
      y[0].set(Interval( ax*cos(th) + ay*sin(th)).inflate(0.01), Interval(t-dt, t)); // acceleration axe i
      y[1].set(Interval(-ax*sin(th) + ay*cos(th)).inflate(0.01), Interval(t-dt, t)); // acceleration axe j
      y[2].set(Interval(1.0).inflate(0.0001), Interval(t-dt, t));                      // vitesse de rotation
      // Mesure de position à chaque cercle complet
      if ((int(t*100) % (2*314)) == 0) {
        cout << "Mesure position \n";
        x[0].set(Interval(1.0).inflate(0.01), t);
        x[1].set(Interval(0.0).inflate(0.01), t);
        x[2].set(Interval(0.0).inflate(0.01), t);
        x[3].set(Interval(1.0).inflate(0.01), t);
      }
      // resample etat
      x.sample(t-dt);
      x.sample(t);
      // MAJ dynamique à 10Hz
      if ((int(t*100) % 10) == 0) {
        // Propagation
        double t0 = std::max(t - 2.0, t0);
        double t1 = std::min(t + timestep, tf);
        Interval domain(t0, t1);
        // Contrainte sur l'etat
        // ctc_fwdbwd->set_restricted_domain(domain);
        ctc_fwdbwd->contract(x);
        // COntrainte en dynamique
        // ctc_dyn->set_restricted_domain(domain);
        ctc_dyn->contract(x[0],  x[4]);
        ctc_dyn->contract(x[1],  x[5]);
        ctc_dyn->contract(x[2], -y[2]*x[3]);
        ctc_dyn->contract(x[3],  y[2]*x[2]);
        ctc_dyn->contract(x[4],  y[0]*x[2] - y[1]*x[3]);
        ctc_dyn->contract(x[5],  y[0]*x[3] + y[1]*x[2]);
      }
    }

  /* =========== GRAPHICS =========== */
    vibes::beginDrawing();
    VIBesFigTube fig_x("Tube x", &x[0]);
    fig_x.add_trajectory(&tx, "traj_x");
    fig_x.set_properties(100, 100, 600, 300);
    fig_x.show();
    VIBesFigTube fig_y("Tube y", &x[1]);
    fig_y.add_trajectory(&ty, "traj_y");
    fig_y.set_properties(100, 100, 600, 300);
    fig_y.show();
    VIBesFigTube fig_ct("Tube cos_th", &x[2]);
    fig_ct.add_trajectory(&tct, "traj_ct");
    fig_ct.set_properties(100, 100, 600, 300);
    fig_ct.show();
    VIBesFigTube fig_st("Tube sin_th", &x[3]);
    fig_st.add_trajectory(&tst, "traj_st");
    fig_st.set_properties(100, 100, 600, 300);
    fig_st.show();
    VIBesFigTube fig_cmd_1("Tube cmd_1", &y[0]);
    fig_cmd_1.set_properties(100, 100, 600, 300);
    fig_cmd_1.show();
    VIBesFigTube fig_cmd_2("Tube cmd_2", &y[1]);
    fig_cmd_2.set_properties(100, 100, 600, 300);
    fig_cmd_2.show();
    VIBesFigTube fig_cmd_3("Tube cmd_3", &y[2]);
    fig_cmd_3.set_properties(100, 100, 600, 300);
    fig_cmd_3.show();
    vibes::endDrawing();

    printf("Time taken: %.2fs\n", (double)(clock() - t_start)/CLOCKS_PER_SEC);

    delete ctc_dyn;
    delete ctc_fwdbwd;
}
