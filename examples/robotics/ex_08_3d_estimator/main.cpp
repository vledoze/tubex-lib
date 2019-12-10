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
#define G               9.81
#define D2R             3.14/180.0

int main()
{
  /* =========== INITIALIZATION =========== */

    double tf = 10.0;
    Interval domain(0.0, tf);
    double timestep = 0.01;

    // Creating tube over the [0,tf] domain with some timestep:
    uint8_t i = 0;
    uint8_t pos_e = i; i++;
    uint8_t pos_n = i; i++;
    uint8_t pos_u = i; i++;
    TubeVector dddx(domain, i);
    // TubeVector dddx(domain, timestep, i);
    uint8_t cos_psi = i; i++;
    uint8_t sin_psi = i; i++;
    uint8_t cos_tht = i; i++;
    uint8_t sin_tht = i; i++;
    uint8_t cos_phi = i; i++;
    uint8_t sin_phi = i; i++;
    TubeVector x(domain, i);
    TubeVector dx(domain, i);
    TubeVector ddx(domain, i);
    // TubeVector x(domain, timestep, i);
    // TubeVector dx(domain, timestep, i);
    // TubeVector ddx(domain, timestep, i);

    // Domaine de l'etat
    x[pos_e] = Tube(x[pos_e], Interval(-999.0, 999.0));
    x[pos_n] = Tube(x[pos_n], Interval(-999.0, 999.0));
    x[pos_u] = Tube(x[pos_u], Interval(-999.0, 999.0));
    x[cos_psi] = Tube(x[cos_psi], Interval(-1.0, 1.0));
    x[sin_psi] = Tube(x[sin_psi], Interval(-1.0, 1.0));
    x[cos_tht] = Tube(x[cos_tht], Interval(-1.0, 1.0));
    x[sin_tht] = Tube(x[sin_tht], Interval(-1.0, 1.0));
    x[cos_phi] = Tube(x[cos_phi], Interval(-1.0, 1.0));
    x[sin_phi] = Tube(x[sin_phi], Interval(-1.0, 1.0));

    // Domaine de la derivée
    dx[pos_e] = Tube(dx[pos_e], Interval(-20.0, 20.0));
    dx[pos_n] = Tube(dx[pos_n], Interval(-20.0, 20.0));
    dx[pos_u] = Tube(dx[pos_u], Interval(-20.0, 20.0));
    dx[cos_psi] = Tube(dx[cos_psi], Interval(-1.0, 1.0));
    dx[sin_psi] = Tube(dx[sin_psi], Interval(-1.0, 1.0));
    dx[cos_tht] = Tube(dx[cos_tht], Interval(-1.0, 1.0));
    dx[sin_tht] = Tube(dx[sin_tht], Interval(-1.0, 1.0));
    dx[cos_phi] = Tube(dx[cos_phi], Interval(-1.0, 1.0));
    dx[sin_phi] = Tube(dx[sin_phi], Interval(-1.0, 1.0));

    // Domaine de la derivée seconde
    ddx[pos_e] = Tube(ddx[pos_e], Interval(-2.0*G, 2.0*G));
    ddx[pos_n] = Tube(ddx[pos_n], Interval(-2.0*G, 2.0*G));
    ddx[pos_u] = Tube(ddx[pos_u], Interval(-2.0*G, 2.0*G));
    ddx[cos_psi] = Tube(ddx[cos_psi], Interval(-1.0, 1.0));
    ddx[sin_psi] = Tube(ddx[sin_psi], Interval(-1.0, 1.0));
    ddx[cos_tht] = Tube(ddx[cos_tht], Interval(-1.0, 1.0));
    ddx[sin_tht] = Tube(ddx[sin_tht], Interval(-1.0, 1.0));
    ddx[cos_phi] = Tube(ddx[cos_phi], Interval(-1.0, 1.0));
    ddx[sin_phi] = Tube(ddx[sin_phi], Interval(-1.0, 1.0));

    // Domaine de la derivée troisieme
    dddx[pos_e] = Tube(dddx[pos_e], Interval(-1.0, 1.0));
    dddx[pos_n] = Tube(dddx[pos_n], Interval(-1.0, 1.0));
    dddx[pos_u] = Tube(dddx[pos_u], Interval(-1.0, 1.0));


  /* =========== CREATING FUNCTIONS ===========*/
    ibex::Function f_matR("fun_matR.txt");
    ibex::Function f_matT("fun_matT.txt");
    ibex::Function f_matTi("fun_matTi.txt");


  /* =========== CREATING CONTRACTORS =========== */

    CtcDeriv ctc_dyn;
    ctc_dyn.set_fast_mode(true);

    bool with_eval = true;
    CtcEval ctc_eval;
    ctc_eval.set_fast_mode(true);
    // ctc_eval.enable_temporal_propagation(false);
    // ctc_eval.preserve_slicing(true);

    tubex::CtcFwdBwd ctc_state(tubex::Function(x.size(), "fctc_state.txt"));
    ctc_state.set_fast_mode(true);

    ibex::CtcFwdBwd* ctc_eul1 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_eul/ctc_eul1.txt")));
    ibex::CtcFwdBwd* ctc_eul2 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_eul/ctc_eul2.txt")));
    ibex::CtcFwdBwd* ctc_eul3 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_eul/ctc_eul3.txt")));
    ibex::CtcFwdBwd* ctc_matR1 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_matR/ctc_matR1.txt")));
    ibex::CtcFwdBwd* ctc_matR2 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_matR/ctc_matR2.txt")));
    ibex::CtcFwdBwd* ctc_matR3 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_matR/ctc_matR3.txt")));
    ibex::CtcFwdBwd* ctc_matR4 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_matR/ctc_matR4.txt")));
    ibex::CtcFwdBwd* ctc_matR5 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_matR/ctc_matR5.txt")));
    ibex::CtcFwdBwd* ctc_matR6 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_matR/ctc_matR6.txt")));
    ibex::CtcFwdBwd* ctc_mag1 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_mag/ctc_mag1.txt")));
    ibex::CtcFwdBwd* ctc_mag2 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_mag/ctc_mag2.txt")));
    ibex::CtcFwdBwd* ctc_mag3 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_mag/ctc_mag3.txt")));
    ibex::CtcFwdBwd* ctc_mag4 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_mag/ctc_mag4.txt")));
    ibex::CtcFwdBwd* ctc_mag5 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_mag/ctc_mag5.txt")));
    ibex::CtcFwdBwd* ctc_mag6 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_mag/ctc_mag6.txt")));
    ibex::CtcFwdBwd* ctc_acc1 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_acc/ctc_acc1.txt")));
    ibex::CtcFwdBwd* ctc_acc2 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_acc/ctc_acc2.txt")));
    ibex::CtcFwdBwd* ctc_acc3 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_acc/ctc_acc3.txt")));
    ibex::CtcFwdBwd* ctc_acc4 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_acc/ctc_acc4.txt")));
    ibex::CtcFwdBwd* ctc_acc5 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_acc/ctc_acc5.txt")));
    ibex::CtcFwdBwd* ctc_acc6 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_acc/ctc_acc6.txt")));
    ibex::CtcFwdBwd* ctc_vispos1 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_vispos/ctc_vispos1.txt")));
    ibex::CtcFwdBwd* ctc_vispos2 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_vispos/ctc_vispos2.txt")));
    ibex::CtcFwdBwd* ctc_vispos3 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_vispos/ctc_vispos3.txt")));
    ibex::CtcFwdBwd* ctc_vispos4 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_vispos/ctc_vispos4.txt")));
    ibex::CtcFwdBwd* ctc_vispos5 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_vispos/ctc_vispos5.txt")));
    ibex::CtcFwdBwd* ctc_vispos6 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_vispos/ctc_vispos6.txt")));
    ibex::CtcFwdBwd* ctc_gps1 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_gps/ctc_gps1.txt")));
    ibex::CtcFwdBwd* ctc_gps2 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_gps/ctc_gps2.txt")));
    ibex::CtcFwdBwd* ctc_gps3 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_gps/ctc_gps3.txt")));
    ibex::CtcFwdBwd* ctc_fluxopt1 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_fluxopt/ctc_fluxopt1.txt")));
    ibex::CtcFwdBwd* ctc_fluxopt2 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_fluxopt/ctc_fluxopt2.txt")));
    ibex::CtcFwdBwd* ctc_fluxopt3 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_fluxopt/ctc_fluxopt3.txt")));
    ibex::CtcFwdBwd* ctc_fluxopt4 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_fluxopt/ctc_fluxopt4.txt")));
    ibex::CtcFwdBwd* ctc_fluxopt5 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_fluxopt/ctc_fluxopt5.txt")));
    ibex::CtcFwdBwd* ctc_fluxopt6 = new ibex::CtcFwdBwd(*(new ibex::NumConstraint("ctc_fluxopt/ctc_fluxopt6.txt")));

    vector<ibex::Ctc*> v_ctc_mag;
    v_ctc_mag.push_back(ctc_mag1);
    v_ctc_mag.push_back(ctc_mag2);
    v_ctc_mag.push_back(ctc_mag3);
    v_ctc_mag.push_back(ctc_mag4);
    v_ctc_mag.push_back(ctc_mag5);
    v_ctc_mag.push_back(ctc_mag6);
    v_ctc_mag.push_back(ctc_eul1);
    v_ctc_mag.push_back(ctc_eul2);
    v_ctc_mag.push_back(ctc_eul3);
    v_ctc_mag.push_back(ctc_matR1);
    v_ctc_mag.push_back(ctc_matR2);
    v_ctc_mag.push_back(ctc_matR3);
    v_ctc_mag.push_back(ctc_matR4);
    v_ctc_mag.push_back(ctc_matR5);
    v_ctc_mag.push_back(ctc_matR6);
    ibex::CtcCompo ctccompo_mag(v_ctc_mag);
    ibex::CtcFixPoint ctcfix_mag(ctccompo_mag);
    IntervalVector vars_ctc_mag(9);

    vector<ibex::Ctc*> v_ctc_acc;
    v_ctc_acc.push_back(ctc_acc1);
    v_ctc_acc.push_back(ctc_acc2);
    v_ctc_acc.push_back(ctc_acc3);
    v_ctc_acc.push_back(ctc_acc4);
    v_ctc_acc.push_back(ctc_acc5);
    v_ctc_acc.push_back(ctc_acc6);
    v_ctc_acc.push_back(ctc_eul1);
    v_ctc_acc.push_back(ctc_eul2);
    v_ctc_acc.push_back(ctc_eul3);
    v_ctc_acc.push_back(ctc_matR1);
    v_ctc_acc.push_back(ctc_matR2);
    v_ctc_acc.push_back(ctc_matR3);
    v_ctc_acc.push_back(ctc_matR4);
    v_ctc_acc.push_back(ctc_matR5);
    v_ctc_acc.push_back(ctc_matR6);
    ibex::CtcCompo ctccompo_acc(v_ctc_acc);
    ibex::CtcFixPoint ctcfix_acc(ctccompo_acc);
    IntervalVector vars_ctc_acc(9);

    vector<ibex::Ctc*> v_ctc_vispos;
    v_ctc_vispos.push_back(ctc_vispos1);
    v_ctc_vispos.push_back(ctc_vispos2);
    v_ctc_vispos.push_back(ctc_vispos3);
    v_ctc_vispos.push_back(ctc_vispos4);
    v_ctc_vispos.push_back(ctc_vispos5);
    v_ctc_vispos.push_back(ctc_vispos6);
    ibex::CtcCompo ctccompo_vispos(v_ctc_vispos);
    ibex::CtcFixPoint ctcfix_vispos(ctccompo_vispos);
    IntervalVector vars_ctc_vispos(15);
    bool first_mes_vispos = true;

    vector<ibex::Ctc*> v_ctc_fluxopt;
    v_ctc_fluxopt.push_back(ctc_fluxopt1);
    v_ctc_fluxopt.push_back(ctc_fluxopt2);
    v_ctc_fluxopt.push_back(ctc_fluxopt3);
    v_ctc_fluxopt.push_back(ctc_fluxopt4);
    v_ctc_fluxopt.push_back(ctc_fluxopt5);
    v_ctc_fluxopt.push_back(ctc_fluxopt6);
    ibex::CtcCompo ctccompo_fluxopt(v_ctc_fluxopt);
    ibex::CtcFixPoint ctcfix_fluxopt(ctccompo_fluxopt);
    IntervalVector vars_ctc_fluxopt(12);

    vector<ibex::Ctc*> v_ctc_gps;
    v_ctc_gps.push_back(ctc_gps1);
    v_ctc_gps.push_back(ctc_gps2);
    v_ctc_gps.push_back(ctc_gps3);
    ibex::CtcCompo ctccompo_gps(v_ctc_gps);
    ibex::CtcFixPoint ctcfix_gps(ctccompo_gps);
    IntervalVector vars_ctc_gps(9);
    bool first_mes_gps = true;

  /* =========== CONDITIONS INITIALES =========== */

    ctc_eval.contract(0.0, Interval(0.0), x[pos_e], dx[pos_e]);     // position e initiale
    ctc_eval.contract(0.0, Interval(0.0), x[pos_n], dx[pos_n]);     // position n initiale
    ctc_eval.contract(0.0, Interval(0.0), x[pos_u], dx[pos_u]);     // position u initiale
    ctc_eval.contract(0.0, Interval(0.0), dx[pos_e], ddx[pos_e]);   // vitesse e initiale
    ctc_eval.contract(0.0, Interval(0.0), dx[pos_n], ddx[pos_n]);   // vitesse n initiale
    ctc_eval.contract(0.0, Interval(0.0), dx[pos_u], ddx[pos_u]);   // vitesse u initiale
    ctc_eval.contract(0.0, Interval(0.0), ddx[pos_e], dddx[pos_e]); // acc e initiale
    ctc_eval.contract(0.0, Interval(0.0), ddx[pos_n], dddx[pos_n]); // acc n initiale
    ctc_eval.contract(0.0, Interval(0.0), ddx[pos_u], dddx[pos_u]); // acc u initiale
    ctc_eval.contract(0.0, Interval(-1.0, 1.0),                        x[cos_psi], dx[cos_psi]);
    ctc_eval.contract(0.0, Interval(-1.0, 1.0),                        x[sin_psi], dx[sin_psi]);
    ctc_eval.contract(0.0, ibex::cos(Interval(0.0).inflate(10.0*D2R)), x[cos_tht], dx[cos_tht]);
    ctc_eval.contract(0.0, ibex::sin(Interval(0.0).inflate(10.0*D2R)), x[sin_tht], dx[sin_tht]);
    ctc_eval.contract(0.0, ibex::cos(Interval(0.0).inflate(10.0*D2R)), x[cos_phi], dx[cos_phi]);
    ctc_eval.contract(0.0, ibex::sin(Interval(0.0).inflate(10.0*D2R)), x[sin_phi], dx[sin_phi]);
    ctc_eval.contract(0.0, Interval(0.0), dx[cos_psi], ddx[cos_psi]);
    ctc_eval.contract(0.0, Interval(0.0), dx[sin_psi], ddx[sin_psi]);
    ctc_eval.contract(0.0, Interval(0.0), dx[cos_tht], ddx[cos_tht]);
    ctc_eval.contract(0.0, Interval(0.0), dx[sin_tht], ddx[sin_tht]);
    ctc_eval.contract(0.0, Interval(0.0), dx[cos_phi], ddx[cos_phi]);
    ctc_eval.contract(0.0, Interval(0.0), dx[sin_phi], ddx[sin_phi]);

    ctc_dyn.contract(ddx[pos_e],  dddx[pos_e]);
    ctc_dyn.contract(ddx[pos_n],  dddx[pos_n]);
    ctc_dyn.contract(ddx[pos_u],  dddx[pos_u]);
    ctc_dyn.contract(dx,  ddx);
    ctc_dyn.contract(x,  dx);

  /* =========== MESURES =========== */

    double t;
    int dtms_imu = 20; // Mesures IMU à 50Hz (toutes les 20ms)
    int dtms_mag = 40; // Mesures IMU a 25Hz (toutes les 40ms)
    int dtms_vispos = 200; // Mesures Visodom à 5Hz (toutes les 200ms)
    int dtms_gps = 1000; // Mesures GPS à 1Hz (toutes les 1000ms)
    int dtms_fluxopt = 1000; // Update Vitesse à 1Hz (toutes les 1000ms)
    int dtms_dyn = 100; // Update dynamique à 10Hz (toutes les 100ms)
    // debug
    clock_t t_start = clock();
    clock_t t_start_imu;
    clock_t t_start_mag;
    clock_t t_start_dyn;
    // boucle principale
    for (int tms=1; tms<=int(tf*1000); tms++) {
      // Instant courant
      t = double(tms)/1000.0;
      // Horizon temporel de propagation
      double t0 = t - 1.0;
      if(t0 < 0.0) t0 = 0.0;
      double t1 = t + 0.01;
      if(t1 > tf) t1 = tf;
      // Test
      x.del_first_slices(t0);
      dx.del_first_slices(t0);
      ddx.del_first_slices(t0);
      dddx.del_first_slices(t0);
      // Mesures IMU
      if ((tms%dtms_imu) == 0)
      {
        // t_start_imu = clock();
        // Horizon de validité des mesures
        double t0_imu = double(tms-dtms_imu)/1000.0;
        if(t0_imu < t0) t0_imu = t0;
        Interval domain = Interval(t0_imu, t);
        ctc_eval.set_restricted_domain(domain);
        // resampling
        x.sample(t);
        dx.sample(t);
        ddx.sample(t);
        dddx.sample(t);
        // Contraction par mesures acc pour orientation
        vars_ctc_acc[0] = Interval(0.0).inflate(0.01);
        vars_ctc_acc[1] = Interval(0.0).inflate(0.01);
        vars_ctc_acc[2] = Interval(G).inflate(0.01);
        vars_ctc_acc[3] = x[cos_psi].interpol(t, dx[cos_psi]);
        vars_ctc_acc[4] = x[sin_psi].interpol(t, dx[sin_psi]);
        vars_ctc_acc[5] = x[cos_tht].interpol(t, dx[cos_tht]);
        vars_ctc_acc[6] = x[sin_tht].interpol(t, dx[sin_tht]);
        vars_ctc_acc[7] = x[cos_phi].interpol(t, dx[cos_phi]);
        vars_ctc_acc[8] = x[sin_phi].interpol(t, dx[sin_phi]);
        ctccompo_acc.contract(vars_ctc_acc);
        // Renvoi des mesures
        if (with_eval) {
          ctc_eval.contract(t, vars_ctc_acc[3], x[cos_psi], dx[cos_psi]);
          ctc_eval.contract(t, vars_ctc_acc[4], x[sin_psi], dx[sin_psi]);
          ctc_eval.contract(t, vars_ctc_acc[5], x[cos_tht], dx[cos_tht]);
          ctc_eval.contract(t, vars_ctc_acc[6], x[sin_tht], dx[sin_tht]);
          ctc_eval.contract(t, vars_ctc_acc[7], x[cos_phi], dx[cos_phi]);
          ctc_eval.contract(t, vars_ctc_acc[8], x[sin_phi], dx[sin_phi]);
       }
       else {
          x[cos_psi].set(vars_ctc_acc[3], t);
          x[sin_psi].set(vars_ctc_acc[4], t);
          x[cos_tht].set(vars_ctc_acc[5], t);
          x[sin_tht].set(vars_ctc_acc[6], t);
          x[cos_phi].set(vars_ctc_acc[7], t);
          x[sin_phi].set(vars_ctc_acc[8], t);
        }
        // Prise en compte des mesures acc pour derivee de la vitesse
        IntervalVector vars_acc(12);
        vars_acc[0] = vars_ctc_acc[0];
        vars_acc[1] = vars_ctc_acc[1];
        vars_acc[2] = vars_ctc_acc[2];
        vars_acc[3] = x[cos_psi].interpol(t, dx[cos_psi]);
        vars_acc[4] = x[sin_psi].interpol(t, dx[sin_psi]);
        vars_acc[5] = x[cos_tht].interpol(t, dx[cos_tht]);
        vars_acc[6] = x[sin_tht].interpol(t, dx[sin_tht]);
        vars_acc[7] = x[cos_phi].interpol(t, dx[cos_phi]);
        vars_acc[8] = x[sin_phi].interpol(t, dx[sin_phi]);
        IntervalVector dx_vit_enu = f_matR.eval_vector(vars_acc);
        if (with_eval) {
          ctc_eval.contract(t, dx_vit_enu[0],   ddx[pos_e], dddx[pos_e]);
          ctc_eval.contract(t, dx_vit_enu[1],   ddx[pos_n], dddx[pos_n]);
          ctc_eval.contract(t, dx_vit_enu[2]-G, ddx[pos_u], dddx[pos_u]);
        }
        else {
          ddx[pos_e].set(dx_vit_enu[0],   t); // acceleration axe e
          ddx[pos_n].set(dx_vit_enu[1],   t); // acceleration axe n
          ddx[pos_u].set(dx_vit_enu[2]-G, t); // acceleration axe u
        }
        // Prise en compte des mesures gyr
        IntervalVector vars_gyr(12);
        vars_gyr[0] = Interval(0.0).inflate(0.005);
        vars_gyr[1] = Interval(0.0).inflate(0.005);
        vars_gyr[2] = Interval(0.0).inflate(0.005);
        vars_gyr[3] = x[cos_psi].interpol(t, dx[cos_psi]);
        vars_gyr[4] = x[sin_psi].interpol(t, dx[sin_psi]);
        vars_gyr[5] = x[cos_tht].interpol(t, dx[cos_tht]);
        vars_gyr[6] = x[sin_tht].interpol(t, dx[sin_tht]);
        vars_gyr[7] = x[cos_phi].interpol(t, dx[cos_phi]);
        vars_gyr[8] = x[sin_phi].interpol(t, dx[sin_phi]);
        IntervalVector dx_psithtphi = f_matT.eval_vector(vars_gyr);
        if (with_eval) {
          ctc_eval.contract(t, -dx_psithtphi[0]*x[sin_psi](t), dx[cos_psi], ddx[cos_psi]);
          ctc_eval.contract(t,  dx_psithtphi[0]*x[cos_psi](t), dx[sin_psi], ddx[sin_psi]);
          ctc_eval.contract(t, -dx_psithtphi[1]*x[sin_tht](t), dx[cos_tht], ddx[cos_tht]);
          ctc_eval.contract(t,  dx_psithtphi[1]*x[cos_tht](t), dx[sin_tht], ddx[sin_tht]);
          ctc_eval.contract(t, -dx_psithtphi[2]*x[sin_phi](t), dx[cos_phi], ddx[cos_phi]);
          ctc_eval.contract(t,  dx_psithtphi[2]*x[cos_phi](t), dx[sin_phi], ddx[sin_phi]);
        }
        else {
          dx[cos_psi].set(-dx_psithtphi[0]*(x[sin_psi].interpol(t, dx[cos_psi])), t);
          dx[sin_psi].set( dx_psithtphi[0]*(x[cos_psi].interpol(t, dx[sin_psi])), t);
          dx[cos_tht].set(-dx_psithtphi[1]*(x[sin_tht].interpol(t, dx[cos_tht])), t);
          dx[sin_tht].set( dx_psithtphi[1]*(x[cos_tht].interpol(t, dx[sin_tht])), t);
          dx[cos_phi].set(-dx_psithtphi[2]*(x[sin_phi].interpol(t, dx[cos_phi])), t);
          dx[sin_phi].set( dx_psithtphi[2]*(x[cos_phi].interpol(t, dx[sin_phi])), t);
        }
        // printf("Time taken IMU: %.2fs\n", (double)(clock() - t_start_imu)/CLOCKS_PER_SEC);
      }
      // Mesures MAG
      if ((tms%dtms_mag) == 0)
      {
        // t_start_mag = clock();
        // Horizon de validité des mesures
        double t0_mag = double(tms-dtms_mag)/1000.0;
        if(t0_mag < t0) t0_mag = t0;
        Interval domain = Interval(t0_mag, t);
        ctc_eval.set_restricted_domain(domain);
        // resampling
        x.sample(t);
        dx.sample(t);
        ddx.sample(t);
        // Contraction
        vars_ctc_mag[0] = Interval(-5).inflate(0.001);
        vars_ctc_mag[1] = Interval(20).inflate(0.001);
        vars_ctc_mag[2] = Interval(-40).inflate(0.001);
        vars_ctc_mag[3] = x[cos_psi].interpol(t, dx[cos_psi]);
        vars_ctc_mag[4] = x[sin_psi].interpol(t, dx[sin_psi]);
        vars_ctc_mag[5] = x[cos_tht].interpol(t, dx[cos_tht]);
        vars_ctc_mag[6] = x[sin_tht].interpol(t, dx[sin_tht]);
        vars_ctc_mag[7] = x[cos_phi].interpol(t, dx[cos_phi]);
        vars_ctc_mag[8] = x[sin_phi].interpol(t, dx[sin_phi]);
        ctccompo_mag.contract(vars_ctc_mag);
        // Renvoi des mesures
        if (with_eval) {
          ctc_eval.contract(t, vars_ctc_mag[3], x[cos_psi], dx[cos_psi]);
          ctc_eval.contract(t, vars_ctc_mag[4], x[sin_psi], dx[sin_psi]);
          ctc_eval.contract(t, vars_ctc_mag[5], x[cos_tht], dx[cos_tht]);
          ctc_eval.contract(t, vars_ctc_mag[6], x[sin_tht], dx[sin_tht]);
          ctc_eval.contract(t, vars_ctc_mag[7], x[cos_phi], dx[cos_phi]);
          ctc_eval.contract(t, vars_ctc_mag[8], x[sin_phi], dx[sin_phi]);
        }
        else {
          x[cos_psi].set(vars_ctc_mag[3], t);
          x[sin_psi].set(vars_ctc_mag[4], t);
          x[cos_tht].set(vars_ctc_mag[5], t);
          x[sin_tht].set(vars_ctc_mag[6], t);
          x[cos_phi].set(vars_ctc_mag[7], t);
          x[sin_phi].set(vars_ctc_mag[8], t);
        }
        // printf("Time taken MAG: %.2fs\n", (double)(clock() - t_start_mag)/CLOCKS_PER_SEC);
      }
      // Mesures Vispos
      if ((tms%dtms_vispos) == 0)
      {
        // Horizon de validité des mesures
        double t0_vispos = double(tms-dtms_vispos)/1000.0;
        if(t0_vispos < t0) t0_vispos = t0;
        ctc_eval.set_restricted_domain(Interval(t0_vispos, t));
        // Contraction
        vars_ctc_vispos[0] = Interval(0.0).inflate(1.0);
        vars_ctc_vispos[1] = Interval(0.0).inflate(1.0);
        vars_ctc_vispos[2] = Interval(0.0).inflate(1.0);
        vars_ctc_vispos[3] = x[pos_e].interpol(t, dx[pos_e]);
        vars_ctc_vispos[4] = x[pos_n].interpol(t, dx[pos_n]);
        vars_ctc_vispos[5] = x[pos_u].interpol(t, dx[pos_u]);
        if (first_mes_vispos) {
          // Domaine de contraction
          Interval domain(t0, t1);
          // Dynamique de positions
          ctc_dyn.set_restricted_domain(domain);
          ctc_dyn.contract(ddx[pos_e],  dddx[pos_e]);
          ctc_dyn.contract(ddx[pos_n],  dddx[pos_n]);
          ctc_dyn.contract(ddx[pos_u],  dddx[pos_u]);
          ctc_dyn.contract(dx,  ddx);
          ctc_dyn.contract(x,  dx);
          // Position de depart
          vars_ctc_vispos[6] = x[cos_psi].interpol(t, dx[cos_psi]);
          vars_ctc_vispos[7] = x[sin_psi].interpol(t, dx[sin_psi]);
          vars_ctc_vispos[8] = x[cos_tht].interpol(t, dx[cos_tht]);
          vars_ctc_vispos[9] = x[sin_tht].interpol(t, dx[sin_tht]);
          vars_ctc_vispos[10] = x[cos_phi].interpol(t, dx[cos_phi]);
          vars_ctc_vispos[11] = x[sin_phi].interpol(t, dx[sin_phi]);
          vars_ctc_vispos[12] = x[pos_e].interpol(t, dx[pos_e]);
          vars_ctc_vispos[13] = x[pos_n].interpol(t, dx[pos_n]);
          vars_ctc_vispos[14] = x[pos_u].interpol(t, dx[pos_u]);
          first_mes_vispos = false;
        }
        ctcfix_vispos.contract(vars_ctc_vispos);
        // Renvoi des mesures
        ctc_eval.contract(t, vars_ctc_vispos[3], x[pos_e], dx[pos_e]);
        ctc_eval.contract(t, vars_ctc_vispos[4], x[pos_n], dx[pos_n]);
        ctc_eval.contract(t, vars_ctc_vispos[5], x[pos_u], dx[pos_u]);
      }
      // Mesures GPS
      if ((tms%dtms_gps) == 0)
      {
        // Horizon de validité des mesures
        double t0_gps = double(tms-dtms_gps)/1000.0;
        if(t0_gps < t0) t0_gps = t0;
        ctc_eval.set_restricted_domain(Interval(t0_gps, t));
        // Contraction
        vars_ctc_gps[0] = Interval(50.0*D2R).inflate(1.0*D2R);
        vars_ctc_gps[1] = Interval(10.0*D2R).inflate(1.0*D2R);
        vars_ctc_gps[2] = Interval(50.0).inflate(5.0);
        vars_ctc_gps[3] = x[pos_e].interpol(t, dx[pos_e]);
        vars_ctc_gps[4] = x[pos_n].interpol(t, dx[pos_n]);
        vars_ctc_gps[5] = x[pos_u].interpol(t, dx[pos_u]);
        if (first_mes_gps) {
          // Position de depart
          vars_ctc_gps[6] = vars_ctc_gps[0];
          vars_ctc_gps[7] = vars_ctc_gps[1];
          vars_ctc_gps[8] = vars_ctc_gps[2];
          first_mes_gps = false;
        }
        ctcfix_gps.contract(vars_ctc_gps);
        // Renvoi des mesures
        ctc_eval.contract(t, vars_ctc_gps[3], x[pos_e], dx[pos_e]);
        ctc_eval.contract(t, vars_ctc_gps[4], x[pos_n], dx[pos_n]);
        ctc_eval.contract(t, vars_ctc_gps[5], x[pos_u], dx[pos_u]);
      }
      // Mesures Vitesse
      if ((tms%dtms_fluxopt) == 0)
      {
        // Horizon de validité des mesures
        double t0_fluxopt = double(tms-dtms_fluxopt)/1000.0;
        if(t0_fluxopt < t0) t0_fluxopt = t0;
        ctc_eval.set_restricted_domain(Interval(t0_fluxopt, t));
        // Contraction
        vars_ctc_fluxopt[0] = Interval(0.0).inflate(0.2);
        vars_ctc_fluxopt[1] = Interval(0.0).inflate(0.2);
        vars_ctc_fluxopt[2] = Interval(0.0).inflate(0.2);
        vars_ctc_fluxopt[3] = x[cos_psi].interpol(t, dx[cos_psi]);
        vars_ctc_fluxopt[4] = x[sin_psi].interpol(t, dx[sin_psi]);
        vars_ctc_fluxopt[5] = x[cos_tht].interpol(t, dx[cos_tht]);
        vars_ctc_fluxopt[6] = x[sin_tht].interpol(t, dx[sin_tht]);
        vars_ctc_fluxopt[7] = x[cos_phi].interpol(t, dx[cos_phi]);
        vars_ctc_fluxopt[8] = x[sin_phi].interpol(t, dx[sin_phi]);
        vars_ctc_fluxopt[9]  = dx[pos_e].interpol(t, ddx[pos_e]);
        vars_ctc_fluxopt[10] = dx[pos_n].interpol(t, ddx[pos_n]);
        vars_ctc_fluxopt[11] = dx[pos_u].interpol(t, ddx[pos_u]);
        ctcfix_fluxopt.contract(vars_ctc_fluxopt);
        // Renvoi des mesures
        ctc_eval.contract(t, vars_ctc_fluxopt[9] , dx[pos_e], ddx[pos_e]);
        ctc_eval.contract(t, vars_ctc_fluxopt[10], dx[pos_n], ddx[pos_n]);
        ctc_eval.contract(t, vars_ctc_fluxopt[11], dx[pos_u], ddx[pos_u]);
      }
      // MAJ dynamique
      if ((tms%dtms_dyn) == 0) {
        #ifdef DEBUG
        t_start_dyn = clock();
        #endif
        // Domaine de contraction
        Interval domain(t0, t1);
        // Contrainte sur l'etat
        ctc_state.set_restricted_domain(domain);
        ctc_state.contract(x);
        // Dynamique de positions
        ctc_dyn.set_restricted_domain(domain);
        ctc_dyn.contract(ddx[pos_e],  dddx[pos_e]);
        ctc_dyn.contract(ddx[pos_n],  dddx[pos_n]);
        ctc_dyn.contract(ddx[pos_u],  dddx[pos_u]);
        ctc_dyn.contract(dx,  ddx);
        ctc_dyn.contract(x,  dx);
        // Publication
        #ifdef DEBUG
        std::cout << "Pose : " << x(t) << "\n";
        IntervalVector vars_dang(12);
        vars_dang[3] = x[cos_psi].interpol(t, dx[cos_psi]);
        vars_dang[4] = x[sin_psi].interpol(t, dx[sin_psi]);
        vars_dang[5] = x[cos_tht].interpol(t, dx[cos_tht]);
        vars_dang[6] = x[sin_tht].interpol(t, dx[sin_tht]);
        vars_dang[7] = x[cos_phi].interpol(t, dx[cos_phi]);
        vars_dang[8] = x[sin_phi].interpol(t, dx[sin_phi]);
        if (!vars_dang[3].contains(0))
          vars_dang[0] = dx[sin_psi].interpol(t, ddx[sin_psi])/vars_dang[3];
        else
          vars_dang[0] = -dx[cos_psi].interpol(t, ddx[cos_psi])/vars_dang[4];
        if (!vars_dang[5].contains(0))
          vars_dang[1] = dx[sin_tht].interpol(t, ddx[sin_tht])/vars_dang[5];
        else
          vars_dang[1] = -dx[cos_tht].interpol(t, ddx[cos_tht])/vars_dang[6];
        if (!vars_dang[7].contains(0))
          vars_dang[2] = dx[sin_phi].interpol(t, ddx[sin_phi])/vars_dang[7];
        else
          vars_dang[2] = -dx[cos_phi].interpol(t, ddx[cos_phi])/vars_dang[8];
        IntervalVector dx_psithtphi = f_matTi.eval_vector(vars_dang);
        std::cout << "Twist : " << dx.subvector(pos_e, pos_u)(t) << dx_psithtphi << "\n";
        printf("Time taken DYN: %.2fs\n", (double)(clock() - t_start_dyn)/CLOCKS_PER_SEC);
        #endif
      }
    }

  /* =========== GRAPHICS =========== */

    vibes::beginDrawing();

    const char* vars_x[x.size()];
    vars_x[pos_e]   = "pos_e";
    vars_x[pos_n]   = "pos_n";
    vars_x[pos_u]   = "pos_u";
    vars_x[cos_psi] = "cos_psi";
    vars_x[sin_psi] = "sin_psi";
    vars_x[cos_tht] = "cos_tht";
    vars_x[sin_tht] = "sin_tht";
    vars_x[cos_phi] = "cos_phi";
    vars_x[sin_phi] = "sin_phi";

    for (i=0; i<x.size(); i++) {
      VIBesFigTube fig(vars_x[i], &x[i]);
      fig.set_properties(100, 100, 600, 300);
      fig.show();
    }

    const char* vars_dx[3];
    vars_x[pos_e] = "vit_e";
    vars_x[pos_n] = "vit_n";
    vars_x[pos_u] = "vit_u";
    for (i=0; i<=3; i++) {
      VIBesFigTube fig(vars_dx[i], &dx[i]);
      fig.set_properties(100, 100, 600, 300);
      fig.show();
    }

    vibes::endDrawing();

    printf("Time taken: %.2fs\n", (double)(clock() - t_start)/CLOCKS_PER_SEC);
    std::cout << "Slices = " << x.nb_slices() <<"\n";
}
