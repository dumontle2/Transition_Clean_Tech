#line 2 "infinite_weave.py"

#define min(a,b) ((a < b) ? a : b)
#define max(a,b) ((a > b) ? a : b)
#define clamp(x,a,b) (x < a ? a : (x > b ? b : x))

class EqSolver
{
  public:
    EqSolver(double delt, int nn, int nt, int n_zero, double v_fact, double mu_fact, double* m_vec, double* subs_vec, double theta, double eta, double L, double F_I, double F_E, double G_I, double G_E, double alpha, double prof, double disc, double ev_ss, double tau_ss, double x_ss, double ws_ss, long int* n_vec, double* dmu_vec, double* flow_d_vec, double* flow_c_vec, double* vfunc_d_new, double* vfunc_c_new, double* rmg_vec, double* mu_vec, double* vfunc_d_vec, double* vfunc_c_vec, double* evd_vec, double* evc_vec, double* xd_vec, double* xc_vec, double* xpd_vec, double* xpc_vec, double* taud_vec, double* tauc_vec, double* ws_vec, double* wu_vec, double* taued_vec, double* tauec_vec, double* retd_vec, double* retc_vec, double* posd_vec, double* posc_vec, double* evcp_vec, double* xe_vec, double* labd_vec, double* labc_vec, double* labec_vec, double* labed_vec, double* alpha_c_vec, double* alpha_d_vec, double* prof_c_vec, double* prof_d_vec, double gamma, double nu, double lam, double lam0, double lbar, double* US_emit_ratio, double growth_other_emit, double* me_vec, double* ce_vec, double* pe_vec, double* R_used_vec, double* Emit_other_vec, double* slamd_vec, double* slamc_vec, double* R_usage_vec, double* dirty_labor_vec) {
this->delt = delt;
this->nn = nn;
this->nt = nt;
this->n_zero = n_zero;
this->v_fact = v_fact;
this->mu_fact = mu_fact;
this->m_vec = m_vec;
this->subs_vec = subs_vec;
this->theta = theta;
this->eta = eta;
this->L = L;
this->F_I = F_I;
this->F_E = F_E;
this->G_I = G_I;
this->G_E = G_E;
this->alpha = alpha;
this->prof = prof;
this->disc = disc;
this->ev_ss = ev_ss;
this->tau_ss = tau_ss;
this->x_ss = x_ss;
this->ws_ss = ws_ss;
this->n_vec = n_vec;
this->dmu_vec = dmu_vec;
this->flow_d_vec = flow_d_vec;
this->flow_c_vec = flow_c_vec;
this->vfunc_d_new = vfunc_d_new;
this->vfunc_c_new = vfunc_c_new;
this->rmg_vec = rmg_vec;
this->mu_vec = mu_vec;
this->vfunc_d_vec = vfunc_d_vec;
this->vfunc_c_vec = vfunc_c_vec;
this->evd_vec = evd_vec;
this->evc_vec = evc_vec;
this->xd_vec = xd_vec;
this->xc_vec = xc_vec;
this->xpd_vec = xpd_vec;
this->xpc_vec = xpc_vec;
this->taud_vec = taud_vec;
this->tauc_vec = tauc_vec;
this->ws_vec = ws_vec;
this->wu_vec = wu_vec;
this->taued_vec = taued_vec;
this->tauec_vec = tauec_vec;
this->retd_vec = retd_vec;
this->retc_vec = retc_vec;
this->posd_vec = posd_vec;
this->posc_vec = posc_vec;
this->evcp_vec = evcp_vec;
this->xe_vec = xe_vec;
this->labd_vec = labd_vec;
this->labc_vec = labc_vec;
this->labec_vec = labec_vec;
this->labed_vec = labed_vec;
this->alpha_c_vec = alpha_c_vec;
this->alpha_d_vec = alpha_d_vec;
this->prof_c_vec = prof_c_vec;
this->prof_d_vec = prof_d_vec;
this->gamma = gamma;
this->nu = nu;
this->lam = lam;
this->lam0 = lam0;
this->lbar = lbar;
this->US_emit_ratio = US_emit_ratio;
this->growth_other_emit = growth_other_emit;
this->me_vec = me_vec;
this->ce_vec = ce_vec;
this->pe_vec = pe_vec;
this->R_used_vec = R_used_vec;
this->Emit_other_vec = Emit_other_vec;
this->slamd_vec = slamd_vec;
this->slamc_vec = slamc_vec;
this->R_usage_vec = R_usage_vec;
this->dirty_labor_vec = dirty_labor_vec;
}

double delt;
int nn;
int nt;
int n_zero;
double v_fact;
double mu_fact;
double* m_vec;
double* subs_vec;
double theta;
double eta;
double L;
double F_I;
double F_E;
double G_I;
double G_E;
double alpha;
double prof;
double disc;
double ev_ss;
double tau_ss;
double x_ss;
double ws_ss;
long int* n_vec;
double* dmu_vec;
double* flow_d_vec;
double* flow_c_vec;
double* vfunc_d_new;
double* vfunc_c_new;
double* rmg_vec;
double* mu_vec;
double* vfunc_d_vec;
double* vfunc_c_vec;
double* evd_vec;
double* evc_vec;
double* xd_vec;
double* xc_vec;
double* xpd_vec;
double* xpc_vec;
double* taud_vec;
double* tauc_vec;
double* ws_vec;
double* wu_vec;
double* taued_vec;
double* tauec_vec;
double* retd_vec;
double* retc_vec;
double* posd_vec;
double* posc_vec;
double* evcp_vec;
double* xe_vec;
double* labd_vec;
double* labc_vec;
double* labec_vec;
double* labed_vec;
double* alpha_c_vec;
double* alpha_d_vec;
double* prof_c_vec;
double* prof_d_vec;
double gamma;
double nu;
double lam;
double lam0;
double lbar;
double* US_emit_ratio;
double growth_other_emit;
double* me_vec;
double* ce_vec;
double* pe_vec;
double* R_used_vec;
double* Emit_other_vec;
double* slamd_vec;
double* slamc_vec;
double* R_usage_vec;
double* dirty_labor_vec;


    // variables
    double mv,subst,rmg,ev_d,ev_c,ws,x_d,x_c,x_e,return_d,return_c,labor_d,labor_c,labor_e,M_e,M_ec,M_ed,M_d,M_c,tev_d,tev_c,tev,L_var,tau_d,tau_c,tau,rnd_cost_d,rnd_cost_c,vss_frac,pos_d,pos_c,xp_d,xp_c;
    int t,i,j;

    int crossover;
    double ent_clean_frac,ent_dirty_frac;
    double ws_lo,ws_hi,ws_mid,fw_lo,fw_hi,fw_mid,lret,FI_min,FI_max,F_slope,Fstar_d,Fstar_c;

    // emissions section
    int me,mv0,mv1,mt0,mt1;
    double mf,mu_0,mu_1,mu_d,mu_c,mu_dp,mu_cn,prof_frac_d,prof_frac_c,tax,itax,Omega,labor_ed,labor_ec,tau_ed,tau_ec,dirty_sum,Omega_d,ce,pe,R_used,R_usage,dirty_labor,R_usage_world;

    // state vectors
    double* mu_base;
    double* mu_next;

    double* vfunc_d_base;
    double* vfunc_c_base;

    double* vfunc_d_next;
    double* vfunc_c_next;

    double* vfunc_d_save;
    double* vfunc_c_save;

    double* vfunc_c_prime;
    double* vfunc_d_prime;

    void init_begin()
    {
      mu_base = mu_vec;
      mu_next = mu_vec + nn;

      vfunc_d_base = vfunc_d_vec;
      vfunc_c_base = vfunc_c_vec;

      vfunc_d_next = vfunc_d_vec + nn;
      vfunc_c_next = vfunc_c_vec + nn;

      vfunc_d_save = vfunc_d_new;
      vfunc_c_save = vfunc_c_new;
    }

    void increment()
    {
      mu_base += nn;
      mu_next += nn;

      vfunc_d_base += nn;
      vfunc_c_base += nn;

      vfunc_d_next += nn;
      vfunc_c_next += nn;

      vfunc_d_save += nn;
      vfunc_c_save += nn;
    }

    void decrement()
    {
      mu_base -= nn;
      mu_next -= nn;

      vfunc_d_base -= nn;
      vfunc_c_base -= nn;

      vfunc_d_next -= nn;
      vfunc_c_next -= nn;

      vfunc_d_save -= nn;
      vfunc_c_save -= nn;
    }

    void get_policy()
    {
      // policy
      mv = m_vec[t];
      me = me_vec[t];
      subst = subs_vec[t];

      mv0 = floor(me);
      mv1 = mv0 + 1;
      mf = 1.0 - (me-mv0);

      mt0 = mv0 + n_zero;
      mt1 = mv1 + n_zero;

      for (i = 0; i < nn; i++) prof_c_vec[i] = prof*clamp(me-n_vec[i],0.0,1.0);
      for (i = 0; i < nn; i++) prof_d_vec[i] = prof*clamp(n_vec[i]-me,0.0,1.0);

      // interest rate
      rmg = rmg_vec[t];
    }

    void get_ev()
    {
      // expected vp
      ev_d = 0.0;
      for (i = 0; i < n_zero; i++) ev_d += mu_base[i]*(alpha*vfunc_d_next[n_zero+1]+(1.0-alpha)*vfunc_d_next[i+1]);
      for (i = n_zero; i < nn-1; i++) ev_d += mu_base[i]*vfunc_d_next[i+1];
      ev_d += mu_base[nn-1]*vfunc_d_next[nn-1];

      ev_c = 0.0;
      ev_c += mu_base[0]*vfunc_c_next[0];
      for (i = 1; i <= n_zero; i++) ev_c += mu_base[i]*vfunc_c_next[i-1];
      for (i = n_zero+1; i < nn; i++) ev_c += mu_base[i]*(alpha*vfunc_c_next[n_zero-1]+(1.0-alpha)*vfunc_c_next[i-1]);

      ev_d = max(0.0,ev_d);
      ev_c = max(0.0,ev_c);
    }

    double labor_obj(double wsp)
    {
      // innovation intensities
      if (wsp == 0.0) {
        x_d = 0.0;
        x_c = 0.0;
      } else {
        x_d = pow(pow(theta,1.0/eta)*eta*ev_d/wsp,eta/(1.0-eta));
        x_c = pow(pow(theta,1.0/eta)*eta*ev_c/(wsp*(1.0-subst)),eta/(1.0-eta));
      }

      // return to innovation
      return_d = x_d*ev_d-wsp*pow(x_d/theta,1.0/eta);
      return_c = x_c*ev_c-wsp*(1.0-subst)*pow(x_c/theta,1.0/eta);

      // fixed cost cut out
      if (FI_max == 0.0) {
        Fstar_d = 0.0;
        Fstar_c = 0.0;
        pos_d = 1.0;
        pos_c = 1.0;
      } else {
        Fstar_d = max(FI_min,min(FI_max,return_d/wsp));
        Fstar_c = max(FI_min,min(FI_max,return_c/((1.0-subst)*wsp)));
        pos_d = F_slope*(Fstar_d-FI_min);
        pos_c = F_slope*(Fstar_c-FI_min);
      }

      // skilled labor usage
      labor_d = pos_d*(pow(x_d/theta,1.0/eta)+0.5*(FI_min+Fstar_d));
      labor_c = pos_c*(pow(x_c/theta,1.0/eta)+0.5*(FI_min+Fstar_c));

      return (labor_d+labor_c-L);
    }

    void solve_eq()
    {
      // skilled wage from entry condition
      ws = max(ev_d,ev_c/(1.0-subst))*theta*pow(eta,eta)*pow((1.0-eta)/F_E,1.0-eta);

      double Q = ((1.0-eta)/eta)*pow(theta*eta/ws,1.0/(1.0-eta))*(pow(ev_c/(1.0-subst),1.0/(1.0-eta))-pow(ev_d,1.0/(1.0-eta)))/F_E;
      if (Q >= G_E) {
        ent_clean_frac = 1.0;
        ent_dirty_frac = 0.0;
      } else if (Q <= -G_E) {
        ent_clean_frac = 0.0;
        ent_dirty_frac = 1.0;
      } else {
        //printf("t = %i: crossover\n",t);
        crossover = 1;
        ent_clean_frac = 0.5*(Q/G_E+1.0);
        ent_dirty_frac = 0.5*(1.0-Q/G_E);
      }

      // innovation rates and returns
      lret = labor_obj(ws);

      // zero entry, solve for wage
      if (lret > 0.0) {
        crossover = 0;

        ws_lo = 0.001;
        ws_hi = 50.0;
        fw_lo = labor_obj(ws_lo);
        fw_hi = labor_obj(ws_hi);
        j = 0;
        do {
          ws_mid = 0.5*(ws_lo+ws_hi);
          fw_mid = labor_obj(ws_mid);
          if (fw_mid > 0.0) {
            ws_lo = ws_mid;
          } else {
            ws_hi = ws_mid;
          }
          if (++j == 128) break;
        } while (fabs(fw_mid)>1e-6);

        ws = ws_mid;
        x_e = 0.0;
        labor_e = 0.0;
        M_e = 0.0;
        M_ed = 0.0;
        M_ec = 0.0;
      } else {
        labor_e = L-labor_d-labor_c;
        M_e = labor_e/(ent_clean_frac*pow(x_c/theta,1.0/eta)+ent_dirty_frac*pow(x_d/theta,1.0/eta)+F_E);
        M_ed = ent_dirty_frac*M_e;
        M_ec = ent_clean_frac*M_e;
      }

      // entry
      labor_ec = M_ec*labor_c;
      labor_ed = M_ed*labor_d;

      tau_ed = x_d*M_ed;
      tau_ec = x_c*M_ec;

      // aggregate rates
      M_d = pos_d + M_ed;
      M_c = pos_c + M_ec;

      xp_d = pos_d*x_d;
      xp_c = pos_c*x_c;

      tau_d = M_d*x_d;
      tau_c = M_c*x_c;
      tau = tau_d + tau_c;
    }

    void solve_energy() {
      // mu aggregates
      mu_0 = mu_base[mt0]; // fractional profits for clean leader
      mu_1 = mu_base[mt1]; // fractional profits for dirty leader
      mu_dp = 0.0; for (i = mt1+1; i < nn; i++) mu_dp += mu_base[i];
      mu_cn = 0.0; for (i = 0; i < mt0; i++) mu_cn += mu_base[i];

      // taxes
      tax = pow(lam,mv)-1.0;
      itax = 1.0/(1.0+tax);

      prof_frac_c = pow(lam,mf); // at mv0
      prof_frac_d = pow(lam,1.0-mf); // at mv1

      // resource usage
      ce = ce_vec[t];
      pe = pe_vec[t];
      Omega = itax*((1.0-nu)+nu*(ce/pe))*((1.0/prof_frac_d)*mu_1+(1.0/lam)*mu_dp) + ((1.0/prof_frac_c)*mu_0+(1.0/lam)*mu_cn);
      Omega_d = (1.0-nu)*((1.0/prof_frac_d)*mu_1+(1.0/lam)*mu_dp)*itax;
      dirty_labor = Omega_d/Omega;
      R_usage = (nu/(1.0-nu))*(dirty_labor/pe);
      R_usage_world = R_usage/US_emit_ratio[t];

      // update resources
      R_used = R_used_vec[t];
      R_used_vec[t+1] = R_used + delt*R_usage_world;

      // store
      dirty_labor_vec[t] = dirty_labor;
      R_usage_vec[t] = R_usage;
    }

    void update_values()
    {
      // instantaneous flows
      rnd_cost_d = ws*labor_d;
      rnd_cost_c = ws*(1.0-subst)*labor_c;
      for (i = 0; i < nn; i++) flow_d_vec[i] = prof_d_vec[i] - rnd_cost_d;
      for (i = 0; i < nn; i++) flow_c_vec[i] = prof_c_vec[i] - rnd_cost_c;

      // update values
      vfunc_d_save[0] = delt*flow_d_vec[0] + (1.0-delt*rmg)*((1.0-delt*tau)*vfunc_d_next[0]+delt*xp_d*ev_d+delt*tau_c*vfunc_d_next[0]);
      for (i = 1; i <= n_zero; i++) vfunc_d_save[i] = delt*flow_d_vec[i] + (1.0-delt*rmg)*((1.0-delt*tau)*vfunc_d_next[i]+delt*xp_d*ev_d+delt*tau_c*vfunc_d_next[i-1]);
      for (i = n_zero+1; i < nn; i++) vfunc_d_save[i] = delt*flow_d_vec[i] + (1.0-delt*rmg)*((1.0-delt*tau)*vfunc_d_next[i]+delt*xp_d*ev_d+delt*tau_c*(1.0-alpha)*vfunc_d_next[i-1]+delt*tau_c*alpha*vfunc_d_next[n_zero-1]);
      for (i = 0; i < nn; i++) vfunc_d_base[i] = (1.0-v_fact)*vfunc_d_base[i] + v_fact*vfunc_d_save[i];

      for (i = 0; i < n_zero; i++) vfunc_c_save[i] = delt*flow_c_vec[i] + (1.0-delt*rmg)*((1.0-delt*tau)*vfunc_c_next[i]+delt*xp_c*ev_c+delt*tau_d*(1.0-alpha)*vfunc_c_next[i+1]+delt*tau_d*alpha*vfunc_c_next[n_zero+1]);
      for (i = n_zero; i < nn-1; i++) vfunc_c_save[i] = delt*flow_c_vec[i] + (1.0-delt*rmg)*((1.0-delt*tau)*vfunc_c_next[i]+delt*xp_c*ev_c+delt*tau_d*vfunc_c_next[i+1]);
      vfunc_c_save[nn-1] = delt*flow_c_vec[nn-1] + (1.0-delt*rmg)*((1.0-delt*tau)*vfunc_c_next[nn-1]+delt*xp_c*ev_c+delt*tau_d*vfunc_c_next[nn-1]);
      for (i = 0; i < nn; i++) vfunc_c_base[i] = (1.0-v_fact)*vfunc_c_base[i] + v_fact*vfunc_c_save[i];
    }

    // iterate mu distribution
    void update_dists()
    {
      dmu_vec[0] = tau_c*mu_base[1] - tau_d*mu_base[0];
      for (i = 1; i < n_zero; i++) dmu_vec[i] = (1.0-alpha)*tau_d*mu_base[i-1] + tau_c*mu_base[i+1] - tau*mu_base[i];
      dmu_vec[n_zero] = (1.0-alpha)*tau_d*mu_base[n_zero-1] + (1.0-alpha)*tau_c*mu_base[n_zero+1] - tau*mu_base[n_zero];
      for (i = n_zero+1; i < nn-1; i++) dmu_vec[i] = tau_d*mu_base[i-1] + (1.0-alpha)*tau_c*mu_base[i+1] - tau*mu_base[i];
      dmu_vec[nn-1] = tau_d*mu_base[nn-2] - tau_c*mu_base[nn-1];
      for (i = 0; i < n_zero; i++) dmu_vec[n_zero+1] += alpha*tau_d*mu_base[i];
      for (i = n_zero+1; i < nn; i++) dmu_vec[n_zero-1] += alpha*tau_c*mu_base[i];
      for (i = 0; i < nn; i++) mu_next[i] = (1.0-mu_fact)*mu_next[i] + mu_fact*(mu_base[i] + delt*dmu_vec[i]);
    }

    // boundary condition (steady state) - because the dirty resource is exhaustible, this will always be clean
    void steady_state()
    {
      // policy
      mv = m_vec[nt-1];
      subst = subs_vec[nt-1];

      for (i = 0; i < nn; i++) prof_c_vec[i] = prof*clamp(me-n_vec[i],0.0,1.0);
      for (i = 0; i < nn; i++) prof_d_vec[i] = prof*clamp(n_vec[i]-me,0.0,1.0);

      for (i = 0; i < nn; i++) vfunc_d_save[i] = 0.0;
      for (i = 0; i < n_zero; i++) vfunc_c_save[i] = ev_ss;
      for (i = n_zero; i < nn; i++) vfunc_c_save[i] = double(ev_ss) - prof/(disc+double(tau_ss));

      x_d = 0.0;
      x_c = x_ss;
      tau_d = 0.0;
      tau_c = tau_ss;
      ev_d = 0.0;
      ev_c = ev_ss;
      ws = ws_ss;

      for (i = 0; i < nn; i++) vfunc_d_base[i] = (1.0-v_fact)*vfunc_d_base[i] + v_fact*vfunc_d_save[i];
      for (i = 0; i < nn; i++) vfunc_c_base[i] = (1.0-v_fact)*vfunc_c_base[i] + v_fact*vfunc_c_save[i];
    }

    // store innovation rates
    void store_vars()
    {
      xd_vec[t] = x_d;
      xc_vec[t] = x_c;
      xe_vec[t] = x_e;
      xpd_vec[t] = xp_d;
      xpc_vec[t] = xp_c;
      taud_vec[t] = tau_d;
      tauc_vec[t] = tau_c;
      ws_vec[t] = ws;
      wu_vec[t] = Omega;
      evd_vec[t] = ev_d;
      evc_vec[t] = ev_c;
      evcp_vec[t] = ev_c/(1.0-subst);
      taued_vec[t] = tau_ed;
      tauec_vec[t] = tau_ec;
      posd_vec[t] = pos_d;
      posc_vec[t] = pos_c;
      retd_vec[t] = return_d;
      retc_vec[t] = return_c;
      labd_vec[t] = labor_d;
      labc_vec[t] = labor_c;
      labec_vec[t] = labor_ec;
      labed_vec[t] = labor_ed;
    }

    void run_path()
    {
      crossover = 0;
      init_begin();

      for (t = 0; t < nt-1; t++) {
        get_policy();
        get_ev();
        solve_eq();
        solve_energy();
        update_dists();
        store_vars();
        increment();
      }

      steady_state();
      store_vars();

      for (t = nt-2; t >= 0; t--) {
        decrement();
        get_policy();
        get_ev();
        solve_eq();
        store_vars();
        update_values();
      }
    }

    void run()
    {
      FI_min = (1.0-G_I)*F_I;
      FI_max = (1.0+G_I)*F_I;
      F_slope = 1.0/(2.0*G_I*F_I);

      run_path();

      if (crossover == 1) printf("crossover\n");
    }
};

EqSolver eqs = EqSolver(delt,nn,nt,n_zero,v_fact,mu_fact,m_vec,subs_vec,theta,eta,L,F_I,F_E,G_I,G_E,alpha,prof,disc,ev_ss,tau_ss,x_ss,ws_ss,n_vec,dmu_vec,flow_d_vec,flow_c_vec,vfunc_d_new,vfunc_c_new,rmg_vec,mu_vec,vfunc_d_vec,vfunc_c_vec,evd_vec,evc_vec,xd_vec,xc_vec,xpd_vec,xpc_vec,taud_vec,tauc_vec,ws_vec,wu_vec,taued_vec,tauec_vec,retd_vec,retc_vec,posd_vec,posc_vec,evcp_vec,xe_vec,labd_vec,labc_vec,labec_vec,labed_vec,alpha_c_vec,alpha_d_vec,prof_c_vec,prof_d_vec,gamma,nu,lam,lam0,lbar,US_emit_ratio,growth_other_emit,me_vec,ce_vec,pe_vec,R_used_vec,Emit_other_vec,slamd_vec,slamc_vec,R_usage_vec,dirty_labor_vec);
eqs.run();
