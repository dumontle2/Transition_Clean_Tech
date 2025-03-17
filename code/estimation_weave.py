import os
import sys

import numpy as np
import weave
import scipy.stats as stats

import load_data as ld
import infinite_weave as iw

#-----------------------------#
# constants                   #
#-----------------------------#

# simulation alg params
r_per = 100
t_tot = 5.0
prod_max = 200
nf_inc = 16384
nf_ent_max = 5*nf_inc
SEED = 578945342
pgen_frac = 1.0
rand_frac = 1.0

delt = 1.0/r_per
n_rep = np.int(np.round(r_per*t_tot))
nf_max = nf_inc + nf_ent_max

#-----------------------------#
# simulation setup            #
#-----------------------------#

# allocate random vectors
pgen_size = np.int(np.round(pgen_frac*nf_max*t_tot))
rand_size = np.int(np.round(rand_frac*nf_max*n_rep))

pgen_clean_vec = np.zeros(pgen_size)
pgen_dirty_vec = np.zeros(pgen_size)
rand_vec = np.zeros(rand_size)

ent_clean_nstep_vec = np.zeros(nf_ent_max)
ent_clean_nprod_vec = np.zeros(nf_ent_max,dtype=np.int)

ent_dirty_nstep_vec = np.zeros(nf_ent_max)
ent_dirty_nprod_vec = np.zeros(nf_ent_max,dtype=np.int)

# starting time (for entrants)
rep0_vec = np.zeros(nf_max,dtype=np.int)

# output arrays
nprod_out = np.zeros(nf_max,dtype=np.int)
nrnd_out = np.zeros(nf_max,dtype=np.int)
pattot_out = np.zeros(nf_max,dtype=np.int)
active_out = np.zeros(nf_max,dtype=np.int)
zeros_out = np.zeros(nf_max,dtype=np.int)
nstep_out = np.zeros((nf_max,prod_max),dtype=np.int)
pgain_out = np.zeros(nf_max,dtype=np.int)
plose_out = np.zeros(nf_max,dtype=np.int)
ptest_out = np.zeros(nf_max,dtype=np.int)

#-----------------------------#
# generate random data        #
#-----------------------------#

gen_alpha = None
gen_gap_dist = None
def gen_rand_prods(alpha,gap_dist):
  global gen_alpha,gen_gap_dist,nd_zero,nf_clean,nf_dirty,mu_0,mu_dp,mu_cn
  global type_vec,sel_inc,sel_clean_inc,sel_dirty_inc,nd_clean_dist,nd_dirty_dist,nd_clean_cdf,nd_dirt_cdf
  global nprod_vec,nprod_tot,nrnd_vec,nstep_vec,zeros_vec,nprod_bin_0,rand_vec,pgen_clean_vec,pgen_dirty_vec,ent_clean_nstep_vec,ent_clean_nprod_vec,ent_dirty_nstep_vec,ent_dirty_nprod_vec

  print 'Loading ew gap dist {}'.format(gap_dist)

  # load data if needed
  if ld.loaded_dist_type != gap_dist:
    ld.load_gap_dist(gap_dist)

  # store alpha used to generate this data
  gen_alpha = alpha
  gen_gap_dist = gap_dist

  # seed rng
  np.random.seed(SEED)

  # product type fractions
  nd_zero = ld.nd_zero
  clean_nd_frac = 0.5*ld.nd_dist[nd_zero] + np.sum(ld.nd_dist[:nd_zero])

  meann_clean = np.sum(ld.prod_vals*ld.prod_clean_pdf)
  meann_dirty = np.sum(ld.prod_vals*ld.prod_dirty_pdf)
  clean_firm_frac = clean_nd_frac*meann_dirty/(clean_nd_frac*meann_dirty+(1.0-clean_nd_frac)*meann_clean)

  nf_clean = np.int(np.round(nf_inc*clean_firm_frac))
  nf_dirty = nf_inc - nf_clean

  type_vec = np.zeros(nf_max,dtype=np.int) # dirty = 0, clean = 1
  type_vec[:nf_clean] = 1

  # selections and counts
  sel_inc = np.zeros(nf_max,dtype=np.bool)
  sel_inc[:nf_inc] = True
  sel_clean_inc = (type_vec == 1) & sel_inc
  sel_dirty_inc = (type_vec == 0) & sel_inc

  # step differentials distributions
  nd_clean_dist = np.r_[ld.nd_dist[:nd_zero],0.5*ld.nd_dist[nd_zero],np.zeros(nd_zero)]
  nd_clean_dist /= np.sum(nd_clean_dist)
  nd_dirty_dist = np.r_[np.zeros(nd_zero),0.5*ld.nd_dist[nd_zero],ld.nd_dist[nd_zero+1:]]
  nd_dirty_dist /= np.sum(nd_dirty_dist)

  nd_clean_cdf = np.cumsum(nd_clean_dist)
  nd_dirty_cdf = np.cumsum(nd_dirty_dist)

  # mu aggregates
  mu_0 = ld.nd_dist[nd_zero]
  mu_dp = np.sum(ld.nd_dist[nd_zero+1:])
  mu_cn = np.sum(ld.nd_dist[:nd_zero])

  # product counts
  nprod_vec = np.zeros(nf_max,dtype=np.int)
  nprod_vec[:nf_clean] = np.digitize(np.random.random(size=nf_clean),ld.prod_clean_cdf)+1
  nprod_vec[nf_clean:nf_inc] = np.digitize(np.random.random(size=nf_dirty),ld.prod_dirty_cdf)+1

  nprod_clean = np.sum(nprod_vec[:nf_clean])
  nprod_dirty = np.sum(nprod_vec[nf_clean:nf_inc])
  nprod_tot = nprod_clean + nprod_dirty

  # distribute these as poisson(nprod), law of large numbers says sum should be about the same
  prod_frac_clean = np.float(nprod_clean)/np.float(nprod_tot)
  prod_frac_dirty = np.float(nprod_dirty)/np.float(nprod_tot)
  prod_frac_type = np.array([prod_frac_dirty,prod_frac_clean])
  nrnd_vec = np.random.poisson(lam=nprod_vec/prod_frac_type[type_vec])

  # draw step differentials for incumbents
  all_clean_prods = np.digitize(np.random.random(size=nprod_clean),nd_clean_cdf)-nd_zero
  all_dirty_prods = np.digitize(np.random.random(size=nprod_dirty),nd_dirty_cdf)-nd_zero

  clean_prod_idx = np.r_[0,np.cumsum(nprod_vec[sel_clean_inc])]
  dirty_prod_idx = np.r_[0,np.cumsum(nprod_vec[sel_dirty_inc])]

  nstep_vec = np.zeros((nf_max,prod_max),dtype=np.int)
  nstep_vec[sel_clean_inc,:] = np.vstack([np.r_[all_clean_prods[i1:i2],np.zeros(prod_max-p)] for (i1,i2,p) in zip(clean_prod_idx[:-1],clean_prod_idx[1:],nprod_vec[sel_clean_inc])])
  nstep_vec[sel_dirty_inc,:] = np.vstack([np.r_[all_dirty_prods[i1:i2],np.zeros(prod_max-p)] for (i1,i2,p) in zip(dirty_prod_idx[:-1],dirty_prod_idx[1:],nprod_vec[sel_dirty_inc])])

  # find firm employment
  zeros_vec = np.array([np.sum(steps[:n]==0) for (n,steps) in zip(nprod_vec,nstep_vec)])

  # find initial size bins
  nprod_median_0 = np.median(nprod_vec[:nf_inc])
  nprod_bin_0 = np.zeros(nf_max,dtype=np.int)
  nprod_bin_0[:nf_inc] = (nprod_vec[:nf_inc]>nprod_median_0)+1 # out = 0, small = 1, large = 2

  # find distribution of innovated on step sizes (before innovation, so leap-frogging draws n=0)
  inno_clean_pdf = ld.nd_dist.copy()
  inno_clean_pdf[nd_zero+1:] *= 1.0-alpha
  inno_clean_pdf[nd_zero] += alpha*mu_dp
  inno_clean_cdf = np.cumsum(inno_clean_pdf)

  inno_dirty_pdf = ld.nd_dist.copy()
  inno_dirty_pdf[:nd_zero] *= 1.0-alpha
  inno_dirty_pdf[nd_zero] += alpha*mu_cn
  inno_dirty_cdf = np.cumsum(inno_dirty_pdf)

  # generate random numbers
  pgen_clean_vec = np.digitize(np.random.random(size=pgen_size),inno_clean_cdf)-nd_zero-1
  pgen_dirty_vec = np.digitize(np.random.random(size=pgen_size),inno_dirty_cdf)-nd_zero+1
  rand_vec = np.random.random(size=rand_size)

  ent_clean_nstep_vec = np.digitize(np.random.random(size=nf_ent_max),inno_clean_cdf)-nd_zero-1
  ent_clean_nprod_vec = np.zeros(nf_ent_max,dtype=np.int)
  ent_clean_nprod_vec[ent_clean_nstep_vec<0] = 1
  ent_clean_nprod_vec[(ent_clean_nstep_vec==0)*(np.random.random(size=nf_ent_max)<0.5)] = 1

  ent_dirty_nstep_vec = np.digitize(np.random.random(size=nf_ent_max),inno_dirty_cdf)-nd_zero+1
  ent_dirty_nprod_vec = np.zeros(nf_ent_max,dtype=np.int)
  ent_dirty_nprod_vec[ent_dirty_nstep_vec>0] = 1
  ent_dirty_nprod_vec[(ent_dirty_nstep_vec==0)*(np.random.random(size=nf_ent_max)<0.5)] = 1

#-----------------------------#
# simulation functions        #
#-----------------------------#

def print_vec(v,prec=5):
  return ', '.join(len(v)*['{}:0<9.{}{}'.format('{',prec,'}')]).format(*v)

def smm_obj(pvec,output=False,best_val=-np.inf,disc=iw.disc0,alpha=iw.alpha0,eta=iw.eta0,nu=iw.nu0,gap_dist=iw.default_gap_dist,emit_model='output',ce_growth=iw.ce_growth0,ce_min=iw.ce_min0):
  if alpha != gen_alpha or gap_dist != gen_gap_dist:
    print 'REGENERATING PRODUCTS FOR ALPHA {} -> {} AND GAP DIST {} -> {}'.format(gen_alpha,alpha,gen_gap_dist,gap_dist)
    gen_rand_prods(alpha,gap_dist)

  if np.any(pvec<=0.0):
    print '({}): {:s}'.format(print_vec(pvec),'Negative param')
    if not output: return np.inf

  (lam0,theta,F_I,F_E,Rmax,ce_0) = pvec
  lam = 1.0+lam0

  if output:
    print 'lam0  = {:.5f}'.format(lam0)
    print 'theta = {:.5f}'.format(theta)
    print 'F_E   = {:.5f}'.format(F_E)
    print 'F_I   = {:.5f}'.format(F_I)
    print 'Rmax  = {:7.0f}'.format(Rmax)
    print 'ce    = {:.5f}'.format(ce_0)

  # zero out policy
  iw.m_vec[:,0] = 0
  iw.subs_vec[:,0] = 0.0

  # call equilibrium solver
  iw.sim_pvec(pvec=pvec,output=output,disc=disc,alpha=alpha,eta=eta,nu=nu,gap_dist=gap_dist,emit_model=emit_model,ce_growth=ce_growth,ce_min=ce_min)

  # use time zero equilibrium values
  x_c = iw.xc_vec[0,0]
  x_d = iw.xd_vec[0,0]
  tau_c = iw.tauc_vec[0,0]
  tau_d = iw.taud_vec[0,0]
  tau_ed = iw.taued_vec[0,0]
  tau_ec = iw.tauec_vec[0,0]
  ws = iw.ws_vec[0,0]

  tau_e = tau_ed + tau_ec

  if tau_e == 0.0:
    print 'No entry: bailing'
    print

  labor_c = np.power(x_c/theta,1.0/eta)+F_I
  labor_d = np.power(x_d/theta,1.0/eta)+F_I

  wu = mu_0 + (1.0/lam)*(1.0-mu_0)

  Gamma_d = (mu_0+mu_dp) + alpha*mu_cn
  Gamma_c = (mu_0+mu_cn) + alpha*mu_dp
  growth = np.log(lam)*(tau_c*Gamma_c+tau_d*Gamma_d)

  if output:
    print 'Eqvars:'
    print '  growth = {}'.format(growth)
    print '  ws = {}'.format(ws)
    print '  wu = {}'.format(wu)
    print '  tau = {}'.format(tau_c+tau_d)
    print '  tau_c = {}'.format(tau_c)
    print '  tau_d = {}'.format(tau_d)
    print '  tau_e = {}'.format(tau_e)
    print '  x_c = {}'.format(x_c)
    print '  x_d = {}'.format(x_d)
    print '  tau_ec = {}'.format(tau_ec)
    print '  tau_ed = {}'.format(tau_ed)
    print '  Gamma_d = {}'.format(Gamma_d)
    print '  Gamma_c = {}'.format(Gamma_c)
    print

  #-----------------------------#
  # calculate emissions targets #
  #-----------------------------#

  emit_delt_back = 30
  emit_delt_fwrd = 10
  emit_delt_bin = np.int(np.floor(emit_delt_fwrd/iw.delt))
  emit_world_zero = iw.Emit_follow_vec[0,0]
  emit_world_data = ld.wem_emit_cut[-1]
  emit_change_fwrd  = (iw.Emit_follow_vec[emit_delt_bin,0]/iw.Emit_follow_vec[0,0]-1.0)/emit_delt_fwrd
  emit_change_back = (ld.wem_emit_cut[-1]/ld.wem_emit_cut[-(emit_delt_back+1)]-1.0)/emit_delt_back

  #-----------------------------#
  # simulate firm panel         #
  #-----------------------------#

  # how many entrants?
  nf_ent_d = np.int(np.round(tau_ed*t_tot*nprod_tot))
  nf_ent_c = np.int(np.round(tau_ec*t_tot*nprod_tot))
  nf_ent = nf_ent_d + nf_ent_c
  nf = nf_inc + nf_ent

  if output:
    print 'nf_ent = {}, nf_ent_max = {}'.format(nf_ent,nf_ent_max)
    print

  if nf > nf_max:
    print 'Too many entrants: nf_ent = {}, nf_ent_max = {}'.format(nf_ent,nf_ent_max)
    return np.inf

  # initial entrant states
  type_vec[nf_inc:nf_inc+nf_ent_d] = 0
  type_vec[nf_inc+nf_ent_d:nf] = 1
  nprod_vec[nf_inc:nf_inc+nf_ent_d] = ent_dirty_nprod_vec[:nf_ent_d]
  nprod_vec[nf_inc+nf_ent_d:nf] = ent_clean_nprod_vec[:nf_ent_c]
  nstep_vec[nf_inc:nf_inc+nf_ent_d,0] = ent_dirty_nstep_vec[:nf_ent_d]
  nstep_vec[nf_inc+nf_ent_d:nf,0] = ent_clean_nstep_vec[:nf_ent_c]
  rep0_vec[nf_inc:nf_inc+nf_ent_d] = np.round(np.linspace(0.0,1.0-1.0/n_rep,nf_ent_d)*n_rep)
  rep0_vec[nf_inc+nf_ent_d:nf] = np.round(np.linspace(0.0,1.0-1.0/n_rep,nf_ent_c)*n_rep)

  # run c++ simulation code with weave
  code = open('simulation.cpp','r').read()
  args_param = ['alpha','delt','n_rep','prod_max','nf','nd_zero','r_per']
  args_eqv = ['x_c','x_d','tau_c','tau_d']
  args_inp = ['nprod_vec','nrnd_vec','nstep_vec','type_vec','rep0_vec']
  args_out = ['nprod_out','nrnd_out','active_out','zeros_out','nstep_out','pgain_out','plose_out','ptest_out']
  args_rand = ['pgen_clean_vec','pgen_dirty_vec','rand_vec','pgen_size','rand_size']
  args = args_param + args_eqv + args_inp + args_out + args_rand
  weave.inline(code,arg_names=args)

  # restrict to our sample
  nprod_samp0 = nprod_vec[:nf]
  nrnd_samp0 = nrnd_vec[:nf]
  zeros_samp0 = zeros_vec[:nf]
  type_samp0 = type_vec[:nf]
  nprod_bin_samp0 = nprod_bin_0[:nf]
  nprod_samp1 = nprod_out[:nf]
  nrnd_samp1 = nrnd_out[:nf]
  active_samp1 = active_out[:nf]
  zeros_samp1 = zeros_out[:nf]
  type_samp1 = type_vec[:nf]

  # selections
  exited = active_samp1 == 0
  sel_posn_1 = nprod_samp1 > 0
  sel_active = ~exited
  sel_ent = (nprod_bin_samp0 == 0) & sel_active
  sel_inc = nprod_bin_samp0 > 0
  sel_growth = sel_inc & sel_active & sel_posn_1

  # counts
  nf_clean_1 = np.sum(active_samp1[:nf_clean])
  nf_dirty_1 = np.sum(active_samp1[nf_clean:nf_inc])
  nf_new_1 = np.sum(active_samp1[nf_inc:])
  nf_1 = nf_clean_1 + nf_dirty_1 + nf_new_1

  # find final size bins
  nprod_median_samp1 = np.median(nprod_samp1[sel_active])
  nprod_bin_samp1 = np.zeros(nf,dtype=np.int)
  nprod_bin_samp1[sel_active] = (nprod_samp1[sel_active]>nprod_median_samp1)+1

  # sales labor statistics
  rnd_labor = np.array([labor_d,labor_c])
  firm_rnd_labor_samp0 = rnd_labor[type_samp0]
  firm_rnd_labor_samp1 = rnd_labor[type_samp1]

  firm_empl_samp0 = (1.0/wu)*(zeros_samp0+(1.0/lam)*(nprod_samp0-zeros_samp0)) + firm_rnd_labor_samp0*(nprod_samp0+nrnd_samp0)
  firm_empl_samp1 = (1.0/wu)*(zeros_samp1+(1.0/lam)*(nprod_samp1-zeros_samp1)) + firm_rnd_labor_samp1*(nprod_samp1+nrnd_samp1)

  firm_sales_samp0 = nprod_samp0
  firm_sales_samp1 = nprod_samp1

  firm_sales_empl_samp0 = firm_sales_samp0/firm_empl_samp0
  firm_sales_empl_samp1 = firm_sales_samp1/firm_empl_samp1

  # entrant employment share
  sum_empl_inc = np.sum(firm_empl_samp1[sel_inc&sel_active])
  sum_empl_ent = np.sum(firm_empl_samp1[sel_ent&sel_active])
  ent_empl_frac = sum_empl_ent/sum_empl_inc
  ent_empl_frac_adj = ent_empl_frac/t_tot

  # R&D to sales ratio
  sel_rnd0 = sel_inc
  sel_rnd1 = sel_active * sel_posn_1
  firm_rnd_sales_0 = ws*firm_rnd_labor_samp0[sel_rnd0]*(nprod_samp0[sel_rnd0]+nrnd_samp0[sel_rnd0])/nprod_samp0[sel_rnd0]
  firm_rnd_sales_1 = ws*firm_rnd_labor_samp1[sel_rnd1]*(nprod_samp1[sel_rnd1]+nrnd_samp1[sel_rnd1])/nprod_samp1[sel_rnd1]
  rnd_sales_0 = np.mean(firm_rnd_sales_0)
  rnd_sales_1 = np.mean(firm_rnd_sales_1)
  mean_rnd_sales = 0.5*(rnd_sales_0+rnd_sales_1)

  # exit rate
  exit_rate = np.mean(exited[sel_inc])
  exit_rate_yr = exit_rate/t_tot

  if output:
    # entry ratios
    mean_empl_inc = np.mean(firm_empl_samp1[sel_inc&sel_active])
    mean_empl_ent = np.mean(firm_empl_samp1[sel_ent&sel_active])
    ent_mean_empl_frac = mean_empl_ent/mean_empl_inc

    mean_empl_inc = np.mean(firm_empl_samp1[sel_inc&sel_active])
    mean_empl_ent = np.mean(firm_empl_samp1[sel_ent&sel_active])
    ent_mean_empl_frac = mean_empl_ent/mean_empl_inc

    mean_sales_inc = np.mean(firm_sales_samp1[sel_inc&sel_active])
    mean_sales_ent = np.mean(firm_sales_samp1[sel_ent&sel_active])
    ent_mean_sales_frac = mean_sales_ent/mean_sales_inc

    firm_sales_empl_inc = firm_sales_samp1[sel_inc&sel_active]/firm_empl_samp1[sel_inc&sel_active]
    firm_sales_empl_ent = firm_sales_samp1[sel_ent&sel_active]/firm_empl_samp1[sel_ent&sel_active]

    mean_sales_empl_inc = np.mean(firm_sales_empl_inc)
    mean_sales_empl_ent = np.mean(firm_sales_empl_ent)
    ent_mean_sales_empl_frac = mean_sales_empl_ent/mean_sales_empl_inc

    # growth deciles
    n_quant = 5
    pct_quantiles = np.linspace(0.0,100.0,n_quant+1)[1:]

    empl_quantiles = np.array([stats.scoreatpercentile(firm_empl_samp0[sel_inc],pr) for pr in pct_quantiles])
    firm_empl_quantile = np.digitize(firm_empl_samp0,empl_quantiles)

    empl_growth_dhs_cond = (firm_empl_samp1[sel_growth]-firm_empl_samp0[sel_growth])/(0.5*(firm_empl_samp0[sel_growth]+firm_empl_samp1[sel_growth]))
    mean_quantile_empl_growth_cond = np.array([np.mean(empl_growth_dhs_cond[firm_empl_quantile[sel_growth]==d]) for d in range(n_quant)])

    # growth scaling
    growth_fracs_down = [0.25,0.50,0.75]
    growth_fracs_up = [1.25,1.50,1.75,2.00]
    nd_gf = len(growth_fracs_down)
    nu_gf = len(growth_fracs_up)

    empl_growth_unc_dprobs = [np.mean(firm_empl_samp1[sel_inc]<=gf*firm_empl_samp0[sel_inc]) for gf in growth_fracs_down]
    empl_growth_unc_uprobs = [np.mean(firm_empl_samp1[sel_inc]>=gf*firm_empl_samp0[sel_inc]) for gf in growth_fracs_up]

  # data moments targets
  (ent_empl_frac_dat,exit_rate_dat,rnd_sales_dat,growth_dat) = ld.mmt_data

  if nu == 0.0:
    new_mmt_sim = np.array([ent_empl_frac_adj,exit_rate_yr,mean_rnd_sales,growth])
    new_mmt_dat = np.array([ent_empl_frac_dat,exit_rate_dat,rnd_sales_dat,growth_dat])
    new_mmt_wgt = np.array([1.0,1.0,1.0,10.0])
    new_mmt_obj = np.mean(new_mmt_wgt*np.abs(new_mmt_sim-new_mmt_dat))
  else:
    new_mmt_sim = np.array([ent_empl_frac_adj,exit_rate_yr,mean_rnd_sales,growth,emit_world_zero,emit_change_fwrd])
    new_mmt_dat = np.array([ent_empl_frac_dat,exit_rate_dat,rnd_sales_dat,growth_dat,emit_world_data,emit_change_back])
    new_mmt_wgt = np.array([1.0,1.0,1.0,10.0,0.01,1.0])
    new_mmt_obj = np.mean(new_mmt_wgt*np.abs(new_mmt_sim-new_mmt_dat))

  if output or (new_mmt_obj<best_val):
    # main summary
    names = ['ent empl frac adj','exit rate','mean R&D/sales','agg sales/worker growth','world emissions','emissions change']
    print 'Estimation moments:'
    print '  {:30s}: {:>9s} {:>9s}'.format('name','simul','data')
    for i in range(len(new_mmt_dat)):
      print '  {:30s}: {:9.5f} {:9.5f}'.format(names[i],new_mmt_sim[i],new_mmt_dat[i])

  if output:
    print
    print 'Main summary:'
    print
    print 'Median entry ratios:'
    print '             {:>9s}  {:>9s}'.format('model','data')
    print 'employment = {:9.5f}  {:9.5f}'.format(ent_mean_empl_frac,0.03)
    print 'sales      = {:9.5f}  {:9.5f}'.format(ent_mean_sales_frac,0.20)
    print 'sales/empl = {:9.5f}  {:9.5f}'.format(ent_mean_sales_empl_frac,1.05)
    print
    print 'Unconditional employment growth bins:'
    growth_dprobs_data = [0.11,0.15,0.26]
    growth_uprobs_data = [0.31,0.20,0.14,0.11]
    print '{:>4s}:  {:>12s}  {:>12s}'.format('pct','model','data')
    for i in range(nd_gf):
      gfp = np.int(np.round(100.0*(growth_fracs_down[i]-1.0)))
      print '{:>3d}%:  {:12.5f}  {:12.5f}'.format(gfp,empl_growth_unc_dprobs[i],growth_dprobs_data[i])
    for i in range(nu_gf):
      gfp = np.int(np.round(100.0*(growth_fracs_up[i]-1.0)))
      print '{:>3d}%:  {:12.5f}  {:12.5f}'.format(gfp,empl_growth_unc_uprobs[i],growth_uprobs_data[i])
    print
    print 'Conditional employment growth deciles'
    empl_growth_data = [0.31,0.14,0.11,-0.01,-0.10]
    print ('{:>2s}:    {:>12s}  {:>12s}').format('q','model','data')
    for d in range(n_quant):
      print ('{:2d}:    {:>12.5f}  {:>12.5f}').format(d,mean_quantile_empl_growth_cond[d],empl_growth_data[d])
    print

  print '[{}]: {:.8f}'.format(print_vec(pvec),new_mmt_obj)
  return new_mmt_obj

def estimation_obj(**kwargs):
  return lambda pol,**ikwargs: smm_obj(pol,**dict(kwargs,**ikwargs))

# ghetto optimizer
def anneal0_est(f,x0,scale=0.25,maxiter=2500,**kwargs):
  n = len(x0)
  xp = x0
  fp = f(x0,best_val=np.inf,**kwargs)
  xmin = xp
  fmin = fp
  stay = 0
  for i in range(maxiter):
    xp = xmin*np.exp(np.random.normal(loc=0.0,scale=scale,size=n))
    fp = f(xp,best_val=fmin,**kwargs)
    if fp < fmin and ~np.isnan(fp):
      xmin = xp
      fmin = fp
      print 'MIN -> ({},{})'.format(xmin,fmin)
      stay = 0
    else:
      stay += 1

    if stay == 50:
      scale *= 0.5
      print 'SCALE -> {}'.format(scale)
      stay = 0

    if scale < 1e-5:
      break

  print iw.quick_print(xmin)+' # {:<8.5f}'.format(-fmin)

  return (xmin,fmin)
