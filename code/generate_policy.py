from os.path import join
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import infinite_weave as iw

plt.ioff()
sns.set(style='whitegrid',rc={'axes.titlesize': 18,
                              'axes.labelsize': 18,
                              'xtick.labelsize': 18,
                              'ytick.labelsize': 18,
                              'legend.fontsize': 14,
                              'figure.autolayout': True,
                              'lines.linewidth': 2.5})

colors = sns.color_palette()
(blue,green) = colors[:2]

img_format = 'png'
figsize0 = (5.5,4.5)
ymax0_s = 1.0
ymax0_t = 1.0

# options
base_dir = 'figures'
policies = iw.policies_opt
ptype0 = 'quartic'
soc_disc0 = iw.disc0
argz = {}

tcut = 200.0
tcut_bin = int(np.ceil((tcut/iw.t_max)*iw.nt))
t_vec = iw.t_vec[:tcut_bin,0]

n_ticks = 6

def generate_policy(name,pol_tag=None,par_tag=None,ptype=ptype0,label=None,ymax_s=ymax0_s,ymax_t=ymax0_t):
  if pol_tag is None: pol_tag = name
  if par_tag is None: par_tag = pol_tag
  if label is None: label = 'Optimal Policies'
  lam = iw.pvecs[par_tag][0] + 1.0
  figsize = (figsize0[0]+0.6,figsize0[1])

  # policy graphs
  (fig,ax_s) = plt.subplots(figsize=figsize)
  iw.calculate_policy(policies['{} {}'.format(pol_tag,ptype)],pol_type=ptype,lam=lam)
  (line_s,) = ax_s.plot(t_vec,iw.subs_vec[:tcut_bin,0],color=green,linestyle='dashed',label='$s$')
  ax_s.set_ylim(0.0,ymax_s)
  ax_s.set_yticks(np.linspace(0.0,ymax_s,n_ticks))
  ax_s.set_ylabel('Subsidy Rate')

  ax_t = ax_s.twinx()
  (line_t,) = ax_t.plot(t_vec,iw.tax_vec[:tcut_bin,0],color=blue,linestyle='solid',label='$\\tau$')
  ax_t.set_yticks(np.linspace(0.0,ymax_t,n_ticks))
  ax_t.set_ylim(0.0,ymax_t)
  ax_t.set_ylabel('Tax Rate')

  ax_s.set_xlabel('Number of Years')
  ax_s.set_title(label,y=1.03)

  lns = [line_s,line_t]
  labs = [l.get_label() for l in lns]
  ax_s.legend(lns,labs,loc='best')

  plt.savefig(join(base_dir,'{}_{}.{}'.format(name,ptype,img_format)))
  plt.close()

#--------------#
# output time! #
#--------------#

if __name__ == "__main__":
  # innovation rates and temperature - no policy
  iw.simulate_type(iw.policies['null'],**argz)

  (fig,ax) = plt.subplots(figsize=figsize0)
  ax.plot(t_vec,iw.taud_vec[:tcut_bin,0],'--',label='Dirty Innovation')
  ax.plot(t_vec,iw.tauc_vec[:tcut_bin,0],'-',label='Clean Innovation')
  ax.set_xlabel('Number of Years')
  ax.set_title('Innovation Rates, Null Policy',y=1.03)
  ax.legend(loc='best')
  plt.savefig(join(base_dir,'null_innovation_rates.{}'.format(img_format)))
  plt.close()

  (fig,ax) = plt.subplots(figsize=figsize0)
  ax.plot(t_vec,iw.Temp_follow_vec[:tcut_bin,0])
  ax.set_xlabel('Number of Years')
  ax.set_title('Temperature, Null Policy')
  plt.savefig(join(base_dir,'null_temperature.{}'.format(img_format)))
  plt.close()

  # initial productivity gap distribution
  (fig,ax) = plt.subplots(figsize=figsize0)
  ax.plot(iw.ld.nd_vals,iw.ld.nd_dist)
  ax.set_xlabel('Initial Productivity Gap')
  ax.set_title('Distribution of Initial Productivity Gaps',y=1.03)
  plt.savefig(join(base_dir,'initial_distribution.{}'.format(img_format)))
  plt.close()

  # innovation rates and temperature - optimal policies
  (_,_,welf_base) = iw.simulate_type(policies['baseline {}'.format(ptype0)],pol_type=ptype0,**argz)

  (fig,ax) = plt.subplots(figsize=figsize0)
  ax.plot(t_vec,iw.taud_vec[:tcut_bin,0],'--',label='Dirty')
  ax.plot(t_vec,iw.tauc_vec[:tcut_bin,0],'-',label='Clean')
  ax.set_xlabel('Number of Years')
  ax.set_title('Innovation Rates, Optimal Policies',y=1.03)
  ax.legend(loc='best')
  plt.savefig(join(base_dir,'optimal_innovation_rates.{}'.format(img_format)))
  plt.close()

  (fig,ax) = plt.subplots(figsize=figsize0)
  ax.plot(t_vec,iw.Temp_follow_vec[:tcut_bin,0])
  ax.set_xlabel('Number of Years')
  ax.set_title('Temperature, Optimal Policies',y=1.03)
  ax.set_ylim(1.0,2.0)
  plt.savefig(join(base_dir,'optimal_temperature.{}'.format(img_format)))
  plt.close()

  # innovation rates and temperature - constant policies
  (_,_,welf_const) = iw.simulate_type(policies['baseline const'],pol_type='const',**argz)

  (fig,ax) = plt.subplots(figsize=figsize0)
  ax.plot(t_vec,iw.taud_vec[:tcut_bin,0],'--',label='Dirty')
  ax.plot(t_vec,iw.tauc_vec[:tcut_bin,0],'-',label='Clean')
  ax.set_xlabel('Number of Years')
  ax.set_title('Innovation Rates, Constant Policies',y=1.03)
  ax.legend(loc='best')
  plt.savefig(join(base_dir,'const_innovation_rates.{}'.format(img_format)))
  plt.close()

  (fig,ax) = plt.subplots(figsize=figsize0)
  ax.plot(t_vec,iw.Temp_follow_vec[:tcut_bin,0])
  ax.set_xlabel('Number of Years')
  ax.set_title('Temperature, Constant Policies',y=1.03)
  ax.set_ylim(1.0,2.0)
  plt.savefig(join(base_dir,'const_temperature.{}'.format(img_format)))
  plt.close()

  # baseline optimal policy
  generate_policy('baseline')

  # carbon only policy
  (_,_,welf_carbon) = iw.simulate_type(policies['carbon {}'.format(ptype0)],pol_type=ptype0+'_carbon',**argz)

  (fig,ax) = plt.subplots(figsize=figsize0)
  ax.plot(t_vec,iw.tax_vec[:tcut_bin,0],label='$\\tau$')
  ax.set_ylim(0.0,3.5)
  ax.set_ylabel('Tax Rate')
  ax.set_xlabel('Number of Years')
  ax.set_title('Optimal Policy, Tax Only',y=1.03)
  ax.legend(loc='best')
  plt.savefig(join(base_dir,'counterfactual_carbon.{}'.format(img_format)))
  plt.close()

  # 50-year delay policies
  (_,_,welf_delay) = iw.simulate_type(policies['delay {}'.format(ptype0)],pol_type=ptype0+'_delay',**argz)

  figsize = (figsize0[0]+0.6,figsize0[1])
  (fig,ax_s) = plt.subplots(figsize=figsize)
  (line_s,) = ax_s.plot(t_vec,iw.subs_vec[:tcut_bin,0],color=green,linestyle='dashed',label='$s$')
  ax_s.set_ylim(0.0,ymax0_s)
  ax_s.set_yticks(np.linspace(0.0,ymax0_s,n_ticks))
  ax_s.set_ylabel('Subsidy Rate')

  ax_t = ax_s.twinx()
  (line_t,) = ax_t.plot(t_vec,iw.tax_vec[:tcut_bin,0],color=blue,linestyle='solid',label='$\\tau$')
  ax_t.set_yticks(np.linspace(0.0,ymax0_t,n_ticks))
  ax_t.set_ylim(0.0,ymax0_t)
  ax_t.set_ylabel('Tax Rate')

  ax_s.set_xlabel('Number of Years')
  ax_s.set_title('Optimal Policies, 50 Year Delay',y=1.03)

  lns = [line_s,line_t]
  labs = [l.get_label() for l in lns]
  ax_s.legend(lns,labs,loc='best')

  plt.savefig(join(base_dir,'counterfactual_delay.{}'.format(img_format)))
  plt.close()

  # business as usual
  (_,_,welf_bau) = iw.simulate_type(iw.policies['bau'],pol_type='const',**argz)

  # constant policies
  lam = iw.pvecs['baseline'][0] + 1.0
  pol_map_st = lambda (t,s): (lam**t-1.0,s)
  pol_map_t = lambda (t,): (lam**t-1.0,)
  print
  print 'Constant policies:'
  print '  {:30s}: t = {}, s = {}'.format('Baseline',*pol_map_st(policies['baseline const']))
  print '  {:30s}: t = {}, s = {}'.format('Delay',*pol_map_st(policies['delay const']))
  print '  {:30s}: t = {}'.format('Carbon',*pol_map_t(policies['carbon const']))

  # consumption equivalent losses
  print
  print 'Consumption equivalent losses:'
  print '  {:30s}: {}'.format('Business as usual',np.exp(soc_disc0*(welf_base-welf_bau))-1.0)
  print '  {:30s}: {}'.format('Constant',np.exp(soc_disc0*(welf_base-welf_const))-1.0)
  print '  {:30s}: {}'.format('Delay',np.exp(soc_disc0*(welf_base-welf_delay))-1.0)
  print '  {:30s}: {}'.format('Carbon',np.exp(soc_disc0*(welf_base-welf_carbon))-1.0)

  # robustness policies
  generate_policy('baseline',ptype='3step',label='Optimal 3-step Policies')
  generate_policy('gamma1',par_tag='baseline',label='Optimal Policies, $\\gamma = 2x$',ymax_t=2.5)
  generate_policy('gamma2',par_tag='baseline',label='Optimal Policies, $\\gamma = 5x$',ymax_t=5.0)
  generate_policy('beta1',par_tag='baseline',label='Optimal Policies, $\\chi = 0\\%$',ymax_t=1.0)
  generate_policy('beta2',par_tag='baseline',label='Optimal Policies, $\\chi = 20\\%$',ymax_t=1.5)
  generate_policy('middisc',par_tag='baseline',label='Optimal Policies, $\\rho = 0.5\%$',ymax_t=2.0)
  generate_policy('lowdisc',par_tag='baseline',label='Optimal Policies, $\\rho = 0.1\%$',ymax_t=2.0)
  generate_policy('eta1',label='Optimal Policies, $\\eta = 0.35$',ymax_t=3.5)
  generate_policy('eta2',label='Optimal Policies, $\\eta = 0.65$',ymax_t=1.5)
  generate_policy('alpha1',label='Optimal Policies, $\\alpha = 0.03$',ymax_t=1.5)
  generate_policy('alpha2',label='Optimal Policies, $\\alpha = 0.05$',ymax_t=1.0)
  generate_policy('altemit_high',par_tag='altemit',label='Optimal Policies, Alternative Emission $\\rho = 1\%$',ymax_t=1.0)
  generate_policy('altemit_low',par_tag='altemit',label='Optimal Policies, Alternative Emission $\\rho = 0.1\%$',ymax_t=1.0)
  generate_policy('infres',label='Optimal Policies, No Energy Resource',ymax_t=1.0)
  generate_policy('sic4dist',label='Optimal Policies, 4-digit SIC',ymax_t=1.0)
