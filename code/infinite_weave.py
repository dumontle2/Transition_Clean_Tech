# equilibrium solver for clean_tech model

# import libraries
import os
import sys
import itertools as it

import numpy as np
import weave
import scipy.stats as stats

import load_data as ld
import weave_class as wc

#-----------------------------#
# constants                   #
#-----------------------------#

# equilibrium alg params
n_min = -100 # min_diff
n_max = 100 # max_diff
t_max = 3000.0
nt = 2048
v_fact = 1.0
mu_fact = 1.0
val_tol = 1.0e-6
nreps0 = 2000

delt = t_max/nt
nn = n_max-n_min+1

alg_fact1 = 0.95
alg_fact2 = 0.85
alg_fact3 = 0.75

# fixed params
beta0 = 0.1
beta1 = 0.0
beta2 = 0.2

alpha0 = 0.04
alpha1 = 0.03
alpha2 = 0.05

eta0 = 0.50
eta1 = 0.35
eta2 = 0.65

disc0 = 0.01
disc1 = 0.001
disc2 = 0.0001

L = 0.055 # 0.134

# new exhaustible resource params
nu0 = 0.04 # energy input elasticity
nu1 = 0.0 # no exhaustible resource

# policy constraints
m_max = n_max-1
subs_max = 0.99

# params -> (lam0,theta,F_I,F_E,Rmax,ce0)
pvecs = {}
pvecs['baseline'] = [0.0631470, 0.9584000, 0.0020372, 0.0404940, 1.3549e+04, 0.0159820]
pvecs['alpha1'] = [0.0637650, 0.9483600, 0.0020269, 0.0405360, 1.1666e+04, 0.0146230]
pvecs['alpha2'] = [0.0632840, 0.9458200, 0.0019174, 0.0407810, 1.1515e+04, 0.0147890] # 0.00265126
pvecs['eta1'] = [0.0599080, 0.6616200, 0.0070802, 0.0597940, 4.3504e+04, 0.0078974] # 0.00405959
pvecs['eta2'] = [0.0633320, 1.4853000, 0.0013452, 0.0242490, 8930.0000, 0.0142730] # 0.00329983
pvecs['altemit'] = [0.0664150, 0.9085800, 0.0021212, 0.0394570, 1.2196e+04, 0.0213460] # 0.00528114
pvecs['infres'] = [0.0638270, 0.9418200, 0.0020753, 0.0402300, 1.4129e+04, 0.0163290] # 0.00360352
pvecs['sic4dist'] = [0.0649200, 0.9064200, 0.0023711, 0.0391890, 1.985e+04, 0.0209410] # 0.00353397

# default param vec
default_gap_dist = 'sic3_mean'

# falling marginal cost of extraction
ce_growth0 = 0.0 # extraction costs fall at this rate
ce_min0 = 0.0 # as a percentage of the initial level

# alternative emissions model ('energy' emissions) - need decreasing cost of extraction
ce_growth1 = 0.02
ce_min1 = 0.01

# fixed cost fuzzing (for stability)
G_I = 0.2 # 20% heterogeneity in incumbent fixed costs
G_E = 0.2 # 20% heterogeneity in entrant fixed costs

# fixed policies
policies = {}
policies['null'] = [0,0.0]
policies['bau'] = [3.36,0.43] # H = 113.37096 / L = 13803.08015
policies['bau_subs'] = [0,0.43]

# optimal policies
policies_opt = {}
policies_opt['baseline quartic'] = [0.748002742694,-0.00965395404453,0.0014692748031,-3.16394199637e-06,-2.47151128335e-08,0.909828322476,-0.00645633952416,9.74849153419e-05,-3.68488364385e-06,-3.06884261087e-08] # 117.63686
policies_opt['delay quartic'] = [2.92445943913,-0.00222497016368,0.000317584049839,-2.9008366749e-07,-5.71860204156e-10,0.0177685668448,-0.000101767065721,5.05261730291e-05,-1.14257856117e-07,-1.93860686399e-10] # 115.98629
policies_opt['carbon quartic'] = [4.927342338,-0.00515977187102,0.00183745946798,-4.42861171462e-06,-1.54325353031e-08] # 115.75916
policies_opt['baseline 3step'] = [1.06405684817,5.95930252251,0.48708334152,62.7410840381,379.7950485,0.807319045755,0.631321757079,0.00646842813033,32.6932134029,23.8870894964] # 117.42681
policies_opt['baseline const'] = [1.95929614958,0.626530999531] # 116.63689
policies_opt['delay const'] = [2.95646726013,0.632683702611] # 114.34590
policies_opt['carbon const'] = [5.9592957486] # 115.61691

policies_opt['beta1 quartic'] = [0.515240677621,-0.00876503454494,0.000354749463165,-2.8444689022e-07,-1.95271160471e-10,4.29361276831,-0.0146341965185,4.52463339528e-05,-7.4441064536e-07,-2.18219788102e-07] # 118.17032
policies_opt['beta2 quartic'] = [0.819393899211,-0.0114125792751,0.00176042993476,-4.24454153752e-06,-2.21833662916e-08,0.892784506791,-0.00950692495606,0.000190488784895,-4.09222689429e-06,-4.4058136451e-08] # 117.21783
policies_opt['gamma1 quartic'] = [0.738552566401,-0.00798719812632,0.0019234814906,-3.3424408422e-06,-2.56643981282e-08,0.912449058922,-0.00703154586875,0.000100946883142,-5.21731579003e-06,-3.69576669898e-08] # 116.17776
policies_opt['gamma2 quartic'] = [0.827569385905,-0.00609829625275,0.00491626498254,-7.33690875562e-06,-1.19314366746e-08,0.934612248485,-0.0155370825587,2.37319968793e-05,-3.08694269679e-06,-3.56399395344e-08] # 112.12917
policies_opt['middisc quartic'] = [0.391171564787,-0.00874346017821,0.0015164607254,-2.76750611976e-06,-1.76982906885e-08,0.964592271011,-0.00909467446079,0.000102331576565,-4.59094031261e-06,-2.43787186961e-08] # 504.21572
policies_opt['lowdisc quartic'] = [1.70539999139,-0.00656875571498,0.00523940262663,-3.78234861112e-06,-1.33858931014e-08,0.922847941238,-0.0190820853697,0.000240853945144,-2.0447916508e-05,-3.49470265073e-08] # 13402.67610

policies_opt['alpha1 quartic'] = [0.614046332064,-0.0141593726506,0.00218659084791,-3.56951268227e-06,-4.73697839178e-08,0.92471722584,-0.00665368801666,7.17535726928e-05,-2.5957131705e-06,-5.33736511224e-08] # 116.23499
policies_opt['alpha2 quartic'] = [0.769818744412,-0.0057680071594,0.00116586141618,-3.26193009751e-06,-1.33138346968e-08,0.888130327344,-0.00582136852982,9.56265683519e-05,-3.45767714098e-06,-3.19306165088e-08] # 117.29894
policies_opt['eta1 quartic'] = [0.527533636461,-0.0114600655723,0.00106896345012,-7.99992534026e-07,-1.97144051518e-09,0.865043339824,-0.00149648384402,0.00011300328312,-6.18069698815e-06,-2.7182914757e-08] # 118.01317
policies_opt['eta2 quartic'] = [0.905952942784,-0.0100910401107,0.00208376743764,-8.63046594185e-06,-1.91140082793e-08,0.952745676058,-0.0144084972744,2.94507357366e-05,-5.03899432955e-06,-3.30536893848e-08] # 117.99208
policies_opt['altemit_high quartic'] = [0.309693274986,-0.00633865628643,0.000216541242828,-1.54278555842e-06,-1.44373586907e-08,1.05296598191e-06,-0.00440754706368,0.000143333503333,-2.43532916573e-05,-2.205897664e-08] # 129.519186987!
policies_opt['altemit_low quartic'] = [1.18532676736,-0.00970883489391,0.00112369412515,-4.73161531129e-06,-1.23414550334e-08,0.945125418251,-0.0113184875756,0.000224666910152,-3.2756634429e-06,-3.80287647429e-08] # 13347.27943
policies_opt['infres quartic'] = [0.719460164163,-0.00823672747501,0.00146198905279,-4.67193396605e-06,-1.42868717896e-08,0.905855387454,-0.00593553703738,9.14901833653e-05,-4.07471929587e-06,-1.80009690835e-08] # 116.58949
policies_opt['sic4dist quartic'] = [0.71715326748,-0.0119018295586,0.00136848339206,-3.90221312545e-06,-2.5781082404e-08,0.914716254746,-0.00656579950221,0.000106214215128,-4.47226089579e-06,-2.20127780347e-08] # 112.75442

policies_opt['altemit_cost_high quartic'] = [1.05865524347,-0.0169101922857,0.000874610795606,-5.2444813379e-07,-7.58952197505e-10,0.837828187762,-0.00148671834515,3.7169987296e-05,-3.20719300915e-06,-1.32605602472e-08] # 116.753288782 (not reestimated!)
policies_opt['altemit_cost_low quartic'] = [2.33441029501,-0.00267059788302,0.00573082380612,-1.181416151e-06,-2.42849357112e-10,0.870795057011,-0.0182456778681,2.25727065263e-05,-1.72681423349e-06,-4.49057384541e-09] # 13345.63846

# carbon cycle params
Emit0_US = 1.5469
Emit0_world = 8.749
Emit0_other = Emit0_world - Emit0_US
US_emit_ratio0 = Emit0_US/Emit0_world
US_emit_ratio_growth0 = 0.0
growth_other_emit = 0.02
Sbar = 587.0 # GtC
Temp_hell = np.inf # catastrophe

# damage model - exponential
gamma0 = 5.3e-5 # low gamma, exponential damage model
gamma1 = 10.6e-5 # high gamma, exponential damage model
gamma2 = 2.65e-4 # super high gamma, exponential damage model

# emissions model - kappa set to match present day emissions
kappa0 = 3.0640 # JPE R&R

#-----------------------------#
# eq solver setup             #
#-----------------------------#

# step sizes
n_zero = -n_min
n_vec = np.arange(n_min,n_max+1).reshape((1,nn))

# time
t_vec = np.linspace(0.0,t_max,nt).reshape((nt,1))

# other world emissions
Emit_other_vec = Emit0_other*np.exp(growth_other_emit*t_vec)

# tax vector allocation
tax_vec = np.zeros((nt,1))
m_vec = np.zeros((nt,1))
me_vec = np.zeros((nt,1))
subs_vec = np.zeros((nt,1))
rmg_vec = np.zeros((nt,1))
ptax_vec = np.zeros((nt,1))

# distributions over n
mu_vec = np.zeros((nt,nn))
lin_Qhat_d_vec = np.zeros((nt,nn))
lin_Qhat_c_vec = np.zeros((nt,nn))
lQbar_vec = np.zeros((nt,1))
dmu_vec = np.zeros((1,nn))
dlin_Qhat_vec = np.zeros((1,nn))
slamd_vec = np.zeros((1,nn))
slamc_vec = np.zeros((1,nn))

# value functions
vfunc_d_vec = np.zeros((nt,nn))
vfunc_c_vec = np.zeros((nt,nn))

vfunc_d_old = np.zeros((nt,nn))
vfunc_c_old = np.zeros((nt,nn))

vfunc_d_new = np.zeros((nt,nn))
vfunc_c_new = np.zeros((nt,nn))

# growing variables
prof_c_vec = np.zeros((1,nn))
prof_d_vec = np.zeros((1,nn))
flow_d_vec = np.zeros((1,nn))
flow_c_vec = np.zeros((1,nn))

# emmissions
S_vec = np.zeros((nt,1))
S_perm_vec = np.zeros((nt,1))
S_trans_vec = np.zeros((nt,1))

S_world_vec = np.zeros((nt,1))
S_world_perm_vec = np.zeros((nt,1))
S_world_trans_vec = np.zeros((nt,1))

S_follow_vec = np.zeros((nt,1))
S_follow_perm_vec = np.zeros((nt,1))
S_follow_trans_vec = np.zeros((nt,1))

ce_vec = np.zeros((nt,1))
pe_vec = np.zeros((nt,1))
netpe_vec = np.zeros((nt,1))
netpe_growth = np.zeros((nt,1))
R_usage_vec = np.zeros((nt,1))
R_used_vec = np.zeros((nt,1))
R_left_vec = np.zeros((nt,1))
US_emit_ratio = np.zeros((nt,1))

# stored variables
xd_vec = np.zeros((nt,1))
xc_vec = np.zeros((nt,1))
xe_vec = np.zeros((nt,1))
xpd_vec = np.zeros((nt,1))
xpc_vec = np.zeros((nt,1))
evd_vec = np.zeros((nt,1))
evc_vec = np.zeros((nt,1))
ws_vec = np.zeros((nt,1))
wu_vec = np.zeros((nt,1))
labd_vec = np.zeros((nt,1))
labc_vec = np.zeros((nt,1))
labec_vec = np.zeros((nt,1))
labed_vec = np.zeros((nt,1))
taued_vec = np.zeros((nt,1))
tauec_vec = np.zeros((nt,1))
taud_vec = np.zeros((nt,1))
tauc_vec = np.zeros((nt,1))
Gamma_d_vec = np.zeros((nt,1))
Gamma_c_vec = np.zeros((nt,1))
tGamma_d_vec = np.zeros((nt,1))
tGamma_c_vec = np.zeros((nt,1))
nbar_vec = np.zeros((nt,1))
posn_vec = np.zeros((nt,1))
growth_vec = np.zeros((nt,1))
Yd_vec = np.zeros((nt,1))
lY_vec = np.zeros((nt,1))
Ygap_vec = np.zeros((nt,1))
Emit_vec = np.zeros((nt,1))
Emit_world_vec = np.zeros((nt,1))
Emit_follow_vec = np.zeros((nt,1))
Temp_vec = np.zeros((nt,1))
Temp_follow_vec = np.zeros((nt,1))
dlab_vec = np.zeros((nt,1))
lSloss_vec = np.zeros((nt,1))
retd_vec = np.zeros((nt,1))
retc_vec = np.zeros((nt,1))
posd_vec = np.zeros((nt,1))
posc_vec = np.zeros((nt,1))
evcp_vec = np.zeros((nt,1))
labord_vec = np.zeros((nt,1))
laborc_vec = np.zeros((nt,1))
labore_vec = np.zeros((nt,1))
lYbase_vec = np.zeros((nt,1))
lSfact_vec = np.zeros((nt,1))
Omega_vec = np.zeros((nt,1))
Lambda_vec = np.zeros((nt,1))
ltax_fact_vec = np.zeros((nt,1))
vdifft_d_vec = np.zeros((nt,1))
vdifft_c_vec = np.zeros((nt,1))
dirty_labor_vec = np.zeros((nt,1))

# restrict distributions to our binning
nd_dist1 = np.zeros((1,nn))
qd_dist1 = np.zeros((1,nn))
qc_dist1 = np.zeros((1,nn))

iw_loaded_dist_type = None
def load_gap_dist(gap_dist):
  global iw_loaded_dist_type

  print 'Loading iw gap dist {}'.format(gap_dist)

  # load data if needed
  if ld.loaded_dist_type != gap_dist:
    ld.load_gap_dist(gap_dist)

  iw_loaded_dist_type = gap_dist

  min_diff = ld.min_diff
  max_diff = ld.max_diff

  nr_min = np.maximum(n_min,min_diff)
  nr_max = np.minimum(n_max,max_diff)

  nd_dist1[0,nr_min-n_min:nr_max-n_min+1] = ld.nd_dist[nr_min-min_diff:nr_max-min_diff+1]
  if min_diff < n_min: nd_dist1[0,0] += np.sum(ld.nd_dist[:n_min-min_diff+1])
  if max_diff > n_max: nd_dist1[0,-1] += np.sum(ld.nd_dist[n_max-min_diff:])

  qd_dist1[0,nr_min-n_min:nr_max-n_min+1] = ld.qd_dist[nr_min-min_diff:nr_max-min_diff+1]
  if min_diff < n_min: qd_dist1[0,0] += np.sum(ld.qd_dist[:n_min-min_diff+1])
  if max_diff > n_max: qd_dist1[0,-1] += np.sum(ld.qd_dist[n_max-min_diff:])

  qc_dist1[0,nr_min-n_min:nr_max-n_min+1] = ld.qc_dist[nr_min-min_diff:nr_max-min_diff+1]
  if min_diff < n_min: qc_dist1[0,0] += np.sum(ld.qc_dist[:n_min-min_diff+1])
  if max_diff > n_max: qc_dist1[0,-1] += np.sum(ld.qc_dist[n_max-min_diff:])

# set real-time up output
def real_time_setup(rt_year,Rmax):
  tplot_bin = int(np.ceil((rt_year/t_max)*nt))

  import matplotlib.pylab as plt
  import seaborn as sns
  sns.set_style('whitegrid')
  plt.ion()
  (fig,((ax1,ax2,ax5),(ax3,ax4,ax6))) = plt.subplots(nrows=2,ncols=3,figsize=(15,10))

  ax1.set_xlim(0.0,rt_year)
  ax1.set_ylim(0.0,1.0)
  line1d, = ax1.plot([],[])
  line1d.set_xdata(t_vec[:tplot_bin,0])
  line1c, = ax1.plot([],[])
  line1c.set_xdata(t_vec[:tplot_bin,0])
  line1a, = ax1.plot([],[])
  line1a.set_xdata(t_vec[:tplot_bin,0])
  ax1.legend(['dirty','clean','clean_adj'])
  ax1.set_title('gain from innovation')

  ax2.set_xlim(0.0,rt_year)
  ax2.set_ylim(0.0,0.4)
  line2d, = ax2.plot([],[])
  line2d.set_xdata(t_vec[:tplot_bin,0])
  line2c, = ax2.plot([],[])
  line2c.set_xdata(t_vec[:tplot_bin,0])
  ax2.legend(['dirty','clean'])
  ax2.set_title('innovation rates')

  ax3.set_xlim(0.0,rt_year)
  ax3.set_ylim(0.0,0.06)
  line3d, = ax3.plot([],[])
  line3d.set_xdata(t_vec[:tplot_bin,0])
  line3c, = ax3.plot([],[])
  line3c.set_xdata(t_vec[:tplot_bin,0])
  ax3.set_title('entry rates')

  ax4.set_xlim(0.0,rt_year)
  ax4.set_ylim(0.0,1.0)
  line4d, = ax4.plot([],[])
  line4d.set_xdata(t_vec[:tplot_bin,0])
  line4c, = ax4.plot([],[])
  line4c.set_xdata(t_vec[:tplot_bin,0])
  ax4.set_title('valfunc error')

  ax5.set_xlim(0.0,rt_year)
  ax5.set_ylim(0.0,5.0)
  line5, = ax5.plot([],[])
  line5.set_xdata(t_vec[:tplot_bin,0])
  ax5.set_title('')

  ax6.set_xlim(0.0,rt_year)
  ax6.set_ylim(0.0,Rmax)
  line6, = ax6.plot([],[])
  line6.set_xdata(t_vec[:tplot_bin,0])
  ax6.set_title('')

  # draw output
  def draw_data():
    line1d.set_ydata(evd_vec[:tplot_bin,0])
    line1c.set_ydata(evc_vec[:tplot_bin,0])
    line1a.set_ydata(evcp_vec[:tplot_bin,0])
    line2d.set_ydata(taud_vec[:tplot_bin,0])
    line2c.set_ydata(tauc_vec[:tplot_bin,0])
    line3d.set_ydata(taued_vec[:tplot_bin,0])
    line3c.set_ydata(tauec_vec[:tplot_bin,0])
    line4d.set_ydata(vdifft_d_vec[:tplot_bin,0])
    line4c.set_ydata(vdifft_c_vec[:tplot_bin,0])
    line5.set_ydata(np.log(pe_vec[:tplot_bin,0]))
    line6.set_ydata(R_left_vec[:tplot_bin,0])
    plt.draw()

  return draw_data

def plot_output(tcut):
  import matplotlib
  try:
    matplotlib.pyplot
  except AttributeError:
    print 'Loading pylab with Agg backend'
    matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  plt.ioff()

  import seaborn as sns
  sns.set(style='whitegrid',rc={'axes.titlesize': 18,
                                'axes.labelsize': 18,
                                'xtick.labelsize': 18,
                                'ytick.labelsize': 18,
                                'legend.fontsize': 18})

  img_format = 'svg'
  tcut_bin = int(np.ceil((tcut/t_max)*nt))

  # plotting to file
  def plot_series(xvar,series,styles=None,xlabel=None,title=None,legend=None,fname_out=None,title_color='#101010',ymax=None):
    if not styles: styles = ''
    if type(series) not in [list,tuple]: series = [series]
    if type(styles) not in [list,tuple]: styles = [styles]*len(series)

    (fig,ax) = plt.subplots(figsize=(5.5,4.5))
    for (ser,sty) in zip(series,styles): ax.plot(xvar,ser,sty)
    if xlabel: ax.set_xlabel(xlabel)
    if title: ax.set_title(title,color=title_color)
    if legend: ax.legend(legend,loc='best')
    if ymax: ax.set_ylim(ymax=ymax)

    fig.tight_layout()
    if fname_out is not None: fig.savefig(fname_out)

  plot_series(t_vec[:tcut_bin,0],tax_vec[:tcut_bin,0],xlabel='Number of Years',title='Production Tax (m)',fname_out='output/prod_tax.{}'.format(img_format),ymax=1.0)
  plot_series(t_vec[:tcut_bin,0],subs_vec[:tcut_bin,0],xlabel='Number of Years',title='Research Subsidy (s)',fname_out='output/subsidy.{}'.format(img_format),ymax=1.0)
  plot_series(t_vec[:tcut_bin,0],[tauc_vec[:tcut_bin,0],taud_vec[:tcut_bin,0]],xlabel='Number of Years',title='Aggregate Innovation Rates (z)',legend=('Clean','Dirty'),fname_out='output/innovation.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],nbar_vec[:tcut_bin,0],xlabel='Number of Years',title='Mean Step Differential',fname_out='output/nbar.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],1.0-posn_vec[:tcut_bin,0],xlabel='Number of Years',title='Clean Production Lines',fname_out='output/clean_lines.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],[evc_vec[:tcut_bin,0],evd_vec[:tcut_bin,0]],xlabel='Number of Years',title='Expected Gain From Innovation (Ev)',legend=('Clean','Dirty'),fname_out='output/ev_dc.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],[retc_vec[:tcut_bin,0],retd_vec[:tcut_bin,0]],xlabel='Number of Years',title='Expected Return to Innovation',legend=('Clean','Dirty'),fname_out='output/return_dc.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],Emit_vec[:tcut_bin,0],xlabel='Number of Years',title='Emissions (K)',fname_out='output/emission.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],S_follow_vec[:tcut_bin,0],xlabel='Number of Years',title='Carbon Stock (S)',fname_out='output/S_environ.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],lY_vec[:tcut_bin,0],xlabel='Number of Years',title='Production (Y)',fname_out='output/Y_output.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],growth_vec[:tcut_bin,0],xlabel='Number of Years',title='Growth Rate',fname_out='output/growth_rate.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],[tauec_vec[:tcut_bin,0],taued_vec[:tcut_bin,0]],xlabel='Number of Years',title='Entry Rates',legend=('Clean','Dirty'),fname_out='output/entry.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],Temp_follow_vec[:tcut_bin,0],xlabel='Number of Years',title='Temperature Change',fname_out='output/temperature.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],R_usage_vec[:tcut_bin,0],xlabel='Number of Years',title='Energy Resource Usage',fname_out='output/energy_usage.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],R_left_vec[:tcut_bin,0],xlabel='Number of Years',title='Energy Resource Left',fname_out='output/energy_left.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],pe_vec[:tcut_bin,0],xlabel='Number of Years',title='Energy Price',fname_out='output/energy_price.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],np.log(pe_vec[:tcut_bin,0]),xlabel='Number of Years',title='Log Energy Price',fname_out='output/log_energy_price.{}'.format(img_format))
  plot_series(t_vec[:tcut_bin,0],[m_vec[:tcut_bin,0],me_vec[:tcut_bin,0]],xlabel='Number of Years',title='Dirty Step Disadvantage',legend=('Tax Only','Tax + Energy'),fname_out='output/dirty_step_disadv.{}'.format(img_format))

  tf_vec = np.r_[ld.wem_year_cut-2008,t_vec[:,0]]
  tcut_comb = len(ld.wem_year_cut) + tcut_bin
  Emit_comb_vec = np.r_[ld.wem_emit_cut,Emit_follow_vec[:,0]]
  plot_series(tf_vec[:tcut_comb],Emit_comb_vec[:tcut_comb],xlabel='Number of Years',title='World Emissions (K)',fname_out='output/world_emission.{}'.format(img_format))

  plt.close('all')

def quick_print(vec):
  return '['+','.join(map(str,vec))+']'

#-----------------------------+
# eq solve functions          |
#-----------------------------+

def sim_pvec(pvec=pvecs['baseline'],output=False,graph=False,real_time=False,disc=disc0,gamma=gamma0,beta=beta0,soc_disc=disc0,alpha=alpha0,eta=eta0,nu=nu0,gap_dist=default_gap_dist,kappa=kappa0,nreps=nreps0,ce_growth=ce_growth0,ce_min=ce_min0,emit_model='output',US_emit_ratio_growth=US_emit_ratio_growth0,graph_rep=25,graph_year=500.0,real_time_year=300.0):
  if iw_loaded_dist_type != gap_dist:
    load_gap_dist(gap_dist)

  (lam0,theta,F_I,F_E,Rmax,ce0) = pvec
  lam = 1.0+lam0
  prof = 1.0-1.0/lam
  lbar = np.log(lam)

  # can it get leapfrogged?
  alpha_d_vec = alpha*(n_vec<0)
  alpha_c_vec = alpha*(n_vec>0)

  # discounting
  disc_vec = np.exp(-disc*t_vec)
  soc_disc_vec = np.exp(-soc_disc*t_vec)

  # steady state values (for either outcome)
  ev_ss = prof/(disc+L*theta*np.power(eta,eta)*np.power((1.0-eta)/F_E,1.0-eta))
  x_ss = theta*np.power(eta*F_E/(1.0-eta),eta)
  tau_ss = theta*np.power(eta,eta)*np.power((1.0-eta)/F_E,1.0-eta)*(L-F_I+F_E)
  ws_ss = ev_ss*theta*np.power(eta,eta)*np.power((1.0-eta)/F_E,1.0-eta)

  # initialize value functions
  taubar = 0.15
  vfunc_d_vec[:,:] = prof*(n_vec>0)/(disc+taubar)
  vfunc_c_vec[:,:] = prof*(n_vec<0)/(disc+taubar)

  # initial energy price level
  netpe_level = 0.0275 # initial guess for level
  nuhat = ((1.0-nu)**(1.0-nu)*nu**nu)**(1.0/nu) if nu > 0.0 else 1.0
  netpe_vec[:,0] = netpe_level*np.exp(disc*t_vec[:,0])
  ce_vec[:,0] = ce0*np.maximum(ce_min,np.exp(-ce_growth*t_vec[:,0]))
  pe_vec[:,0] = netpe_vec[:,0] + ce_vec[:,0]
  me_vec[:,0] = m_vec[:,0] + nu*(-np.log(nuhat)+np.log(pe_vec[:,0]))/lbar
  US_emit_ratio[:,0] = US_emit_ratio0*np.exp(-US_emit_ratio_growth*t_vec[:,0])

  # initial interest rate
  rmg_vec[:,0] = disc

  # initialize distribution vectors
  mu_vec[:,:] = nd_dist1[:,:]

  ###############################
  ## set up climate and policy ##
  ###############################

  # aggregate productivity (not tax-adjusted)
  lQbar_vec[0,0] = 0.0

  # lin_Qhat - direct linear integral of q's used in production - from the patent count data
  # assumption - all product lines start at q=1, then patents modify them
  lin_Qhat_c_vec[0,:] = np.power(lam,qc_dist1[0,:])
  lin_Qhat_c_vec[0,:] /= np.sum(lin_Qhat_c_vec[0,:])
  lin_Qhat_d_vec[0,:] = np.power(lam,qd_dist1[0,:])
  lin_Qhat_d_vec[0,:] /= np.sum(lin_Qhat_d_vec[0,:])

  # initialize emissions
  S_perm_vec[0,0] = ld.em_perm0
  S_trans_vec[0,0] = ld.em_trans0

  S_world_perm_vec[0,0] = ld.em_perm0
  S_world_trans_vec[0,0] = ld.em_trans0

  S_follow_perm_vec[0,0] = ld.em_perm0
  S_follow_trans_vec[0,0] = ld.em_trans0

  Emit_other_vec[0,0] = Emit0_other

  R_used_vec[0,0] = 0.0

  ##############################
  ## set up internal solver   ##
  ##############################

  # localize carbon cycle parameters
  phi_e = ld.phi_e
  phi_L = ld.phi_L
  phi_T = ld.phi_T

  # weave setup
  args_alg = ['delt','nn','nt','n_zero','v_fact','mu_fact']
  args_pol = ['m_vec','subs_vec']
  args_par = ['theta','eta','L','F_I','F_E','G_I','G_E','alpha','prof','disc']
  args_eqv = ['ev_ss','tau_ss','x_ss','ws_ss']
  args_inp = ['n_vec','dmu_vec','flow_d_vec','flow_c_vec','vfunc_d_new','vfunc_c_new','rmg_vec']
  args_out = ['mu_vec','vfunc_d_vec','vfunc_c_vec','evd_vec','evc_vec','xd_vec','xc_vec','xpd_vec','xpc_vec','taud_vec','tauc_vec','ws_vec','wu_vec','taued_vec','tauec_vec','retd_vec','retc_vec','posd_vec','posc_vec','evcp_vec','xe_vec','labd_vec','labc_vec','labec_vec','labed_vec']
  args_ext = ['alpha_c_vec','alpha_d_vec','prof_c_vec','prof_d_vec']
  args_emt = ['gamma','nu','lam','lam0','lbar','US_emit_ratio','growth_other_emit','me_vec','ce_vec','pe_vec','R_used_vec','Emit_other_vec','slamd_vec','slamc_vec','R_usage_vec','dirty_labor_vec']
  args = args_alg + args_pol + args_par + args_eqv + args_inp + args_out + args_ext + args_emt
  code = wc.generate_code('dynamics_class.cpp',args,locals(),globals(),class_name='EqSolver',inst_name='eqs',store='last_code.cpp')

  if real_time:
    real_time_update = real_time_setup(real_time_year,Rmax)

  # solve value functions with forward-backward iteration
  for rep in range(nreps):
    # to check convergence later
    vfunc_d_old[:,:] = vfunc_d_vec[:,:]
    vfunc_c_old[:,:] = vfunc_c_vec[:,:]
    netpe_level_old = netpe_level

    # this is an art
    if rep < 10:
      v_fact = 1.0
      mu_fact = 1.0
    elif rep < 100:
      v_fact = alg_fact1
      mu_fact = alg_fact1
    elif rep < 500:
      v_fact = alg_fact2
      mu_fact = alg_fact2
    elif rep < 750:
      v_fact = alg_fact3
      mu_fact = alg_fact3
    else:
      v_fact = 0.5
      mu_fact = 0.5

    # run c++ compiled code
    weave.inline(code,arg_names=args,extra_compile_args=['-O2'])

    # update effective interest rate - 0.5 bound is for stability, not even close to holding in real scenarios
    lSloss_vec[:,0] = np.log(1.0-np.minimum(0.5,ws_vec[:,0]*(labc_vec[:,0]+labec_vec[:,0])*subs_vec[:,0]*beta))
    rmg_vec[:,0] = disc + np.r_[np.diff(lSloss_vec[:,0])/np.diff(t_vec[:,0]),0.0]

    # update energy price
    R_used = R_used_vec[nt-1,0]
    netpe_level *= (R_used/Rmax)**(0.25) # exponent is smooth adjustment parameter
    wu_growth = np.r_[np.diff(np.log(wu_vec[:,0]))/np.diff(t_vec[:,0]),0.0]
    netpe_growth[:,0] = rmg_vec[:,0] + wu_growth
    netpe_vec[:,0] = netpe_level*np.exp(netpe_growth[:,0]*t_vec[:,0])
    pe_vec[:,0] = netpe_vec[:,0] + ce_vec[:,0]
    me_vec[:,0] = m_vec[:,0] + nu*(-np.log(nuhat)+np.log(pe_vec[:,0]))/lbar
    R_left_vec[:,0] = Rmax - R_used_vec[:,0]

    # check for convergence
    vdiff_d = np.abs(vfunc_d_old-vfunc_d_new)
    vdiff_c = np.abs(vfunc_c_old-vfunc_c_new)
    vdifft_d_vec[:,0] = np.max(vdiff_d,axis=1)
    vdifft_c_vec[:,0] = np.max(vdiff_c,axis=1)
    verr_d = np.max(vdifft_d_vec[:,0])
    verr_c = np.max(vdifft_c_vec[:,0])
    verr = np.maximum(verr_d,verr_c)
    netpe_err = np.abs(netpe_level-netpe_level_old)/netpe_level_old
    Rerr = np.abs(R_used-Rmax)/Rmax

    if output and rep%25 == 0:
      print '{:4d}: {: >17} ({:4.0f}), {: >17} ({:4.0f}), {: >17}'.format(rep,verr_d,t_vec[np.argmax(vdifft_d_vec[:,0]),0],verr_c,t_vec[np.argmax(vdifft_c_vec[:,0]),0],netpe_level)

    if real_time and rep%graph_rep == 0:
      real_time_update()

    if verr < val_tol and (((netpe_err < val_tol) and (Rerr < 1e-6)) or (netpe_level < 1e-15)):
      if output: print '{:4d}: {: >17} ({:4.0f}), {: >17} ({:4.0f}), {: >17}'.format(rep,verr_d,t_vec[np.argmax(vdifft_d_vec[:,0]),0],verr_c,t_vec[np.argmax(vdifft_c_vec[:,0]),0],netpe_level)
      if real_time: real_time_update()
      break

  ##############################
  ## evaluate climate outcome ##
  ##############################

  for t in range(nt):
    # policy
    mv = m_vec[t,0]
    me = me_vec[t,0]
    subst = subs_vec[t,0]

    mv0 = np.floor(me)
    mv1 = mv0 + 1
    mf = 1.0 - (me-mv0)

    mt0 = int(mv0 + n_zero)
    mt1 = int(mv1 + n_zero)

    # mu aggregates
    mu_0 = mu_vec[t,mt0] # fractional profits for clean leader
    mu_1 = mu_vec[t,mt1] # fractional profits for dirty leader
    mu_dp = np.sum(mu_vec[t,mt1+1:])
    mu_cn = np.sum(mu_vec[t,:mt0])
    mu_d = mu_dp + mu_1
    mu_c = mu_cn + mu_0

    # probability of productivity improvement
    mu0_d = np.sum(mu_vec[t,n_zero:])
    mu0_c = np.sum(mu_vec[t,:n_zero+1])
    Gamma_d = mu0_d + alpha*(1.0-mu0_d)
    Gamma_c = mu0_c + alpha*(1.0-mu0_c)

    # taxes
    tax = (lam**mv)-1.0
    itax = 1.0/(1.0+tax)

    prof_frac_c = lam**mf # at mv0
    prof_frac_d = lam**(1.0-mf) # at mv1

    # innovation
    tau_d = taud_vec[t,0]
    tau_c = tauc_vec[t,0]
    tau = tau_d + tau_c

    # aggregates
    ce = ce_vec[t,0]
    pe = pe_vec[t,0]
    Omega = itax*((1.0-nu)+nu*(ce/pe))*((1.0/prof_frac_d)*mu_1+(1.0/lam)*mu_dp) + ((1.0/prof_frac_c)*mu_0+(1.0/lam)*mu_cn)
    Omega_d = (1.0-nu)*((1.0/prof_frac_d)*mu_1+(1.0/lam)*mu_dp)*itax;
    Lambda = ((1.0/prof_frac_c)**mu_0)*((1.0/prof_frac_d)**mu_1)*((1.0/lam)**(mu_cn+mu_dp))
    growth = lbar*(tau_c*Gamma_c+tau_d*Gamma_d)
    nbar = np.sum(n_vec[0,:]*mu_vec[t,:])
    wu = Omega

    if emit_model == 'output':
      # emissions model - dirty production
      lin_Qd = (1.0/prof_frac_d)*lin_Qhat_d_vec[t,mt1] + (1.0/lam)*np.sum(lin_Qhat_d_vec[t,mt1+1:])
      Y_d = lin_Qd/(wu*lam**me)
      Emit = kappa*Y_d
    elif emit_model == 'energy':
      # emissions model - energy resource
      Emit = R_usage_vec[t,0] # energy resource = emissions

    # carbon law of motion
    Emit_other = Emit_other_vec[t,0]
    Emit_world = Emit + Emit_other
    Emit_follow = Emit/US_emit_ratio[t,0]

    S_perm = S_perm_vec[t,0]
    S_trans = S_trans_vec[t,0]
    S = S_perm + S_trans

    S_world_perm = S_world_perm_vec[t,0]
    S_world_trans = S_world_trans_vec[t,0]
    S_world = S_world_perm + S_world_trans

    S_follow_perm = S_follow_perm_vec[t,0]
    S_follow_trans = S_follow_trans_vec[t,0]
    S_follow = S_follow_perm + S_follow_trans

    Temp = 3.0*np.log2(S_world/Sbar)
    Temp_follow = 3.0*np.log2(S_follow/Sbar)

    # damage model - exponential
    lS_fact = -gamma*(S_follow-Sbar)
    S_fact = np.exp(lS_fact)

    # subsidy losses
    ws = ws_vec[t,0]
    lSloss = lSloss_vec[t,0]
    Sloss = np.exp(lSloss)

    lQbar = lQbar_vec[t,0] # not tax-adjusted
    if me < 0.0:
      ltax_factor = lbar*np.sum(np.minimum(-me,-n_vec[0,:n_zero])*mu_vec[t,:n_zero])
    else:
      ltax_factor = lbar*np.sum(np.minimum(me,n_vec[0,n_zero+1:])*mu_vec[t,n_zero+1:])
    lQbar_tadj = lQbar - ltax_factor
    Qbar = np.exp(lQbar)
    Qbar_tadj = np.exp(lQbar_tadj)
    Ybase = Sloss*(Lambda/Omega)*Qbar_tadj
    Y = S_fact*Ybase
    lY_base = np.log(Lambda/Omega) + lSloss + lQbar_tadj
    lY = lY_base + lS_fact

    # extra stats
    Q_gap = np.exp(-np.sum(np.minimum(mv,n_vec[0,n_zero+1:])*mu_vec[t,n_zero+1:]))
    Omega_null = (1.0/prof_frac_c)*mu_0 + (1.0/prof_frac_d)*mu_1 + (1.0/lam)*(mu_cn+mu_dp)
    Y_gap = Q_gap*(Omega/Omega_null)

    # storage
    nbar_vec[t,0] = nbar
    posn_vec[t,0] = mu_d
    Gamma_c_vec[t,0] = Gamma_c
    Gamma_d_vec[t,0] = Gamma_d
    growth_vec[t,0] = growth
    wu_vec[t,0] = wu
    lY_vec[t,0] = lY
    Ygap_vec[t,0] = Y_gap
    S_vec[t,0] = S
    S_world_vec[t,0] = S_world
    S_follow_vec[t,0] = S_follow
    Emit_vec[t,0] = Emit
    Emit_world_vec[t,0] = Emit_world
    Emit_follow_vec[t,0] = Emit_follow
    Temp_vec[t,0] = Temp
    Temp_follow_vec[t,0] = Temp_follow
    lYbase_vec[t,0] = lY_base
    lSfact_vec[t,0] = lS_fact
    Omega_vec[t,0] = Omega
    Lambda_vec[t,0] = Lambda
    ltax_fact_vec[t,0] = ltax_factor
    R_left_vec[t,0] = Rmax - R_used_vec[t,0]

    if t < nt-1:
      # deterministic growth in aggregate effective quality lQbar
      lQbar_vec[t+1,0] = lQbar + delt*growth

      # update lin_Qhat_d distribution
      dlin_Qhat_vec[0,0] = tau_c*lin_Qhat_d_vec[t,0] + tau_c*lin_Qhat_d_vec[t,1] - tau*lin_Qhat_d_vec[t,0]
      dlin_Qhat_vec[0,1:-1] = (1.0-alpha_d_vec[0,:-2])*tau_d*lam*lin_Qhat_d_vec[t,:-2] + (1.0-alpha_c_vec[0,2:])*tau_c*lin_Qhat_d_vec[t,2:] - tau*lin_Qhat_d_vec[t,1:-1]
      dlin_Qhat_vec[0,-1] = tau_d*lam*lin_Qhat_d_vec[t,-1] + tau_d*lam*lin_Qhat_d_vec[t,-2] - tau*lin_Qhat_d_vec[t,-1]
      dlin_Qhat_vec[0,n_zero-1] += alpha*tau_c*np.sum(lin_Qhat_d_vec[t,n_zero+1:])
      dlin_Qhat_vec[0,n_zero+1] += alpha*tau_d*lam*np.sum(lin_Qhat_c_vec[t,:n_zero])
      lin_Qhat_d_vec[t+1,:] = lin_Qhat_d_vec[t,:] + delt*dlin_Qhat_vec[0,:]

      # update lin_Qhat_c distribution
      dlin_Qhat_vec[0,0] = tau_c*lam*lin_Qhat_c_vec[t,0] + tau_c*lam*lin_Qhat_c_vec[t,1] - tau*lin_Qhat_c_vec[t,0]
      dlin_Qhat_vec[0,1:-1] = (1.0-alpha_d_vec[0,:-2])*tau_d*lin_Qhat_c_vec[t,:-2] + (1.0-alpha_c_vec[0,2:])*tau_c*lam*lin_Qhat_c_vec[t,2:] - tau*lin_Qhat_c_vec[t,1:-1]
      dlin_Qhat_vec[0,-1] = tau_d*lin_Qhat_c_vec[t,-1] + tau_d*lin_Qhat_c_vec[t,-2] - tau*lin_Qhat_c_vec[t,-1]
      dlin_Qhat_vec[0,n_zero-1] += alpha*tau_c*lam*np.sum(lin_Qhat_d_vec[t,n_zero+1:])
      dlin_Qhat_vec[0,n_zero+1] += alpha*tau_d*np.sum(lin_Qhat_c_vec[t,:n_zero])
      lin_Qhat_c_vec[t+1,:] = lin_Qhat_c_vec[t,:] + delt*dlin_Qhat_vec[0,:]

      # growth of environmental stock
      S_perm_vec[t+1,0] = S_perm + delt*phi_L*Emit
      S_trans_vec[t+1,0] = (1.0-delt*phi_e)*S_trans + delt*phi_T*Emit

      S_world_perm_vec[t+1,0] = S_world_perm + delt*phi_L*Emit_world
      S_world_trans_vec[t+1,0] = (1.0-delt*phi_e)*S_world_trans + delt*phi_T*Emit_world

      S_follow_perm_vec[t+1,0] = S_follow_perm + delt*phi_L*Emit_follow
      S_follow_trans_vec[t+1,0] = (1.0-delt*phi_e)*S_follow_trans + delt*phi_T*Emit_follow

      # deterministic growth of rest of the world emissions
      Emit_other_vec[t+1,0] = Emit_other*(1.0+delt*growth_other_emit)

  # continuation welfare
  welfare_0 = delt*np.sum(lY_vec[:,0]*soc_disc_vec[:,0])

  # welfare decomposition
  welfare_lYprod = delt*np.sum((np.log(Lambda_vec/Omega_vec)-ltax_fact_vec)[:,0]*soc_disc_vec[:,0])
  welfare_lSloss = delt*np.sum(lSloss_vec[:,0]*soc_disc_vec[:,0])
  welfare_lQgrow = delt*np.sum(lQbar_vec[:,0]*soc_disc_vec[:,0])
  welfare_lSfact = delt*np.sum(lSfact_vec[:,0]*soc_disc_vec[:,0])

  if verr >= val_tol:
    conv = 'noconv'
    welfare_0 = -np.inf
  elif Temp_vec[nt-1,0] >= Temp_hell:
    conv = 'catastrophe'
    welfare_0 = -np.inf
  else:
    conv = 'clean'
    welfare_0 += soc_disc_vec[nt-1,0]*(lYbase_vec[nt-1,0]+growth_vec[nt-1,0]/soc_disc-gamma*(S_perm_vec[nt-1,0]+S_trans_vec[nt-1,0]*(soc_disc/(soc_disc+ld.phi_e))-Sbar))/soc_disc

  # final aggregates
  tax_vec[:,0] = np.power(lam,m_vec[:,0]) - 1.0
  ptax_vec[:,0] = tax_vec[:,0]/(1.0+tax_vec[:,0])

  #########################
  ## output tons of info ##
  #########################

  # textual output
  if output:
    print 'status = {}'.format(conv)
    print 'soc_disc = {}'.format(soc_disc)
    print 'disc = {}'.format(disc)
    print 'gamma = {}'.format(gamma)
    print 'beta = {}'.format(beta)
    print 'alpha = {}'.format(alpha)
    print 'eta = {}'.format(eta)
    print 'delt = {}'.format(delt)
    print 'sum(mu_n) = {}'.format(np.sum(mu_vec[-1,:]))
    print 'lQbar = {}'.format(lQbar_vec[-1,0])
    print 'mean n = {}'.format(nbar_vec[-1,0])
    print 'pos_n = {}'.format(posn_vec[-1,0])
    print 'Y = {}'.format(np.exp(lY_vec[-1,0]))
    print 'S = {}'.format(S_vec[-1,0])
    print 'netpe_level = {}'.format(netpe_level)
    print 'R_used = {}/{}'.format(R_used_vec[-1,0],Rmax)
    print 'Emit0 = {}/{}'.format(Emit_vec[0,0],Emit0_US)
    print 'lQbar0 = {}'.format(lQbar_vec[0,0])
    print 'Temp = {}'.format(Temp_vec[-1,0])
    print 'x   = {:<8.5f} {:<8.5f}'.format(xd_vec[0,0],xc_vec[0,0])
    print 'tau = {:<8.5f} {:<8.5f}'.format(taud_vec[0,0],tauc_vec[0,0])
    print 'entry = {:<8.5f} {:<8.5f}'.format(taued_vec[0,0],tauec_vec[0,0])
    print 'growth = {:.4f}'.format(growth_vec[0,0])
    print 'ws = {:.4f}'.format(ws_vec[0,0])
    print 'wu = {:.4f}'.format(wu_vec[0,0])
    print 'welfare = {:15.10f}'.format(welfare_0)
    print 'W_prod = {:.5f}, W_subs = {:.5f}, W_growth = {:.5f}, W_clim = {:.5f}'.format(welfare_lYprod,welfare_lSloss,welfare_lQgrow,welfare_lSfact)
    print

  # graphical output
  if graph:
    plot_output(graph_year)

  return (rep,conv,welfare_0)

def calculate_policy(pol_in,pol_type='const',pol_data={},graph=False,graph_year=300.0,lam=None):
  if pol_type == 'const':
    (m,subs) = pol_in

    m_vec[:,0] = m
    subs_vec[:,0] = subs
  elif pol_type == 'delay': # null policy for fifty years, then constant
    (m,subs) = pol_in

    delay_bin = int(np.ceil((50.0/t_max)*nt))

    m_vec[0:delay_bin,0] = 0
    m_vec[delay_bin:,0] = m

    subs_vec[0:delay_bin,0] = 0.0
    subs_vec[delay_bin:,0] = subs
  elif pol_type == 'carbon': # only carbon tax
    (m,) = pol_in

    m_vec[:,0] = m
    subs_vec[:,0] = 0.0
  elif pol_type == '3step': # three stage policy
    (m1,m2,m3,m_time1,m_time2,subs1,subs2,subs3,subs_time1,subs_time2) = pol_in

    m_bin1 = int(np.ceil((m_time1/t_max)*nt))
    m_bin2 = int(np.ceil(((m_time1+m_time2)/t_max)*nt))
    subs_bin1 = int(np.ceil((subs_time1/t_max)*nt))
    subs_bin2 = int(np.ceil(((subs_time1+subs_time2)/t_max)*nt))

    m_vec[:m_bin1,0] = m1
    m_vec[m_bin1:m_bin2,0] = m2
    m_vec[m_bin2:,0] = m3

    subs_vec[:subs_bin1,0] = subs1
    subs_vec[subs_bin1:subs_bin2,0] = subs2
    subs_vec[subs_bin2:,0] = subs3
  elif pol_type == 'quartic':
    (m0,m1,m2,m3,m4,subs0,subs1,subs2,subs3,subs4) = pol_in

    t2_vec = t_vec[:,0]*t_vec[:,0]
    t3_vec = t2_vec*t_vec[:,0]
    t4_vec = t3_vec*t_vec[:,0]

    m_vec[:,0] = np.maximum(0.0,m0+m1*t_vec[:,0]+m2*t2_vec+m3*t3_vec+m4*t4_vec)
    subs_vec[:,0] = np.maximum(0.0,subs0+subs1*t_vec[:,0]+subs2*t2_vec+subs3*t3_vec+subs4*t4_vec)
  elif pol_type == 'quartic_delay': # quartic policy + 50 year delay
    (m0,m1,m2,m3,m4,subs0,subs1,subs2,subs3,subs4) = pol_in

    delay_bin = int(np.ceil((50.0/t_max)*nt))
    d1_vec = t_vec[:-delay_bin,0]
    d2_vec = d1_vec*d1_vec
    d3_vec = d2_vec*d1_vec
    d4_vec = d3_vec*d1_vec

    m_vec[0:delay_bin,0] = 0
    subs_vec[0:delay_bin,0] = 0.0

    m_vec[delay_bin:,0] = np.maximum(0.0,m0+m1*d1_vec+m2*d2_vec+m3*d3_vec+m4*d4_vec)
    subs_vec[delay_bin:,0] = np.maximum(0.0,subs0+subs1*d1_vec+subs2*d2_vec+subs3*d3_vec+subs4*d4_vec)
  elif pol_type == 'quartic_carbon': # quartic policy + carbon only
    (m0,m1,m2,m3,m4) = pol_in

    t2_vec = t_vec[:,0]*t_vec[:,0]
    t3_vec = t2_vec*t_vec[:,0]
    t4_vec = t3_vec*t_vec[:,0]

    m_vec[:,0] = np.maximum(0.0,m0+m1*t_vec[:,0]+m2*t2_vec+m3*t3_vec+m4*t4_vec)
    subs_vec[:,0] = 0.0
  else:
    print 'Policy not {} found.'.format(pol_type)
    raise Exception

  # impose bounds
  m_vec[m_vec[:,0]>m_max,0] = m_max
  m_vec[m_vec[:,0]<0,0] = 0.0

  subs_vec[subs_vec[:,0]>subs_max,0] = subs_max
  subs_vec[subs_vec[:,0]<0.0,0] = 0.0

  # possibly round
  if pol_data.get('round',False): m_vec[:,0] = np.round(m_vec[:,0])

  # possibly cut off
  if 'cutoff' in pol_data:
    icut = int(np.ceil((pol_data['cutoff']/t_max)*nt))
    m_vec[icut:] = 0.0
    subs_vec[icut:] = 0.0

  # possibly calculate real rate
  if lam:
    tax_vec[:,0] = np.power(lam,m_vec[:,0]) - 1.0
    ptax_vec[:,0] = tax_vec[:,0]/(1.0+tax_vec[:,0])

def simulate_type(pol_in,pol_type='const',pol_data={},rep=0,**kwargs):
  calculate_policy(pol_in,pol_type=pol_type,pol_data=pol_data)
  return sim_pvec(**kwargs)

# wrap objective function
def sim_robust_type(pol,pol_type='const',pol_data={},**kwargs):
  if pol_type == 'const':
    lower = np.array([0.0,0.0])
    upper = np.array([m_max,0.95])
  elif pol_type == 'delay':
    lower = np.array([0.0,0.0])
    upper = np.array([m_max,0.95])
  elif pol_type == 'carbon':
    lower = np.array([0.0])
    upper = np.array([m_max])
  elif pol_type == '3step':
    lower = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])
    upper = np.array([m_max,m_max,m_max,t_max,t_max,subs_max,subs_max,subs_max,t_max,t_max])
  elif pol_type == 'quartic':
    lower = np.array([-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf])
    upper = np.array([np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])
  elif pol_type == 'quartic_delay':
    lower = np.array([-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf])
    upper = np.array([np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])
  elif pol_type == 'quartic_carbon':
    lower = np.array([-np.inf,-np.inf,-np.inf,-np.inf,-np.inf])
    upper = np.array([np.inf,np.inf,np.inf,np.inf,np.inf])
  else:
    print 'Policy not {} found.'.format(pol_type)
    raise Exception

  if np.any((pol<lower)|(pol>upper)):
    welf = np.nan
  else:
    try:
      (rep,conv,welf) = simulate_type(pol,pol_type=pol_type,pol_data=pol_data,**kwargs)
    except Exception as e:
      (etype,evalue,etraceback) = sys.exc_info()
      print etype
      print evalue
      print etraceback
      welf = np.nan

  return -welf

def sim_robust_obj(**kwargs):
  return lambda pol,**ikwargs: sim_robust_type(pol,**dict(kwargs,**ikwargs))

# simple annealing optimizer
def anneal0_pol(f,x0,scale=0.25,end_scale=1.0e-5,down_scale=0.5,stay_steps=50,maxiter=1000000,seed=None,mint=[],graph=False,graph_year=200.0,**kwargs):
  rs = np.random.RandomState()
  if seed:
    rs.seed(seed)

  n = len(x0)
  xp = np.array(x0)
  fp = f(x0,**kwargs)
  if graph: plot_output(graph_year)
  xmin = xp
  if mint != None: xmin[mint] = np.round(xmin[mint]) # to deal with discrete m
  fmin = fp
  stay = 0
  print 'SCALE -> {}'.format(scale)
  for i in range(maxiter):
    xp = xmin*np.exp(np.random.normal(loc=0.0,scale=scale,size=n))
    if mint != None:
      for mi in mint:
        if xmin[mi] <= 1 and np.random.rand() < scale/2.0:
          xp[mi] = 1-xmin[mi]
        elif xmin[mi] == 1 and np.random.rand() < scale/2.0:
          xp[mi] = 2
        elif xmin[mi] == 2 and np.random.rand() < scale/2.0:
          xp[mi] = 1
    fp = f(xp,rep=i+1,**kwargs)
    if fp < fmin and ~np.isnan(fp):
      xmin = xp
      if mint != None: xmin[mint] = np.round(xmin[mint]) # to deal with discrete m
      fmin = fp
      print 'MAX -> {} # {}'.format(quick_print(xmin),-fmin)
      if graph: plot_output(graph_year)
      stay = 0
    else:
      if ~np.isnan(fp):
        stay += 1

    if stay == stay_steps:
      scale *= down_scale
      print 'SCALE -> {}'.format(scale)
      stay = 0

    if scale < end_scale:
      break

  print quick_print(xmin)+' # {:<8.5f}'.format(-fmin)

  return (xmin,fmin)

def gridsearch_const(m_min=0,m_max=10,n_subs=10,pol_type='const',**kwargs):
  m_grid = np.arange(m_min,m_max+1)
  subs_grid = np.linspace(0.0,0.95,n_subs)

  fmin = np.inf

  if pol_type == 'const' or pol_type == 'delay':
    pol_grid = it.product(m_grid,subs_grid)
  elif pol_type == 'carbon':
    pol_grid = zip(m_grid)

  for pol in pol_grid:
    fp = sim_robust_type(pol,pol_type=pol_type,**kwargs)
    if fp < fmin:
      fmin = fp
      pol_opt = pol

  print pol_opt
  print fmin

  return (pol_opt,fmin)
