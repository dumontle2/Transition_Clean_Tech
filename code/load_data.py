from os.path import join
import numpy as np

# paths
data_dir = join('..','data')

# constants
gtc_per_ppmv = 2.12

# load us emmisions data
usem_fname = join(data_dir,'us_emissions.csv')
usem_dat = np.genfromtxt(usem_fname,skip_header=1,delimiter=',')
usem_year = usem_dat[:,0]
usem_emit = usem_dat[:,1]/1e6 # GtC
usem_t0 = usem_year[0]
usem_t1 = usem_year[-1]

# load world emissions data
wem_fname = join(data_dir,'world_emissions.csv')
wem_dat = np.genfromtxt(wem_fname,skip_header=1,delimiter=',')
wem_year = wem_dat[:,0].astype(np.int)
wem_emit = wem_dat[:,1]/1e3 # GtC
wem_t0 = wem_year[0]
wem_t1 = wem_year[1]

# us fraction of emissions
usem_icut = np.nonzero(wem_year==usem_t0)[0][0]
usem_frac = usem_emit/wem_emit[usem_icut:]

# load siple atmospheric carbon data
siple_fname = join(data_dir,'siple_co2.csv')
siple_dat = np.genfromtxt(siple_fname,skip_header=1,delimiter=',')
siple_year = siple_dat[:,0]
siple_co2 = gtc_per_ppmv*siple_dat[:,1]

# load mauna loa atmospheric carbon data
maunaloa_fname = join(data_dir,'maunaloa_co2.csv')
maunaloa_dat = np.genfromtxt(maunaloa_fname,skip_header=1,delimiter=',')
maunaloa_year = maunaloa_dat[:,0]
maunaloa_month_co2 = maunaloa_dat[:,1:]
maunaloa_good = maunaloa_month_co2 != -99.99
maunaloa_month_co2[~maunaloa_good] = 0.0
maunaloa_co2 = gtc_per_ppmv*np.sum(maunaloa_month_co2,axis=1)/np.sum(maunaloa_good,axis=1)

co2_year = np.concatenate([siple_year,maunaloa_year])
co2_vals = np.concatenate([siple_co2,maunaloa_co2])

# construct initial climate state
cidx_base = 14
cyr_base = co2_year[cidx_base]
cval_base = co2_vals[cidx_base]
wem_cidx = np.nonzero(wem_year==cyr_base)[0][0]
wem_emit_cut = wem_emit[wem_cidx:]
wem_year_cut = wem_year[wem_cidx:]
wem_year_min = wem_year_cut[0]
wem_year_max = wem_year_cut[-1]
c_years = len(wem_emit_cut)

# carbon model - golosov et al
phi_L = 0.2
phi_T = 0.7411
phi_e = 0.0313

em_cum_perm = cval_base + phi_L*np.cumsum(wem_emit_cut)
em_cum_trans = phi_T*np.cumsum(np.exp(-phi_e*(wem_year_max-wem_year_cut))*wem_emit_cut)
em_cum_tot = em_cum_perm + em_cum_trans

em_perm0 = cval_base + phi_L*np.sum(wem_emit_cut)
em_trans0 = phi_T*np.sum(np.exp(-phi_e*(wem_year_max-wem_year_cut))*wem_emit_cut)
em_tot0 = em_perm0 + em_trans0

# new patent count distribution
innov_base = 0.247

# load technology gap distribution
loaded_dist_type = None
def load_gap_dist(gap_dist_type):
  global loaded_dist_type,min_diff,max_diff,nd_zero,nd_vals,nd_dist,qd_dist,qc_dist

  loaded_dist_type = gap_dist_type

  print 'Loading {} gap distribution'.format(gap_dist_type)
  if gap_dist_type == 'sic3_mean':
    # based on sic3 data, mean normalized
    pat_fname = join(data_dir,'gap_sample_sic3.csv')
    pat_norm = 42.97/innov_base
  elif gap_dist_type == 'sic4_mean':
    # based on sic4 data, mean normalized
    pat_fname = join(data_dir,'gap_sample_sic4.csv')
    pat_norm = 27.0/innov_base

  # load gap dist data
  pat_data = np.genfromtxt(pat_fname,skip_header=1,delimiter=',').astype(np.int)/pat_norm
  pat_clean = pat_data[:,0]
  pat_dirty = pat_data[:,1]
  pat_diff = np.round(pat_dirty-pat_clean).astype(np.int)

  n_samp = len(pat_diff)
  w_unif = (1.0/n_samp)*np.ones(n_samp)

  min_diff = np.min(pat_diff)
  max_diff = np.max(pat_diff)
  nd_zero = -min_diff

  nd_vals = np.arange(min_diff,max_diff+1)
  nd_dist = np.bincount(pat_diff-min_diff,weights=w_unif)
  qd_dist = np.bincount(pat_diff-min_diff,weights=w_unif*pat_dirty)
  qc_dist = np.bincount(pat_diff-min_diff,weights=w_unif*pat_clean)

# load product count distribution
prod_fname = join(data_dir,'product_counts.csv')
prod_data = np.genfromtxt(prod_fname,skip_header=1,delimiter=',')
prod_counts = np.sum(prod_data[:,1:],axis=0)
prod_data[:,1:] /= prod_counts
prod_binmins = prod_data[:,0].astype(np.int)
prod_dirty_pmf = prod_data[:,1]
prod_clean_pmf = prod_data[:,2]

prod_max = 70 # top of the 60+ bin
prod_bins = np.r_[prod_binmins,prod_max]
prod_binwidth = prod_bins[1:]-prod_bins[:-1]
prod_vals = np.arange(1,prod_max)
prod_clean_pdf = np.concatenate([(v/w)*np.ones(w) for (v,w) in zip(prod_clean_pmf,prod_binwidth)])
prod_dirty_pdf = np.concatenate([(v/w)*np.ones(w) for (v,w) in zip(prod_dirty_pmf,prod_binwidth)])

prod_clean_cdf = np.cumsum(prod_clean_pdf)
prod_dirty_cdf = np.cumsum(prod_dirty_pdf)

clean_prod_frac = prod_counts[0]/np.sum(prod_counts)

# load moment targets
mmt_fname = join(data_dir,'moment_targets.csv')
mmt_data = np.genfromtxt(mmt_fname,delimiter=',')
