#line 2 "estimation_weave.py"

#include <stdio.h>

#define dsamp() (*(pgen_pos++))
#define rand() (*(rand_pos++))

int np,nr,t,active,p,d,rep0,i,rep,zeros,pgain,plose;
double x,r,cut_inno,cut_tauc_p,cut_taud_p,cut_tauc_r,cut_taud_r,cut_test;
int ndv[prod_max];
int* pgen_clean_pos = (int*)pgen_clean_vec;
int* pgen_dirty_pos = (int*)pgen_dirty_vec;
int* pgen_pos;
double* rand_pos = rand_vec;

for (int f = 0; f < nf; f++) {
  active = 1;
  np = nprod_vec[f];
  nr = nrnd_vec[f];
  t = type_vec[f];
  for (i = 0; i < prod_max; i++) ndv[i] = nstep_vec[f*prod_max+i];
  rep0 = rep0_vec[f]; // for entrants

  // dirty = 0, clean = 1
  if (t == 0) {
    x = x_d;
    pgen_pos = pgen_dirty_pos;
  } else {
    x = x_c;
    pgen_pos = pgen_clean_pos;
  }

  pgain = 0;
  plose = 0;

  for (rep = rep0; rep < n_rep; rep++) {
    cut_inno = delt*x*(nr+np);
    cut_tauc_p = cut_inno + delt*double(tau_c)*np;
    cut_taud_p = cut_tauc_p + delt*double(tau_d)*np;
    cut_tauc_r = cut_taud_p + delt*double(tau_c)*nr;
    cut_taud_r = cut_tauc_r + delt*double(tau_d)*nr;

    r = rand();
    if (r < cut_inno) {
      if (np<prod_max) {
        d = dsamp();
        if (((d==0)&&(rand()<0.5)) || ((t==0)&&(d>0)) || ((t==1)&&(d<0))) {
            pgain++;
            ndv[np] = d;
            np++;
        }
      }
    } else if (r < cut_tauc_p) {
      p = (int)floor(rand()*np);
      d = ndv[p];
      if (((d>0)&&(rand()<alpha)) || ((d==1)&&(rand()<0.5*(1.0-alpha))) || (d<=0)) {
        plose++;

        if (t == 0) {
          nr++;
        }
        np--;
        ndv[p] = ndv[np];
      } else {
        ndv[p] = d - 1;
      }
    } else if (r < cut_taud_p) {
      p = (int)floor(rand()*np);
      d = ndv[p];
      if (((d<0)&&(rand()<alpha)) || ((d==-1)&&(rand()<0.5*(1.0-alpha))) || (d>=0)) {
        plose++;

        if (t == 1) {
          nr++;
        }
        np--;
        ndv[p] = ndv[np];
      } else {
        ndv[p] = d + 1;
      }
    } else if (r < cut_tauc_r) {
      if (t == 1) {
        nr--;
      }
    } else if (r < cut_taud_r) {
      if (t == 0) {
        nr--;
      }
    }

    if ((np==0) && (nr==0)) {
      active = 0;
      break;
    }
  }

  nprod_out[f] = np;
  nrnd_out[f] = nr;
  active_out[f] = active;
  zeros = 0;
  for (i = 0; i < np; i++){
    d = ndv[i];
    nstep_out[f*prod_max+i] = d;
    zeros += (d==0);
  }
  zeros_out[f] = zeros;
  pgain_out[f] = pgain;
  plose_out[f] = plose;

  if (t == 0) {
    pgen_dirty_pos = pgen_pos;
  } else {
    pgen_clean_pos = pgen_pos;
  }
}

int pgen_clean_used = (int)(pgen_clean_pos-(int*)pgen_clean_vec);
int pgen_dirty_used = (int)(pgen_dirty_pos-(int*)pgen_dirty_vec);
int rand_used = (int)(rand_pos-rand_vec);

if (pgen_clean_used > pgen_size) {
  printf("pgen_clean_used = %i, pgen_size = %i\n",pgen_clean_used,pgen_size);
}
if (pgen_dirty_used > pgen_size) {
  printf("pgen_dirty_used = %i, pgen_size = %i\n",pgen_dirty_used,pgen_size);
}
if (rand_used > rand_size) {
  printf("rand_used = %i, rand_size = %i\n",rand_used,rand_size);
}
