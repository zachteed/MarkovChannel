#include "cost.hpp"
#include "intel_ode.h"
#include <math.h>

using namespace MarkovChannel;
using namespace std;


void print_matrix(double* A, int n, int m) {
  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++) {
      printf("%8.4f\t", A[i*m+j]);
    }
    printf("\n");
  }
}


inline double peak(int n, double* x)
{
  int idx = cblas_idamax(n, x, 1);
  double res = x[idx]; return res;
}

inline double tau(int n, double* x, double v1, double v2)
{
  int idx = cblas_idamax(n, x, 1);
  cblas_dscal(n, 1/x[idx], x, 1);
  double res;

  if (v2 > v1) {
    int i1=0, i2=0; for(int t=0; t<idx; t++) {
      if (!i1 && x[t+1] > v1 && x[t] <= v1) i1=t;
      if (!i2 && x[t+1] > v2 && x[t] <= v2) i2=t;
    }
    res = i2 - i1;
  } else {
    int i1=0, i2=0; for(int t=idx; t<(n-1); t++) {
      if (!i1 && x[t+1] < v1 && x[t] >= v1) i1=t;
      if (!i2 && x[t+1] < v2 && x[t] >= v2) i2=t;
    }
    res = i2 - i1;
  }

  return res;
}

void rhs(int *n, double *t, double *y, double *f) {
  cblas_dgemv(CblasRowMajor, CblasNoTrans,
    *n, *n, 1.0, &y[*n], *n, y, 1, 0.0, f, 1);
}

void jac(int *n, double *t, double *y, double *a) {
  mkl_domatcopy('C', 'T', *n, *n, 0.0, &y[*n], *n, a, *n);
}



int step_ode(Model::Model& m, vector<Step>& steps,
  double *y, vector<double>& out)
{
  int N = m.n_states();
  double *Q = &y[N];

  int ipar[128];
  int *kd = (int*) malloc(N*sizeof(int));

  int n, ierr, cnt = max(13*N, (7+2*N)*N)*sizeof(double);
  double *dpar = (double*) malloc(cnt), t, t_end, h, hm, ep, tr;

  /*************************** ode params *****************************/

  hm=1.e-12; /* minimal step size for the methods */
  ep=1.e-6;  /* relative tolerance */
  tr=1.e-3;  /* absolute tolerance */
  h=1.e-7;

  void* dm = (void *) rhs;
  void* jm = (void *) jac;

  /*************************** ode params *****************************/

  for ( int i=0; i<steps.size(); i++ ) {

    double dt=steps[i].dt, vm=steps[i].vm;
    Model::transition_matrix(m, vm, Q);

    for(int j=0; j<128; j++) ipar[j] = 0;

    if ( steps[i].stype == NONE ) {

      t=0.0; t_end=dt;
      dodesol(ipar, &N, &t, &t_end, y,
        dm, jm, &h, &hm, &ep, &tr, dpar, kd, &ierr);

      if ( ierr ) {
        free(Q); free(dpar); free(kd); return -1;
      }
    }

    else {

      int n_steps = ceil(dt / steps[i].stepsize);
      double *c_mat = (steps[i].dtype == CONDUCTANCE) ? m.C : m.F;

      double *vals = (double*) malloc (n_steps*sizeof(double));
      vals[0] = cblas_ddot(N, y, 1, c_mat, 1);

      for ( int j=1; j<n_steps; j++ ) {
        t=0.0; t_end=steps[i].stepsize;

        dodesol(ipar, &N, &t, &t_end, y,
          dm, jm, &h, &hm, &ep, &tr, dpar, kd, &ierr);

        if ( ierr ) {
          free(Q); free(dpar); free(kd); return -1;
        }

        vals[j] = cblas_ddot(N, y, 1, c_mat, 1);
      }

      if ( steps[i].stype == TAU ) {
        double arg1=steps[i].args[0], arg2=steps[i].args[1];
        double res = tau(n_steps, vals, arg1, arg2);
        out.push_back(res*steps[i].stepsize);
      }

      else if ( steps[i].stype == PEAK ) {
        double res = peak(n_steps, vals);
        out.push_back(res);
      }

      else if ( steps[i].stype == TRACE ) {
        for ( int j=0; j<n_steps; j++ ) {
          out.push_back(vals[j]);
        }
      }

      free(vals);
    }
  }
  free(dpar); free(kd); return 1;
}

int step_exp(Model::Model& m, vector<Step>& steps,
  double *y, vector<double>& out) {/* TODO */};


double cost(Model::Model& m, ChannelProtocol& proto,
  SolverParameter sparam)
{
  std::vector<double> output;
  int ierr, N = m.n_states();

  double *y0 = (double*) malloc(N*sizeof(double));
  double *y  = (double*) malloc((N+N*N)*sizeof(double));
  Model::initial_state(m, proto.v0, y0);

  for ( int i=0; i<proto.n_traces; i++ ) {
    memcpy(y, y0, N*sizeof(double));

    if ( sparam.simulation_mode() == SolverParameter::ODE ) {
      ierr = step_ode(m, proto.traces[i], y, output);
      if ( ierr ) {free(y); free(y0); return 1e6; }
    }
    else {
      ierr= step_exp(m, proto.traces[i], y, output);
      if ( ierr ) {free(y); free(y0); return 1e6; }
    }
  }

  int idx, n = output.size();

  if ( proto.params.normalize() ) {
    idx = cblas_idamax(n, &output[0], 1);
    cblas_dscal(n, 1/output[idx], &output[0], 1);
  }

  vdSub(n, &output[0], &proto.data[0], &output[0]);
  double err = (1.0/n) * cblas_ddot(n, &output[0], 1, &output[0], 1);
  free(y0); free(y); return err;
}


double model_penality(Model::Model& m, SolverParameter& sparam)
{
  double penality = 0;
  penality += sparam.node_penality() * m.n_states();
  penality += sparam.edge_penality() * m.n_edges();
  return penality;
}


double cost(Model::Model& m, vector<ChannelProtocol>& protos,
  SolverParameter& sparam)
{
  double error = 0;
  for ( int i = 0; i < protos.size(); i++ ) {
    error += protos[i].params.weight() * cost(m, protos[i], sparam);
  }
  return  error + model_penality(m, sparam);
}
