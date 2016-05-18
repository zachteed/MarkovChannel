#include "cost.hpp"
#include "intel_ode.h"
#include <math.h>

using namespace MarkovChannel::SimulationParameter;
using namespace MarkovChannel::ProtocolStep;

void print_matrix(double* A, int n, int m) {
  for (int i=0; i<n; i++) {
    for (int j=0; j<m; j++) {
      printf("%8.4f\t", A[i*m+j]);
    }
    printf("\n");
  }
}


inline double peak(int n, double* x, double& peak)
{
  int idx; cblas_idamax(n, x, &idx);
  peak = x[idx]; return 1;
}

inline double tau(int n, double* x, double v1, double v2, double& tau)
{
  int idx; Math::idamax(n, x, &idx);
  Math::scal(N, 1/x[idx], x);

  if (v2 > v1) {
    int i1=0; i2=0; for(int t=0, t<idx; t++) {
      if (!i1 && x[t+1] > v1 && x[t] <= v1) i1=t;
      if (!i2 && x[t+1] > v2 && x[t] <= v2) it=t;
    }
    tau = i2 - i1; return (i1&&i2) ? 1 : -1;
  } else {
    int i1=0; i2=0; for(int t=idx, t<n-1; t++) {
      if (!i1 && x[t+1] < v1 && x[t] >= v1) i1=t;
      if (!i2 && x[t+1] < v2 && x[t] >= v2) it=t;
    }
    tau = i2 - i1; return (i1&&i2) ? 1 : -1;
  }
}


// int exp_step() {
//   if (solver_param.simulation_mode == SimulationMode::EXPM) {
//
//     double* expm = (double*) malloc(N*N*sizeof(double));
//     double* y = (double*) malloc(N*sizeof(double));
//     memcpy(y, y0, N*sizeof(double));
//
//     if (step.step_type == MarkovChannel::ProtocolStep::StepType::NONE) {
//       Math::expm(transition, dt, expm);
//
//       cblas_dgemv(CblasColMajor, CblasNoTrans,
//           N, N, 1.0, expm, N, y, 1, 0.0, y0, 1);
//
//       free(Q); free(y); free(expm); return 1;
//     }
//
//     else {
//
//       Math::expm(transition, stepsize, exmp);
//       int n_steps = int(ceil(dt / stepsize));
//
//       double* yt = (double*) malloc(n_steps*sizeof(double));
//       double* dtype = step.conductance ? model.G : model.F;
//       *(yt++) = cblas_ddot(N, y, 1, dtype, 1);
//
//       for (int t=1; t<n_steps; t++) {
//
//         cblas_dgemv(CblasColMajor, CblasNoTrans,
//             N, N, 1.0, expm, N, y0, 1, 0.0, y, 1);
//
//         *(yt++) = cblas_ddot(N, y, 1, dtype, 1);
//
//         cblas_dswap(N, y0, 1, y, 1);
//       }
//
//     }
//   }
// }



double step_ode(int N, double *Q, double *y, double dt, int n_steps, double stepsize, 
  MarkovChannel::ProtocolStep& step_params, double *vals, double *args, double& out);
{

  struct ode_struct {
    double *q_mat;
    inline void rhs(int *n, double *t, double *y, double *f) {
      cblas_dgemv(CblasColMajor, CblasNoTrans,
        N, N, 1.0, Q, N, y, 1, 0.0, f, 1);
    }
    inline void jac(int *n, double *t, double *y, double *a) {
      mkl_domatcopy('C', 'T', N, N, 0.0, Q, N, a, N);
    }
  }

  int ipar[128]; for(int i=0; i<128; i++) ipar[i] = 0;
  int *kd = (int*) malloc(N*sizeof(int));

  int n, ierr, i, j, cnt = max(13*N, (7+2*N)*N)*sizeof(double);
  double *dpar = (double*) malloc(cnt), t, t_end, h, hm, ep, tr;

  /*************************** ode params *****************************/

  hm=1.e-12; /* minimal step size for the methods */
  ep=1.e-6;  /* relative tolerance */
  tr=1.e-3;  /* absolute tolerance */
  h=1.e-7;

  /*************************** ode params *****************************/

  ode_struct ode; ode.q_mat = Q;
  void (*rhs) (int*, double*, double*, double*) = ode.rhs;
  void (*jac) (int*, double*, double*, double*) = ode.jac;

  if (step_params.stype() == StepType::NONE) {
    t=0.0; t_end=dt;
    dodesol(ipar, &n, &t, &t_end, y, rhs, jac, &h, &hm, &ep, &tr, dpar, kd, &ierr);
    free(kd); free(dpar); return ierr ? -1 : 1;
  }

  vals[0] = cblas_ddot(N, args, 1, y, 1);
  for ( j=1; j < n_steps; j++ ) {
    t=0.0; t_end = stepsize;
    dodesol(ipar, &n, &t, &t_end, y, rhs, jac, &h, &hm, &ep, &tr, dpar, kd, &ierr);
    vals[j] = cblas_ddot(N, args, 1, y, 1); 
    if(ierr) {free(kd); free(dpar); return -1;}
  }
  free(kd); free(dpar);

  double prm1 = step_params.prm1();
  double prm2 = step_params.prm2();

  switch(step_params.stype()) {
   
    case StepType::PEAK:
      out = peak(n_steps, vals); break;

    case StepType::TAU:
      out = tau(n_steps, vals, prm1, prm2); break;

    default:
      break;

  }

  return 1;
}


double cost(Model& model, MarkovChannel::ChannelProtocol& proto,
  MarkovChannel::SolverParameter& solver_param)
{
  SimulationMode sim_param = solver_param.SimulationMode();
  size_t index=0, splits=1, nxt_splits=1, ierr;
  int N=model.n_states(), E=model.n_edges();

  double* y = (double*) malloc(N*sizeof(double));
  double* Q = (double*) malloc(N*N*sizeof(double));
  Model::init_state(model, proto.params.v0(), y);

  for (size_t i = 0; i < proto.step_size(); i++) {

    MarkovChannel::ProtocolStep& step = proto.params.step(i);
    double *vals = NULL, out;

    if ( !step.has_dt() || !step.has_vm() ) {
       nxt_splits = proto.y[i].size();
    }

    if ( splits != nxt_splits ) {
      double *y_next = (double*) malloc(N*nxt_splits*sizeof(double));
      for ( int j=0; j<nxt_splits; j++ ) {
        memcpy(y_next[j*N], y, N*sizeof(double));
      }
      free(y); y = y_next; splits = nxt_splits;
    }

    if (proto.stype == StepType::NONE) {
      for ( int j=0; j<splits; j++ ) {
        double dt = step.has_dt() ? step.dt() : proto.y[i][j];
        double vm = step.has_vm() ? step.vm() : proto.y[i][j];
        double stepsize = step.stepsize(); int n_steps=1;
        Model::transition_matrix(model, vm, Q);
        ierr = step_ode(N, Q, &y[j*N], dt, n_steps, stepsize, step, NULL, NULL, out);
    }
    else {
       vals = (double*) malloc(splits*sizeof(double));
       for ( int j=0; j<splits; j++ ) {
         double dt = step.has_dt() ? step.dt() : proto.y[i][j];
         double vm = step.has_vm() ? step.vm() : proto.y[i][j];
         double stepszie = step.stepsize(); int n_steps = ceil(dt/stepsize);
         Model::transition_matrix(model, vm, Q);
         ierr = step_ode(N, Q, &y[j*N], dt, n_steps, stepsize



      std::vector<double> output;

      for ( int j=0; j<splits; j++ ) {
        double dt=proto.x[j], vm=step.vm(), step_sz=dt;
        if ( step.step_type() != StepType::NONE ) {
          step_sz = step.step_sz();
        }
        Markov::transition_matrix(model, vm, Q);
        ode_step(N, Q, y[j*N], dt, step_sz, yt);

        switch ( step.step_type() ) {

          case StepType::PEAK :
            output.push_back(find_peak(yt, step_sz)); break;

          case StepType::TAU :
            output.push_back(find_tau(yt, step_sz)); break;

          case default:
            break;

        }

        if ( step.step_type() != StepType::None ) {

        }
      }
    }




double model_penality(Model& model, MarkovChannel::SolverParameter& solver_param)
{
  double penality = 0;
  penality += solver_mode.node_penality() * model.n_states();
  penality += solver_mode.edges_penality() * model.n_edges();

  return penality;
}


double cost(Model& model, std::vector<MarkovChannel::ChannelProtocol>& protocols,
   MarkovChannel::SolverParameter& solver_param)
{
  size_t n_protocols = protocols.size()
  MarkovModel markov = Model::initialize(model);

  double error = 0;

  for (size_t i = 0; i < n_protocols; i++) {
    err[i] = protocols[i].params().weight() *
       cost(markov, protocols[i], solver_param);
  }

  error += model_penality(model, solver_param);
  return error;
}
