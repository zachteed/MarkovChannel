#include "../include/cost.hpp"
#include "intel_ode.h"

using namespace MarkovChannel::SimulationParameter


int void rhs(int *n, double *t, double *y, double *f);
inline void jac_mat(int *n, double *t, double *y, double *a);


inline double peak(int n, double* x, double& peak)
{
  int idx; Math::idamax(n, x, &idx);
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


int step(Model& model, double* y0, double* out, double vm, double dt, double T, double step_size,
  MarkovChannel::ProtocolStep& step, MarkovChannel::SolverParameter& solver_param)
{

  size_t N=model.n_states(), E=model.n_edges();
  double* Q = (double*) malloc(N*N*sizeof(double));
  Markov::transition_matrix(model, vm, T, Q);

  if (solver_param.simulation_mode == SimulationMode::EXPM) {

    double* expm = (double*) malloc(N*N*sizeof(double));
    double* y = (double*) malloc(N*sizeof(double));
    memcpy(y, y0, N*sizeof(double));

    if (step.step_type == MarkovChannel::ProtocolStep::StepType::NONE) {
      Math::expm(transition, dt, expm);

      cblas_dgemv(CblasColMajor, CblasNoTrans,
          N, N, 1.0, expm, N, y, 1, 0.0, y0, 1);

      free(Q); free(y); free(expm); return 1;
    }

    else {

      Math::expm(transition, stepsize, exmp);
      int n_steps = int(ceil(dt / stepsize));

      double* yt = (double*) malloc(n_steps*sizeof(double));
      double* dtype = step.conductance ? model.G : model.F;
      *(yt++) = cblas_ddot(N, y, 1, dtype, 1);

      for (int t=1; t<n_steps; t++) {

        cblas_dgemv(CblasColMajor, CblasNoTrans,
            N, N, 1.0, expm, N, y0, 1, 0.0, y, 1);

        *(yt++) = cblas_ddot(N, y, 1, dtype, 1);

        cblas_dswap(N, y0, 1, y, 1);
      }

    }
  }

  else {

    double* Q = (double*) malloc(N*N*sizeof(double));
    Markov::transition_matrix(model, vm, T, Q);
    mkl_dimatcopy('R', 'T', N, N, Q, 1.0, N, N);

    int *ipar = (int*) calloc(128, sizeof(int));
    int *kd = (int*) malloc(N*sizeof(int));
    double *dpar = (double*) malloc(sizeof(double)*max(13*N,(7+2*N)*N));

    int n, ierr, i;
    double t, t_end, h, hm, ep, tr;

    /*************************** ode params *****************************/

    hm=1.e-12; /* minimal step size for the methods */

    ep=1.e-6;  /* relative tolerance. The code cannot guarantee
           the requested accuracy for ep<1.d-9 */

    tr=1.e-3;  /* absolute tolerance */

    t=0e0;
    h=1.e-7;

    /*************************** ode params *****************************/

    int void rhs(int *n, double *t, double *y, double *f) {
      cblas_dgemv(CblasColMajor, CblasNoTrans, N, N, 1.0, Q, N, y, 1, 0.0, f, 1);
    }
    inline void jacmat(int *n, double *t, double *y, double *a, double *Q) {
      memcpy(a, Q, N*N*sizeof(double));
    }

    if (step.step_type == MarkovChannel::ProtocolStep::StepType::NONE) {
      t_end = dt;
      dodesol(ipar,&n,&t,&t_end,y,rhs,jacmat,&h,&hm,&ep,&tr,dpar,kd,&ierr);
    }

    else {
      int n_steps = int(ceil(dt / stepsize));
      double* yt = (double*) malloc(N*n_steps*sizeof(double));

      for (int i=0; i<n_steps; i++) {
        t=0; t_end=step_size;
        dodesol(ipar,&n,&t,&t_end,y,rhs,jacmat,&h,&hm,&ep,&tr,dpar,kd,&ierr);
        memcpy(yt[i*N], y, N*sizeof(double));
      }
    }
    free(Q); free(ipar); free(kd); free(dpar);
  }

  double output;

  switch (step) {

    case ProtocolStep::StepType::Peak:
    find_peak(N, yt, output);
    *out = output; break;

  case ProtocolStep::StepType::Tau:
    find_tau(N, yt, out, step.val1, step.val2, output);
    *out = output; break;

  case ProtocolSte::StepType::Trace:
    std::memcpy(out, yt, n_steps*sizeof(double); break;

  default:
    break;

  }
}


double cost(Model& model, MarkovChannel::ChannelProtocol& proto,
  MarkovChannel::SolverParameter& solver_param)
{
  MarkovChannel::SolverParamer::SimulationMode sim_param = solver_param.SimulationMode;
  size_t index=0, splits=1, nxt_splits=1, N=model.n_states, E=model.n_edges;

  double* state = (double*) malloc(N*sizeof(double));
  Markov::init_state(model, state, proto.v0, proto.temp);

  for (size_t i = 0; i < proto.steps.size(); i++) {
    MarkovChannel::ProtocolStep& step = proto.steps[i];
    double* dt_vec=NULL, vm_vec=NULL;

    if (!step.has_dt()) {
      nxt_splits = proto.data[index].size()
      dt_vec = &proto.data[index][0];
      index = index + 1;
    }

    if (!step.has_vm()) {
      nxt_splits = proto.data[index].size()
      vm_vec = &proto.data[index][0];
      index = index + 1;
    }

    if (splits != nxt_splits) {
      size_t idx=0, cnt=model.n_states*sizeof(double);
      double* nxt_state = (double*) malloc(nxt*splits*count);
      for (int j=0; j<nxt_splits; j++) {
        std::memcpy(nxt_state[idx], state, count);
        idx += model.n_states;
      }
      std::free(state); state = nxt_state;
      splits = nxt_splits;
    }

    int idx = 0;
    for(int j=0; j<splits; j++) {
      double dt = dt_vec ? dt_vec[j] : step[i].dt;
      double vm = vm_vec ? vm_vec[j] : step[i].vm;

      step(model, state[idx], dt, vm, step[i].step_type, sim_param);
    }
  }
}


double model_penality(Model& model, MarkovChannel::SolverParameter& solver_param)
{
}


double cost(Model& model, std::vector<MarkovChannel::ChannelProtocol>& protocols,
   MarkovChannel::SolverParameter& solver_param)
{
  size_t n_protocols = protocols.size()
  MarkovModel markov = Model::initialize(model);

  double* err = (double *) malloc(n_protocols*sizeof(double));
  double* wts = (double *) malloc(n_protocols*sizeof(double));
  for (size_t i = 0; i < n_protocols; i++) {
    err[i] = cost(markov, protocols[i], solver_param);
	wts[i] = protocols[i].weight;
  }

  double error =


}
