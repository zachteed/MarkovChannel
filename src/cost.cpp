#include "../include/cost.hpp"


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


int step(Model& model, double* s0, double* out, double vm, double dt, double T, double step_size,
  MarkovChannel::ProtocolStep& step, MarkovChannel::SolverParameter& solver_param)
{
  size_t N=model.n_states, E=model.n_edges;
  if (solver_param.simulation_mode == MarkovChannel::SimulationParameter::SimulationMode::EXPM) {
    double* transition = (double*) std::malloc(N*N*sizeof(double));
    double* expm = (double*) std::malloc(N*N*sizeof(double));
    double* state = (double*) std::malloc(N*sizeof(double));
    Markov::transition_matrix(model, vm, T, transition);
    
    if (step.step_type == MarkovChannel::ProtocolStep::StepType::NONE) {
      Math::expm(transition, dt, expm);
      Math::gemv(CblasNoTrans, N, N, 1.0, expm, 0.0, s0, state);
      std::memcpy(s0, state); return 1;
    }
    
    else {
      Math::expm(transition, stepsize, exmp);
      int n_steps = int(ceil(dt / stepsize));
      
      double* vals = (double*) std::malloc(n_steps*sizeof(double));
      double* dtype = step.conductance ? model.G : model.F;
      
      vals[0] = Math::dot(N, s0, dtype);
      for (int t=1; t<n_steps; t++) {
        Math::gemv(CblasNoTrans, N, N, 1.0, expm, 0.0, s0, state);
        vals[t] = Math::dot(N, s0, dtype);
        double* temp = s0; s0 = state; state = temp;
      }
    }
  } 
  
  else {
    double* transition = (double*) malloc(N*N*sizeof(double));
    Markov::transition_matrix(model, vm, T, transition);
    void *(dxdt)(const double*, double*, const double)
    void dxdt(const double* x, double* dxdt, const double /* t */)
    
  }
  
  double output; 
  switch (step) {
    case ProtocolStep::StepType::Peak:
    find_peak(N, vals, output);
    *out = output; break;
  case ProtocolStep::StepType::Tau:
    find_tau(N, vals, out, step.val1, step.val2, output); 
    *out = output; break;
  case ProtocolS::StepType::Trace:
    std::memcpy(out, vals, n_steps*sizeof(double); break;
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
