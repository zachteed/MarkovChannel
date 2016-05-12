#include "../include/cost.hpp"



void step(Model& model, MarkovChannel::)


double cost(Model& model, MarkovChannel::ChannelProtocol& proto, 
  MarkovChannel::SolverParameter& solver_param)
{
  for ()  
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
