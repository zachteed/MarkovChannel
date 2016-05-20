#ifndef SIMULATED_ANNEALING_HPP_
#define SIMULATED_ANNEALING_HPP_


#include "MarkovChannel.pb.h"
#include "Model.hpp"
#include "math_functions.hpp"



typedef double (*cost_function)(Model::Model*, bool);

namespace SimulatedAnnealing
{
  void solve(cost_function cost, MarkovChannel::SAParameter& params);
}


#endif
