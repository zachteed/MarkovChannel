#ifndef COST_HPP_
#defin COST_HPP_

#include <vector>
#include <iostream>
#include <cstdlib>

#include "math_fucntions.hpp"
#include "graph_functions.hpp"
#include "solver_functions.hpp"

#include "Model.hpp"
#include "MarkovChannel.pb.h"

double cost(Model& model, std::vector<MarkovChannel::ChannelProtocol>& protocols, 
   MarkovChannel::SolverParameter& solver_param);

#endif
