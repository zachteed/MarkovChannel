#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <string>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cblas.h>
#include <iostream>

#include "math_functions.hpp"


namespace Markov {

  struct Edge {
    int V1, V2;
  };

  struct Model {
    std::vector<Edge> edges;
    double *rs, *rk, *G, *F;
    double *rates;
    int n_states, n_edges, n_prms;
  };

  int initial_state(Model& m, double vm, double* s);

  int transition_matrix(Model& m, double vm, double* Q);

  void neighbors();

}

#endif
