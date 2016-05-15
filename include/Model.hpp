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
#include "graph_functions.hpp"


namespace Model {

  MarkovChannel::ModelParameter prms;

  struct Model {
    Graph::Graph graph;
    double *rs, *rk;
    double *G, *F, *rates;
  };

  int initial_state(Model& m, double vm, double* s);

  int transition_matrix(Model& m, double vm, double* Q);

  int neighbor(Model& m, int n=1);

  int mutate(Model& m, int n=1);

}

#endif
