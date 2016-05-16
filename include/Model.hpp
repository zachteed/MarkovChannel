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
    Model(int N, double p=0.2);
    ~Model();

    Graph::Graph G;

    double *rs, *rk;
    double *C, *F;
    double *r_vec;
  };

  int neighbor(Model& m, int n=1);

  int mutate(Model& m, int n=1);

  double* initial_state(Model& m, double vm, double* s);

  double* transition_matrix(Model& m, double vm, double* Q);

}

#endif
