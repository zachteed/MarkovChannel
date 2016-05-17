#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <vector>
#include <string>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <mkl.h>

#include "math_functions.hpp"
#include "graph_functions.hpp"
#include "MarkovChannel.pb.h"

namespace Model {

  MarkovChannel::ModelParameter prms;

  struct Model {

    static int count;
    int id;

    Model(int N, double p=0.2);
    Model();
    ~Model();

    Graph::Graph G;

    double *rs, *rk;
    double *C, *F;
    double *r_vec;

  };

  Model* neighbor(Model& m, int num=1);

  double* initial_state(Model& m, double vm, double* s);

  double* transition_matrix(Model& m, double vm, double* Q);

  std::ostream& operator<<(std::ostream& os, const Model& m);

}

#endif
