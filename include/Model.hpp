#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <stdio>
#include <string>
#include <math>
#include <vector>

struct Edge {
  int V1, V2;
}

struct Model {
  std::vector<Edge> edges;
  double* rs, rk, G, F;
  double* transition_matrix;
  int n_states, n_edges;
}

namespace Model {

  init_state(Model&, double vm);

  transition_matrix(Model&, double vm);

}

#endif
