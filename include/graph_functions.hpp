#ifndef GRAPH_FUNCTIONS_HPP_
#define GRAPH_FUNCTIONS_HPP_

#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>

#include "math_functions.hpp"

namespace Graph {

  struct Edge {
    Edge(int v1, int v2) : V1(v1), V2(v2) {};
    int V1, V2;
  };

  struct Graph {
    Graph() : E(0), N(0) {};
    Graph(int n, double p=0.2) : E(0), N(n) {
      random_connected_graph(n, *this, p);
    }
    std::vector<Edge> edges;
    int E, N;
  }

  int random_connected_graph(int N, Graph& G, double p);

  int adj_list(Graph& G, std::vector<std::vector<int> >& lst);

  int depth_first_search(std::vector<std::vector<int> >& lst,
    int n0, std::vector<int>& nodes);

  int depth_first_search(Graph& G, int n0, std::vector<int>& nodes);

  int connect(Graph& G);

  int add_edge(Graph& G, bool force);

  int add_node(Graph& G);

  int rm_edge(Graph& G, int*& idx, bool reconnect);

  int rm_node(Graph& G, int*& idx, bool reconnect);

  std::ostream& operator<< (std::ostream& os, const Graph& G);

}

#endif
