#include "Model.hpp"
#include "cost.hpp"

#include "math_functions.hpp"
#include "graph_functions.hpp"
#include "ChannelProtocol.hpp"
#include "SimulatedAnnealing.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdio.h>
#include <math.h>

#include <fcntl.h>
#include <google/protobuf/io/coded_stream.h>
#include <google/protobuf/io/zero_copy_stream_impl.h>
#include <google/protobuf/text_format.h>

using namespace std;
using google::protobuf::io::FileInputStream;
using google::protobuf::io::FileOutputStream;
using google::protobuf::io::ZeroCopyInputStream;
using google::protobuf::io::CodedInputStream;
using google::protobuf::io::ZeroCopyOutputStream;
using google::protobuf::io::CodedOutputStream;
using google::protobuf::Message;


int print_ic (Graph::Graph& G)
{
  double* ic = (double*) calloc(2*G.E*G.N, sizeof(double));
  for(int i=0; i<G.E; i++) {
    Graph::Edge& e = G.edges[i];
    ic[2*G.E*e.V1 + 2*i] = -1;
    ic[2*G.E*e.V2 + 2*i] = 1;
    ic[2*G.E*e.V1 + 2*i+1] = 1;
    ic[2*G.E*e.V2 + 2*i+1] = -1;
  }
  for (int i=0; i<G.N; i++) {
    for (int j=0; j<2*G.E; j++) {
      printf("%8.4f\t", ic[i*2*G.E+j]);
    }
    printf("\n");
  }
  free(ic);
}

MarkovChannel::SolverParameter solver_param;
vector<ChannelProtocol> protos;


double cost_f(Model::Model* m, bool print)
{
  return cost(*m, protos, solver_param, print);
}



int main(int argc, char* argv[])
{
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  Math::init_stream(time(NULL));

  int fd = open(argv[1], O_RDONLY);
  FileInputStream* input = new FileInputStream(fd);

  google::protobuf::TextFormat::Parse(input, &solver_param);
  Model::prms = solver_param.model_param(); delete input;


  string s1("data/protocols/gv.prototxt");
  string s2("data/protocols/inac.prototxt");

  protos.push_back(ChannelProtocol(s1));
  protos.push_back(ChannelProtocol(s2));

  MarkovChannel::SAParameter sa_param = solver_param.sa_param();
  SimulatedAnnealing::solve(cost_f, sa_param);


  google::protobuf::ShutdownProtobufLibrary();
  exit(0);




  // Model::Model *m = new Model::Model(5);
  // Model::Model *n = Model::neighbor(m);
  // Model::Model *p = new Model::Model(m);
  //
  // cout << *m << endl;
  // cout << *n << endl;
  // cout << *p << endl;

  // const int k_max = 10000;
  // const double T0 = 0.002;
  // const double gamma = 0.5;
  //
  // double T = T0;
  //
  // Model::Model *m = new Model::Model(6);
  // double f = cost_f(*m);
  //
  // Model::Model *argmin = new Model::Model(m);
  // double f_min = f;
  //
  // for ( int k=0; k<k_max; k++ ) {
  //   Model::Model *n = Model::neighbor(m);
  //
  //   if (n->n_states() != 6 || m->n_states() != 6) {
  //     cout << *n << endl;
  //     cout << *m << endl;
  //   }
  //   double fn = cost_f(*n);
  //   if ( Math::rng_uniform() < exp(-(fn-f)/T) ) {
  //     delete m; f = fn; m = n;
  //   }
  //   else {
  //     delete n;
  //   }
  //   if ( fn < f_min ) {
  //     delete argmin; f_min = f;
  //     argmin = new Model::Model(m);
  //   }
  //
  //   if ( Math::rng_uniform() < 0.001 ) {
  //     delete m; f = f_min;
  //     m = new Model::Model(argmin);
  //   }
  //
  //   if ( k  % 100 == 0 ) {
  //     T *= gamma;
  //     std::cout << f_min << std::endl;
  //     cost_f(*m, true);
  //     std::cout << std::endl;
  //     // std::cout << f_min << std::endl;
  //     // std::cout << *m << std::endl;
  //
  //   }
  // }
  //
  // delete m;

  // delete m; delete n;


  // for ( int i=0; i < 100; i++ ) {
  //   Model::Model m(6);
  //   cout << m << endl;
  //   cout << cost_f(m) << endl;
  // }
  //
  // exit(0);
  //
  //
  //
  // // Model::Model model(5);
  // // cost(model, protos, solver_param);
  // // cout << model << endl;
  //
  // MarkovChannel::SAParameter sa_prm = solver_param.sa_param();
  //
  // SimulatedAnnealing::solve(cost_f, sa_prm);


  // for (int i=0; i < 10000; i++) {
  //   Model::Model model(8, 0.4);
  //   cost(model, protos, solver_param);
  // }

}
