#include "Model.hpp"
#include "math_functions.hpp"
#include "MarkovChannel.pb.h"

#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>

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


int main(int argc, char* argv[]) {

  Math::init_stream(time(NULL));

  int fd = open(argv[1], O_RDONLY);
  FileInputStream* input = new FileInputStream(fd);

  MarkovChannel::SolverParameter solver_param;
  google::protobuf::TextFormat::Parse(input, &solver_param);
  Model::prms = solver_param.model_param(); delete input;


  Model::Model* n;

  Model::Model model(5);
  n = Model::neighbor(model);
  delete n;



  return 1;
}



  // std::vector<Edge> edges;
  // edges.push_back({0, 1});
  // edges.push_back({0, 2});
  // edges.push_back({1, 3});
  // edges.push_back({2, 3});
  //
  // double rs[] = { 0.5377,    1.8339,
  //                -2.2588,    0.8622,
  //                 0.0159,   -0.0654,
  //                -0.0217,    0.0171};
  //
  // double rk[] = { 3.5784,    0.0363,
  //                 2.7694,   -0.0032,
  //                -1.3499,    0.0357,
  //                 3.0349,   -0.0102};
  //
  // Model model;
  // model.edges = edges;
  // model.rs = rs;
  // model.rk = rk;
  // model.n_states = 4;
  // model.n_edges = 4;
  // model.n_prms = 2;
  // model.rates = NULL;
  //
  // double *s0 = (double*) malloc(4*sizeof(double));
  // Markov::initial_state(model, 0, s0);
  //
  // double *Q = (double*) malloc(4*4*sizeof(double));
  // Markov::transition_matrix(model, 0, Q);
  //
  // for (int i = 0; i < 4; i++) {
  //   for (int j = 0; j < 4; j++) {
  //     std::cout << Q[4*i+j] << "\t";
  //   }
  //   std::cout << std::endl;
  // }
