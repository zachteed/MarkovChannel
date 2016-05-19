#include "Model.hpp"
#include "cost.hpp"

#include "math_functions.hpp"
#include "graph_functions.hpp"
#include "ChannelProtocol.hpp"

#include <fstream>
#include <iostream>
#include <vector>
#include <cstdlib>
#include <stdio.h>

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


int main(int argc, char* argv[])
{
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  Math::init_stream(time(NULL));

  int fd = open(argv[1], O_RDONLY);
  FileInputStream* input = new FileInputStream(fd);

  MarkovChannel::SolverParameter solver_param;
  google::protobuf::TextFormat::Parse(input, &solver_param);
  Model::prms = solver_param.model_param(); delete input;


  string protobuf(argv[2]);
  ChannelProtocol proto(protobuf);

  vector<ChannelProtocol> protos;
  protos.push_back(proto);

  // Model::Model model(5);
  // cost(model, protos, solver_param);
  // cout << model << endl;


  for (int i=0; i < 10000; i++) {
    Model::Model model(8, 0.4);
    cost(model, protos, solver_param);
  }

  google::protobuf::ShutdownProtobufLibrary();
  exit(0);
}
