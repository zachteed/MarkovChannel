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

void load_protocols(std::string& protolst)
{
  std::ifstream protofile;
  protofile.open(protolst.c_str());

  std::string line;
  while ( protofile >> line ) {
    protos.push_back(ChannelProtocol(line));
  }
}



int main(int argc, char* argv[])
{

  GOOGLE_PROTOBUF_VERIFY_VERSION;
  Math::init_stream(time(NULL));

  int fd = open(argv[1], O_RDONLY);
  FileInputStream* input = new FileInputStream(fd);

  google::protobuf::TextFormat::Parse(input, &solver_param);
  Model::prms = solver_param.model_param(); delete input;

  std::string proto_list = solver_param.protocol_list();
  cout << proto_list << endl;
  load_protocols(proto_list);


  MarkovChannel::SAParameter sa_param = solver_param.sa_param();
  SimulatedAnnealing::solve(cost_f, sa_param);


  google::protobuf::ShutdownProtobufLibrary();
  exit(0);


}
