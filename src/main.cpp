#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

#include "MarkovChannel.pb.h"

typedef std::vector<std::vector<double> > dmatrix_t;

struct ChannelProtocol {
  MarkovChannel::ChannelProtocol params;
  dmatrix_t data;
}


void init_proto(ChannelProtocol& proto)
{
  std::ifstream datfile;
  datfile.open(proto.source);
  
  double x; std:string line;
  dmatrix_t& data = proto.data;

  while (std::getline(datfile, line)) {
    data.push_back(std::vector<double>);
    std::vector<double>& tmp = data.back();
    std::istringstream iss(line);
    while (iss>>x) {tmp.push_back(x);}
  }
} 


int main(int argc, char* argv[])
{
  GOOGLE_PROTOBUF_VERIFY_VERSION;
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " PROTOCOL_LIST" << std::endl;
    return -1;
  }

  std::vector<ChannelProtocol> protocols;
  MarkovChannel::ChannelProtocol prms;
  std:ifstream protolist;
  protolist.open(argv[1]);
  
  std:string prototxt;
  while(protlist >> prototxt) {
    std:fstream input(prototxt, ios::in | ios::binary);
    prms = MarkovChannel::ChannelProtocol.ParseFromIstream(&input));
    protocols.push_back(ChannelProtocol());
    ChannelProtocol& proto = protocols.back();
    proto.params = prms; init_proto(proto);
  }
}


