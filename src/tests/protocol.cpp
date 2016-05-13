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
    std::cerr << "Usage: " << argv[0] << " CHANNEL_PROTOCOL " << std::endl;
    return -1;
  }

  std::fstream input(argv[0], ios::in | ios::binary);
  MarkovChannel::ChannelProtocol cp =
    MarkovChannel::ChannelProtocol.ParseFromIstream(&input));

}
