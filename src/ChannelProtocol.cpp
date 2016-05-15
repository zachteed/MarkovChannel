#include "ChannelProtocol.hpp"



protocol::protocol(std::string& prototxt)
{
  std::fstream input(prototxt, ios::in | ios::binary);
  params = MarkovChannel::ChannelProtocol.ParseFromIstream(&input);

  std::ifstream datfile;
  datfile.open(proto.source);

  double x; std:string line;
  data = dmatrix_t();

  while (std::getline(datfile, line)) {
    data.push_back(std::vector<double>());
    std::vector<double>& tmp = data.back();
    std::istringstream iss(line);
    while (iss>>x) {tmp.push_back(x);}
  }
}

std::ostream& operator<<(std::ostream& os, const protocol& proto)
{
  os << proto.params.name << "\n";
  for (int i = 0; i < proto.data.size(); i++) {
    for (int j = 0; j < proto.data[i].size(); j++) {
      os << proto.data[i][j] << "\t";
    }
    os << "\n";
  }
  return os << std::endl;
}
