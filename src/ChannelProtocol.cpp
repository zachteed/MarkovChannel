#include "ChannelProtocol.hpp"
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


ChannelProtocol::ChannelProtocol(std::string& prototxt)
{
  int fd = open(prototxt.c_str(), O_RDONLY);
  FileInputStream* input = new FileInputStream(fd);
  google::protobuf::TextFormat::Parse(input, &params);

  std::ifstream datfile;
  datfile.open(params.source().c_str());

  double x; std:string line;
  data = dmatrix_t();

  while (std::getline(datfile, line)) {
    data.push_back(std::vector<double>());
    std::vector<double>& tmp = data.back();
    std::istringstream iss(line);
    while (iss>>x) {tmp.push_back(x);}
  }
  delete input;
}

std::ostream& operator<<(std::ostream& os, const ChannelProtocol& proto)
{
  os << proto.params.name() << "\n";
  for (int i = 0; i < proto.data.size(); i++) {
    for (int j = 0; j < proto.data[i].size(); j++) {
      os << proto.data[i][j] << "\t";
    }
    os << "\n";
  }
  return os << std::endl;
}
