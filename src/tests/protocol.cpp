#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "MarkovChannel.pb.h"
#include "ChannelProtocol.hpp"


using namespace std;

int main(int argc, char* argv[])
{

  GOOGLE_PROTOBUF_VERIFY_VERSION;

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " CHANNEL_PROTOCOL " << endl;
    return -1;
  }

  string protobuf(argv[1]);
  ChannelProtocol proto(protobuf);
  cout << proto << endl;

  google::protobuf::ShutdownProtobufLibrary();

  return 0;
}
