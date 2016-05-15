#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "ChannelProtocol.hpp"

using namespace std;


int main(int argc, char* argv[])
{
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " CHANNEL_PROTOCOL " << endl;
    return -1;
  }

  protocol proto(argv[1]);
  cout << proto << endl;
  //
  // std::cout << argv[1] << std::endl;
  return 0;
}
