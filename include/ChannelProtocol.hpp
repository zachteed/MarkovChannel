#ifndef CHANNEL_PROTOCOL_HPP_
#define CHANNEL_PROTOCOL_HPP_

#include "MarkovChannel.pb.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>


typedef std::vector<std::vector<double> > dmatrix_t;

struct ChannelProtocol {
  ChannelProtocol(std::string& prototxt);
  MarkovChannel::ProtocolParameter params;
  dmatrix_t data;
};

std::ostream& operator<<(std::ostream& os, const ChannelProtocol& proto);

#endif
