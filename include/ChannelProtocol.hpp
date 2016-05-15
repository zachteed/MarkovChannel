#ifndef CHANNEL_PROTOCOL_HPP_
#define CHANNEL_PROTOCOL_HPP_

#include "MarkovChannel.pb.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>


typedef std::vector<std::vector<double> > dmatrix_t;

struct protocol {
  protocol(std::string& prototxt);
  MarkovChannel::ChannelProtocol params;
  dmatrix_t data;
};

std::ostream& operator<<(std::ostream& os, const protocol& proto);

#endif
