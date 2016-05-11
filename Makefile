SRC_DIR := src
DST_DIR := src

protoc -I=$(SRC_DIR) --cpp_out=$(DST_DIR) $(SRC_DIR)/MarkovChannel.proto
