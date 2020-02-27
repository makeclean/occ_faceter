#include "read_metadata.hh"

#include <fstream>
#include <iostream>

#include "UniqueId/json.hpp"

using json = nlohmann::json;

void read_metadata(std::string json_file, MaterialsMap &mat_map) {
  std::ifstream json_stream(json_file);
  if (json_stream.fail()) {
    std::cerr << "Warning: Failed to read file " + json_file << std::endl;
  } else {
    json j = json::parse(json_stream);

    for (const auto &p : j) {
      uint64_t uniqueID = p["uniqueID"].get<uint64_t>();
      std::string material = p["material"].get<std::string>();
      mat_map[uniqueID] = material;
    }
  }
}
