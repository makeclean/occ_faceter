#ifndef READ_METADATA_HH
#define READ_METADATA_HH 1

#include <string>
#include <unordered_map>

typedef std::unordered_map<uint64_t, std::string> MaterialsMap;

void read_metadata(std::string json_file, MaterialsMap &mat_map);

#endif // READ_METADATA_HH
