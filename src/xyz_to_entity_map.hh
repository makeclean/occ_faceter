#ifndef XYZ_TO_ENTITY_MAP_HH
#define XYZ_TO_ENTITY_MAP_HH

#include <array>
#include <unordered_map>

#include "moab/EntityHandle.hpp"

struct xyz_coords {
  std::array<double, 3> coords;

  xyz_coords(const std::array<double, 3> &values) : coords(values) {}
};

inline bool operator==(const xyz_coords& lhs, const xyz_coords& rhs) {
    return lhs.coords == rhs.coords;
}

// custom specialization of std::hash can be injected in namespace std
template<>
struct std::hash<xyz_coords>
{
    std::size_t operator()(xyz_coords const& c) const noexcept
    {
        std::size_t h1 = std::hash<double>{}(c.coords[0]);
        std::size_t h2 = std::hash<double>{}(c.coords[1]);
        std::size_t h3 = std::hash<double>{}(c.coords[2]);
        return h1 ^ (h2 << 1) ^ (h3 << 2); // or use boost::hash_combine
    }
};

typedef std::unordered_map<xyz_coords, moab::EntityHandle> coordinates_to_entity_map;

#endif // XYZ_TO_ENTITY_MAP_HH
