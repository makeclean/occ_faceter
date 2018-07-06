#include <map>
#include "moab/Core.hpp"

typedef merge_pairs_t std::pair<moab::EntityHandle,moab::EntityHandle>;

class DAGMCTopology {
  public:

    DAGMCTopology(moab::Core *MBI = NULL);
    ~DAGMCTopology();

    // load the file
    moab::ErrorCode load_file(const std::string filename);

    // setup tags
    moab::ErrorCode setup_tags();

    // identify coincident surfaces
    moab::ErrorCode identify_coincident(std::vector<merge_pairs_t> &coincident);

    // merge coincident surfaces
    moab::ErrorCode merge_coincident(const std::vector<merge_pairs_t> coincident, bool delete = true);

  private:

  private:
    moab::Core *mbi; // moab pointer
    bool own_moab; // does the class own its own moab instance?
    moab::EntityHandle input_set; // the meshset corresponding to input data
    // tag data 
    moab::Tag geometry_dimension_tag, id_tag;
    moab::Tag faceting_tol_tag, geometry_resabs_tag;
    moab::Tag category_tag;

}
