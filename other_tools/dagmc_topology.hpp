#include <map>
#include "moab/Core.hpp"
#include "moab/GeomTopoTool.hpp"

typedef std::pair<moab::EntityHandle,moab::EntityHandle> merge_pairs_t;

class DAGMCTopology {
  public:
    DAGMCTopology(moab::Core *MBI = NULL);
   ~DAGMCTopology();

    // load the file
    moab::ErrorCode load_file(const std::string filename);

    // setup tags
    moab::ErrorCode setup_tags();

    // do the merge
    moab::ErrorCode perform_merge();

    // save the file
    moab::ErrorCode save_file(const std::string filename);
  
  private:
    // for the curves passed in determine the number of edges in the curve
    moab::ErrorCode get_curve_lengths(moab::Range curves,
				      std::map<int,moab::Range> &curves_and_counts);

     moab::ErrorCode find_curve_pairs(const std::map<int,moab::Range>  curves,
				      std::vector<merge_pairs_t> &pairs);
  
    // identify coincident surfaces
    moab::ErrorCode identify_coincident(std::vector<merge_pairs_t> &coincident);

    // merge coincident surfaces
    moab::ErrorCode merge_coincident_curves(const std::vector<merge_pairs_t> coincident,
				     const bool del = true);

    moab::ErrorCode compare_edges(moab::Range edges1, moab::Range edges2, bool &same);
  
    // given a comparitor curve, and a range of curves to compare against if
    // there are any matches, add them to the pairs vector
    moab::ErrorCode compare_curves(moab::EntityHandle curve1, moab::Range curves,
				   std::vector<merge_pairs_t> &pairs);
  
    // compare two curve meshsets for equivalence
    moab::ErrorCode compare_curve(moab::EntityHandle curve1, moab::EntityHandle curve2,
				   bool &match);
    // remove the duplciate  
    moab::ErrorCode remove_duplicate_curves(const std::vector<merge_pairs_t> &coincident,
					    std::vector<merge_pairs_t> &reduced_list);

    // return the vector of all matching surfaces
    moab::ErrorCode identify_coincident_surfaces(std::vector<merge_pairs_t> &coincident_surfaces);

    // given a surface compare it with all surfaces in the problem and if it matches
    // add it to the pairs vector
    moab::ErrorCode compare_surface_curves(const moab::EntityHandle surface,
					   const moab::Range surface_set,
					   std::vector<merge_pairs_t> &surface_pairs);
  
     // given two surface handles, compare their child curves, if the same set true 
     moab::ErrorCode compare_surfaces(const moab::EntityHandle surface1,
				      const moab::EntityHandle surface2,
				      bool &same);

    // merge duplicate surfaces;
    moab::ErrorCode merge_duplicate_surfaces(const std::vector<merge_pairs_t> duplicates,
					     const bool remove);
					     
  
  private:
    moab::Core *mbi; // moab pointer
    bool own_moab; // does the class own its own moab instance?
    moab::EntityHandle input_set; // the meshset corresponding to input data
    // tag data 
    moab::Tag geometry_dimension_tag, id_tag;
    moab::Tag faceting_tol_tag, geometry_resabs_tag;
    moab::Tag category_tag;

};
