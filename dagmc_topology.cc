#include "dagmc_topology.hpp"
#include "MBTagConventions.hpp"
#include <iostream>

// constructor
DAGMCTopology::DAGMCTopology(moab::Core *MBI) {
  if (MBI == NULL) {
    mbi = new moab::Core();
    own_moab = true;
  } else {
    mbi = MBI;
    own_moab = false;
  }
}

// destructor
DAGMCTopology::~DAGMCTopology() {
  if(own_moab) delete mbi;
}

// load the file
moab::ErrorCode DAGMCTopology::load_file(const std::string filename) {
  moab::ErrorCode rval = moab::MB_FAILURE;
  rval = mbi->create_meshset(moab::MESHSET_SET, input_set);
  MB_CHK_SET_ERR(rval,"failed to create meshset");
  rval = mbi->load_file(filename.c_str(),&input_set);
  MB_CHK_SET_ERR(rval,"failed to load the file");
  
  if(rval != moab::MB_SUCCESS) return rval;

  return setup_tags();
}

// setup the tags needed for later query
moab::ErrorCode DAGMCTopology::setup_tags() {
  moab::ErrorCode rval = moab::MB_FAILURE;
  rval = mbi->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR(rval, "Couldnt get geom dim tag");
  rval = mbi->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR(rval, "Couldnt get id tag");
  rval = mbi->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE, faceting_tol_tag,
			       moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR(rval, "Error creating faceting_tol_tag");
  rval = mbi->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE, 
			     geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR(rval, "Error creating geometry_resabs_tag");
  /*
  rval = mbi->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE, 
			     category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR(rval, "Error creating category_tag");
  */
  return rval;
}

// identify and fix up coincident surfaces
moab::ErrorCode DAGMCTopology::perform_merge() {
  std::vector<merge_pairs_t> coincident;
  moab::ErrorCode rval = identify_coincident(coincident);
  if ( rval != moab::MB_SUCCESS ) return rval;
  if ( coincident.size() == 0 ) std::cout << "no curve pairs found" << std::endl;

  // sort the list into a unified list where the curve identical(?ness)
  // ie. curve 1 same as curve 2, curve 3 same as 2 and curve 4 same as 3
  // i.e. curve1 === curve2 ; curve1 === curve3 ; curve1 === curve4
  std::vector<merge_pairs_t> reduced_list;
  rval = remove_duplicate_curves(coincident,reduced_list);
  
  return merge_coincident(reduced_list,true);
}

// loop through the vector and remove duplicates, and identify 
moab::ErrorCode DAGMCTopology::remove_duplicate_curves(const std::vector<merge_pairs_t> &coincident,
						       std::vector<merge_pairs_t> &reduced_list) {
  moab::ErrorCode rval = moab::MB_FAILURE;
  // map of eh aliases
  std::map<moab::EntityHandle,moab::EntityHandle> alias_list;

  // loop over the pairs
  for ( merge_pairs_t merge_pair : coincident ) {
    // the std::min and max are such that we prefer the lowest valued curve to be the
    // reference curve
    moab::EntityHandle curve1 = std::min(merge_pair.first,merge_pair.second);
    moab::EntityHandle curve2 = std::max(merge_pair.first,merge_pair.second);
    // insert the curve into the list
    if ( alias_list.count(curve1) == 0 )  
      alias_list[curve2] = curve1;
    else // point the second curve to the already existant one
      alias_list[curve1] = alias_list[curve2];
  }

  // turn the map back into a vector for use downstream
  for ( std::pair<moab::EntityHandle,moab::EntityHandle> pairs : alias_list ) {
    merge_pairs_t merge_pairs(pairs.second,pairs.first);
    reduced_list.push_back(merge_pairs);
    //    std::cout << pairs.second << " " << pairs.first << std::endl;
  }
     
  return rval;
}

// identify coincident surfaces
moab::ErrorCode DAGMCTopology::identify_coincident(std::vector<merge_pairs_t> &coincident) {
  // get the list of curves present and compare them for equivlance
  moab::ErrorCode rval = moab::MB_FAILURE;
  moab::Range curve_set;

  // maybe use type and tag
  const int dim = 1;
  const void *tag_value[] = {&dim};
  rval = mbi->get_entities_by_type_and_tag(0,moab::MBENTITYSET,&geometry_dimension_tag,
					   tag_value,1,curve_set);

  MB_CHK_SET_ERR(rval,"couldnt get entities by type and tag");
  // first sort by entity set size - i.e. curves of equal length

  std::map<int,moab::Range> curves_and_counts;
  // how many curves
  std::cout << "There are " << curve_set.size() << " curves" <<  std::endl;
  rval = get_curve_lengths(curve_set, curves_and_counts);
  MB_CHK_SET_ERR(rval,"error from get_curve_lengths");
  
  // then compare by vertex handle - these are already shared
  
  rval = find_curve_pairs(curves_and_counts,coincident);
  // list of curves that are equivalent - their parents are surfaces that are the identical
  
  return moab::MB_SUCCESS;
}

// count the number of edges in each curve set and store in the map
moab::ErrorCode DAGMCTopology::get_curve_lengths(moab::Range curves,
						 std::map<int,moab::Range> &curves_and_counts) {
  // 
  moab::ErrorCode rval = moab::MB_FAILURE;
  // loop over the curves and get their lengths
  for ( moab::EntityHandle curve : curves ) {
    // get the id
    int id[1];
    rval = mbi->tag_get_data(id_tag,&(curve),1,id);

    // curve sets
    moab::Range edges;
    edges.clear();
    rval = mbi->get_entities_by_type(curve,moab::MBEDGE,edges);
    int count = edges.size();
    std::cout << count << std::endl;
    if( count == 0 ) {
      std::cout << "curve " << id[0] << " " << curve << " has no edges." << std::endl;
    } else {
      std::cout << "Curve " << id[0] << " " << curve << " has " << edges.size() << " edges." << std::endl;
    }
    // 
    moab::Range existing_range = curves_and_counts[count];
    existing_range.insert(curve);

    // add to the range;
    curves_and_counts[count] = existing_range;
  }

  // print the pair and curve counts
  for ( std::pair<int,moab::Range> pair : curves_and_counts ) {
    std::cout << "count: " << pair.first << " num_edges:" << pair.second.size() << std::endl;
  }
  return moab::MB_SUCCESS;
}

// compare the edges 
moab::ErrorCode DAGMCTopology::compare_edges(moab::Range edges1,
					     moab::Range edges2,
					     bool &same) {
  moab::ErrorCode rval = moab::MB_FAILURE;
  moab::Range connectivity1,connectivity2;
  // get the eh connectivity
  rval = mbi->get_connectivity(edges1,connectivity1);
  rval = mbi->get_connectivity(edges2,connectivity2);

  moab::Range merged = connectivity1;
  merged.merge(connectivity2);
  if ( merged.size() == connectivity1.size() ) {
    same = true;
  }
  return rval;
}

// compare two curves for their connectivity and return if they are comparible
moab::ErrorCode DAGMCTopology::compare_curve(moab::EntityHandle base_curve,
					     moab::EntityHandle other_curve,
					     bool &pair) {
  if ( base_curve == other_curve ) return moab::MB_SUCCESS;
  // the edge sets are the children of the curve sets
  moab::ErrorCode rval = moab::MB_FAILURE;
  moab::Range base_edges,other_edges;
  rval = mbi->get_entities_by_type(base_curve,moab::MBEDGE,base_edges);
  MB_CHK_SET_ERR(rval,"failed to get_entities_by_type for base_edges");
  rval = mbi->get_entities_by_type(other_curve,moab::MBEDGE,other_edges);
  MB_CHK_SET_ERR(rval,"failed to get_entities_by_type for other_edges");
 
  // they cannot possibly match
  if ( base_edges.size() != other_edges.size() ) return moab::MB_SUCCESS;

  // curves are already identical
  if(base_edges.contains(other_edges)) {
    pair = true;
    return rval;
  }

  // compare the edges - true if same
  rval = compare_edges(base_edges,other_edges,pair);
  MB_CHK_SET_ERR(rval,"failed to compare_edges");
  if (pair) {
    return rval;
  }
  
  return rval;
}

// compare a single curve with all curve and generate all matching pairs
moab::ErrorCode DAGMCTopology::compare_curves(moab::EntityHandle single_curve, moab::Range curves,
					      std::vector<merge_pairs_t> &pairs) {
  moab::ErrorCode rval = moab::MB_FAILURE;
  // loop over the curves and compare with the single curve
  for ( moab::EntityHandle curve : curves ) {
    bool match = false;
    // compare the single curve with all the curves
    rval = compare_curve(single_curve,curve,match);
    // if there is a match, add it to the pairs list
    if ( match ) {
      merge_pairs_t pair(single_curve,curve);
      pairs.insert(pairs.end(),pair);
    }
  }
  return rval;
}

// find all the matching curve pairs
moab::ErrorCode DAGMCTopology::find_curve_pairs(const std::map<int,moab::Range> curves,
						std::vector<merge_pairs_t> &pairs) {
  moab::Range curve_list;
  // build the master curve list
  for ( std::pair<int,moab::Range> pair : curves ) {
    moab::Range curve_tmp = pair.second;
    curve_list.merge(curve_tmp);
  }
  moab::ErrorCode rval = moab::MB_FAILURE;
  // loop over each curve and compare
  for ( moab::EntityHandle curve : curve_list ) {
    std::vector<merge_pairs_t> matches;
    // compare all the curves and generate a list 
    rval = compare_curves(curve,curve_list,matches);
    pairs.insert(pairs.end(),matches.begin(),matches.end());
  }
  
  return moab::MB_SUCCESS;
}

// merge surfaces that have been identified as being identical remove one of them
moab::ErrorCode DAGMCTopology::merge_coincident(const std::vector<merge_pairs_t> coincident, const bool del) {
  // the coincident list defines curves which have been idendtified to be identical
  // in terms of connectivity. procede through the list and create new links for the curves
  // and if selected delete the old links
  std::cout << "merge_coincident" << std::endl;
  
  moab::ErrorCode rval = moab::MB_FAILURE;
  for ( merge_pairs_t curves : coincident ) {

    // our preprocessing have already ensured that the master curve is the first entry
    // in the vector the second entry is the one to replace connectivity with
    
    // due to the unordered nature of moab lists 
    moab::Range parents1,parents2;

    rval = mbi->get_parent_meshsets(curves.first, parents1);
    rval = mbi->get_parent_meshsets(curves.second,parents2);

    // get the id
    //int id1[1],id2[1];
    //rval = mbi->tag_get_data(id_tag,&(curves.first),1,id1);
    //rval = mbi->tag_get_data(id_tag,&(curves.second),1,id2);

    /*
    std::cout << id1[0] << " " << id2[0] << std::endl;
    std::cout << curves.first << " " << curves.second << std::endl;
    std::cout << parents1[0] << " " << parents2[0] << std::endl;
    */
    
    // curves can have several parents, at points where 4 surfaces meet
    // for example - add pc relationship between parents2[0] and the master curve
    rval = mbi->add_parent_child(parents2[0],curves.first);
    rval = mbi->add_parent_child(parents1[0],curves.second);
    rval = mbi->add_parent_child(curves.first,curves.second);
    // delete the second curve
    if ( del ) {
      rval = mbi->delete_entities(&(curves.second),1);
    }
    
    // 
    MB_CHK_SET_ERR(rval,"Failed to add parent child relationship");
  }
  
  return moab::MB_SUCCESS;
}

moab::ErrorCode DAGMCTopology::save_file(const std::string filename) {
  moab::ErrorCode rval = moab::MB_FAILURE;
  rval = mbi->write_mesh(filename.c_str());
  return rval;
}
