#ifndef GOAL_INDEXER_HPP
#define GOAL_INDEXER_HPP

/** \file goal_indexer.hpp */

#include "goal_data_types.hpp"

/** \cond */
namespace apf {
struct Node;
class Mesh;
class Field;
class MeshEntity;
template <class T>
class NumberingOf;
typedef NumberingOf<int> Numbering;
typedef NumberingOf<long> GlobalNumbering;
}
/** \endcond */

namespace goal {

using Teuchos::RCP;

/** \cond */
class Field;
class Discretization;
/** \endcond */

/** \brief The DOF indexer interface.
  * \details The indexer structure is responsible for building the maps and
  * graphs that describe the parallel distrubution of local to global data
  * in the linear algebra data structures. It does this by numbering the
  * DOFs associated with an input vector of DOF fields. This class is also
  * responsible for building node sets that are associated with the
  * specification of Dirichlet boundary conditions. */
class Indexer {
 public:
  /** \brief The indexer constructor.
    * \param d The relevant \ref goal::Discretization structure.
    * \param f A vector of DOF \ref goal::Field s. */
  Indexer(RCP<Discretization> d, std::vector<RCP<Field> > f);

  /** \brief Destroy the indexer.
    * \details This destroys all DOF numberings, maps, and graphs that
    * the indexer has constructed. */
  ~Indexer();

  /** \brief Set the current element block to operate over.
    * \details This is useful in scenarios where different element topologies
    * are associated with different element blocks when element-level index
    * information is queried. */
  void set_elem_block(const int idx);

  /** \brief Get the owned DOF map. */
  RCP<const Map> get_owned_map() { return owned_map; }

  /** \brief Get the ghost DOF map. */
  RCP<const Map> get_ghost_map() { return ghost_map; }

  /** \brief Get the owned CRS element connectivity graph. */
  RCP<const Graph> get_owned_graph() { return owned_graph; }

  /** \brief Get the ghost CRS element connectivity graph. */
  RCP<const Graph> get_ghost_graph() { return ghost_graph; }

  /** \brief Get the ith \ref goal::Field associated with this indexer. */
  RCP<Field> get_field(const int idx) { return fields[idx]; }

  /** \brief Get the underlying \ref goal::Discretization of this indexer. */
  RCP<Discretization> get_discretization() { return disc; }

  /** \brief Get the current element block index.
    * \details Can be changed using \ref goal::Indexer::set_elem_block. */
  int get_elem_block_idx() const { return elem_block_idx; }

  /** \brief Get the total number of DOF fields. */
  int get_num_dof_fields() const;

  /** \brief Get the total number of elemental DOF fields.
    * \details For elements in the current operating element block. */
  int get_num_total_elem_dofs() const;

  /** \brief Get the DOF index associated with \ref goal::Field::get_name. */
  int get_dof_idx(std::string const& name) const;

  /** \brief Get the offset into the elemental DOFs.
    * \param idx The \ref goal::Field DOF index.
    * \param n The elemental node associated with the field.
    * \param c The field component associated with the node. */
  int get_elem_dof_offset(const int idx, const int n, const int c);

  /** \brief Get the ghost local row ID.
    * \param idx The \ref goal::Field DOF index.
    * \param e The relevant mesh entity.
    * \param n The local node associated with the mesh entity.
    * \param c The field component associated with local node. */
  LO get_ghost_lid(const int idx, apf::MeshEntity* e, const int n, const int c);

  /** \brief Get the owned local row ID.
    * \param idx The \ref goal::Field DOF index.
    * \param n An APF node data strcture (a MeshEntity* + a local node).
    * \param c The field component associated with the APF node. */
  LO get_owned_lid(const int idx, apf::Node const& n, const int c);

  /** \brief Get the elemental ghost local row IDs.
    * \param e The relevant mesh entity.
    * \param lids The returned local IDs. */
  void get_ghost_lids(apf::MeshEntity* e, std::vector<LO>& lids);

  /** \brief Get the number of node sets. */
  int get_num_node_sets() const;

  /** \brief Get the nodes in a given node set for a given DOF field.
    * \param set The relevant node set name.
    * \param idx The relevant \ref goal::Field DOF index. */
  std::vector<apf::Node> const& get_node_set_nodes(
      std::string const& set, const int idx);

  /** \brief Get the underlying APF mesh used for this discretization.
    * \details This is mainly for convenience. */
  apf::Mesh* get_apf_mesh();

  /** \brief Get the owned APF local numbering for a given DOF field.
    * \param idx The relevant \ref goal::Field DOF index. */
  apf::Numbering* get_apf_numbering(const int idx);

  /** \brief Get the ghost APF local numbering for a given DOF field.
    * \param idx The relevant \ref goal::Field DOF index. */
  apf::Numbering* get_ghost_apf_numbering(const int idx);

 private:
  void compute_owned_map();
  void compute_ghost_map();
  void compute_graphs();
  void compute_node_sets();
  RCP<Discretization> disc;
  std::vector<RCP<Field> > fields;
  std::vector<apf::Field*> apf_fields;
  std::vector<apf::Numbering*> owned_numberings;
  std::vector<apf::Numbering*> ghost_numberings;
  std::vector<apf::GlobalNumbering*> global_numberings;
  std::map<std::string, std::vector<std::vector<apf::Node> > > node_sets;
  RCP<const Comm> comm;
  RCP<const Map> owned_map;
  RCP<const Map> ghost_map;
  RCP<Graph> owned_graph;
  RCP<Graph> ghost_graph;
  int elem_block_idx;
};

}  // namesapce goal

#endif
