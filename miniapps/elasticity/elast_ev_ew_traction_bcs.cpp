#include <apf.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <goal_control.hpp>
#include <goal_discretization.hpp>
#include <goal_field.hpp>
#include <goal_indexer.hpp>
#include <goal_solution_info.hpp>
#include <goal_traits.hpp>
#include <goal_workset.hpp>
#include <Phalanx_DataLayout_MDALayout.hpp>

#include "elast_ev_ew_traction_bcs.hpp"

namespace elast {

using Teuchos::rcp;
using Teuchos::rcpFromRef;
using Teuchos::ParameterEntry;
using Teuchos::getValue;
using Teuchos::Array;

template <typename EVALT, typename TRAITS>
void TractionEWBCs<EVALT, TRAITS>::validate_params() {
  for (auto it = params->begin(); it != params->end(); ++it) {
    auto entry = params->entry(it);
    auto a = getValue<Array<std::string> >(entry);
    assert(a.size() == num_dims + 1);
    auto set = a[0];
    disc->get_sides(set);
  }
}

template <typename EVALT, typename TRAITS>
TractionEWBCs<EVALT, TRAITS>::TractionEWBCs(
    RCP<goal::Indexer> i, RCP<const ParameterList> p,
    RCP<goal::Field> z_, RCP<goal::Field> zfine_) {
  params = p;
  indexer = i;
  assert(indexer->get_num_dof_fields() == 1);
  disc = indexer->get_discretization();
  validate_params();
  z = z_;
  zfine = zfine_;
  num_dims = zfine->get_num_dims();
  auto name = "Traction BCs";
  PHX::Tag<ScalarT> op(name, rcp(new PHX::MDALayout<Dummy>(0)));
  this->addEvaluatedField(op);
  this->setName(name);
}

template <typename EVALT, typename TRAITS>
void TractionEWBCs<EVALT, TRAITS>::postRegistrationSetup(
    SetupData d, PHX::FieldManager<TRAITS>& fm) {
  (void)d;
  (void)fm;
}

template <typename EVALT, typename TRAITS>
void TractionEWBCs<EVALT, TRAITS>::preEvaluate(PreEvalData i) {
  info = rcpFromRef(i);
  assert(info->ghost->R != Teuchos::null);
}

template <typename EVALT, typename TRAITS>
void TractionEWBCs<EVALT, TRAITS>::apply_bc(
    EvalData workset, Teuchos::Array<std::string> const& a) {
  int idx = 0;
  auto t = workset.t_current;
  auto R = info->ghost->R;
  auto set = a[0];
  auto sides = disc->get_sides(set);
  auto u = indexer->get_field(idx);
  auto q_degree = u->get_q_degree();
  auto basis = apf::getLagrange(1);
  auto mesh = indexer->get_apf_mesh();
  apf::Vector3 xi(0, 0, 0);
  apf::Vector3 x(0, 0, 0);
  apf::Vector3 traction(0, 0, 0);
  apf::Vector3 z_val;
  apf::Vector3 zfine_val;
  apf::Vector3 zfine_minus_z;
  apf::NewArray<double> BF;
  std::vector<goal::LO> numbers;
  auto apf_z = z->get_apf_field();
  auto apf_zfine = zfine->get_apf_field();

  for (size_t i = 0; i < sides.size(); ++i) {
    auto f = sides[i];
    auto me = apf::createMeshElement(mesh, f);
    indexer->get_ghost_lids(f, numbers);
    int num_ips = apf::countIntPoints(me, q_degree);
    auto es = basis->getEntityShape(mesh->getType(f));
    int num_nodes = es->countNodes();
    auto z_elem = apf::createElement(apf_z, me);
    auto zfine_elem = apf::createElement(apf_zfine, me);
    for (int ip = 0; ip < num_ips; ++ip) {
      apf::getIntPoint(me, q_degree, ip, xi);
      apf::mapLocalToGlobal(me, xi, x);
      apf::getBF(basis, me, xi, BF);
      apf::getVector(z_elem, xi, z_val);
      apf::getVector(zfine_elem, xi, zfine_val);      
      double w = apf::getIntWeight(me, q_degree, ip);
      double dv = apf::getDV(me, xi);
      for (int i = 0; i < num_dims; ++i) {
        traction[i] = goal::eval(a[i + 1], x[0], x[1], x[2], t);
      }
      zfine_minus_z = zfine_val - z_val; 
      int dof = 0;
      for (int node = 0; node < num_nodes; ++node) {
        for (int dim = 0; dim < num_dims; ++dim) {
          goal::LO row = numbers[dof++];
          double val = -BF[node] * zfine_minus_z[dim]  * traction[dim] * w * dv;
          R->sumIntoLocalValue(row, val);
        }
      }
    }
    apf::destroyElement(zfine_elem);
    apf::destroyElement(z_elem);    
    apf::destroyMeshElement(me);
  }
}

template <typename EVALT, typename TRAITS>
void TractionEWBCs<EVALT, TRAITS>::evaluateFields(EvalData workset) {
  for (auto i = params->begin(); i != params->end(); ++i) {
    ParameterEntry const& entry = params->entry(i);
    Array<std::string> a = getValue<Array<std::string> >(entry);
    apply_bc(workset, a);
  }
}

template class TractionEWBCs<goal::Traits::Residual, goal::Traits>;
template class TractionEWBCs<goal::Traits::Jacobian, goal::Traits>;

}  /* namespace elast */
