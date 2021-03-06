#include "goal_control.hpp"
#include "goal_ev_utils.hpp"
#include "goal_ev_gather_scalar.hpp"
#include "goal_ev_scalar_shape.hpp"
#include "goal_ev_interpolate_scalar.hpp"
#include "goal_ev_scatter_scalar.hpp"
#include "goal_ev_gather_vector.hpp"
#include "goal_ev_vector_shape.hpp"
#include "goal_ev_interpolate_vector.hpp"
#include "goal_ev_scatter_vector.hpp"
#include "goal_ev_dirichlet_bcs.hpp"
#include "goal_ev_qoi_ks.hpp"
#include "goal_ev_qoi_pnorm.hpp"
#include "goal_ev_qoi_scalar_point.hpp"
#include "goal_ev_qoi_vector_point.hpp"
#include "goal_ev_scatter_functional.hpp"
#include "goal_ev_dual_scalar_weight.hpp"
#include "goal_ev_dual_vector_weight.hpp"
#include "goal_ev_scalar_error.hpp"
#include "goal_ev_vector_error.hpp"
#include "goal_field.hpp"
#include "goal_indexer.hpp"

namespace goal {

using Teuchos::rcp;

template <typename EvalT>
void register_dof(RCP<Field> u, RCP<Indexer> i, FieldManager fm) {
  auto type = u->get_value_type();
  RCP<PHX::Evaluator<goal::Traits> > gather;
  RCP<PHX::Evaluator<goal::Traits> > shape;
  RCP<PHX::Evaluator<goal::Traits> > interpolate;
  if (type == SCALAR) {
    gather = rcp(new GatherScalar<EvalT, goal::Traits>(u, i));
    shape = rcp(new ScalarShape<EvalT, goal::Traits>(u));
    interpolate = rcp(new InterpolateScalar<EvalT, goal::Traits>(u));
  } else if (type == VECTOR) {
    gather = rcp(new GatherVector<EvalT, goal::Traits>(u, i));
    shape = rcp(new VectorShape<EvalT, goal::Traits>(u));
    interpolate = rcp(new InterpolateVector<EvalT, goal::Traits>(u));
  } else
    fail("cannot register dof");
  fm->registerEvaluator<EvalT>(gather);
  fm->registerEvaluator<EvalT>(shape);
  fm->registerEvaluator<EvalT>(interpolate);
}

template <typename EvalT>
void require_primal_scatter(RCP<Field> u, RCP<Indexer> i, FieldManager fm) {
  auto type = u->get_value_type();
  RCP<PHX::Evaluator<goal::Traits> > scatter;
  if (type == SCALAR)
    scatter = rcp(new ScatterScalar<EvalT, goal::Traits>(u, i, false));
  else if (type == VECTOR)
    scatter = rcp(new ScatterVector<EvalT, goal::Traits>(u, i, false));
  else
    fail("cannot require primal scatter");
  fm->registerEvaluator<EvalT>(scatter);
  auto op = scatter->evaluatedFields()[0];
  fm->requireField<EvalT>(*op);
}

template <typename EvalT>
void register_dual(RCP<Field> z, RCP<Field> z_fine, FieldManager fm) {
  auto z_type = z->get_value_type();
  auto z_fine_type = z_fine->get_value_type();
  assert (z_type == z_fine_type);
  RCP<PHX::Evaluator<goal::Traits> > weight;
  if (z_type == SCALAR)
    weight = rcp(new DualScalarWeight<EvalT, goal::Traits>(z, z_fine));
  else if (z_type == VECTOR)
    weight = rcp(new DualVectorWeight<EvalT, goal::Traits>(z, z_fine));
  else
    fail("cannot register dual weight");
  fm->registerEvaluator<EvalT>(weight);
}

template <typename EvalT>
void require_adjoint_scatter(RCP<Field> u, RCP<Indexer> i, FieldManager fm) {
  auto type = u->get_value_type();
  RCP<PHX::Evaluator<goal::Traits> > scatter;
  if (type == SCALAR)
    scatter = rcp(new ScatterScalar<EvalT, goal::Traits>(u, i, true));
  else if (type == VECTOR)
    scatter = rcp(new ScatterVector<EvalT, goal::Traits>(u, i, true));
  else
    fail("cannot require adjoint scatter");
  fm->registerEvaluator<EvalT>(scatter);
  auto op = scatter->evaluatedFields()[0];
  fm->requireField<EvalT>(*op);
}

template <typename EvalT>
void require_dbc(RCP<const ParameterList> p, RCP<Indexer> i, bool condense,
    bool adj, FieldManager fm) {
  auto ev =
    rcp(new goal::DirichletBCs<EvalT, goal::Traits>(p, i, condense, adj));
  fm->registerEvaluator<EvalT>(ev);
  fm->requireField<EvalT>(*ev->evaluatedFields()[0]);
}

static RCP<const ParameterList> get_valid_p_qoi_params() {
  auto p = rcp(new ParameterList);
  p->set<std::string>("name", "");
  p->set<std::string>("scalar name", "");
  p->set<double>("p", 0);
  p->set<double>("m", 0);
  return p;
}

void require_qoi_ks(RCP<const ParameterList> p, RCP<Field> u,
    RCP<Indexer> i, FieldManager fm) {
  using T = goal::Traits;
  using J = goal::Traits::Jacobian;
  p->validateParameters(*get_valid_p_qoi_params(), 0);
  auto qoi_n = p->get<std::string>("scalar name");
  auto qoi_p = p->get<double>("p");
  auto qoi_m = p->get<double>("m");
  auto qoi = rcp(new QoIKS<J, T>(u, qoi_n, qoi_p, qoi_m));
  auto scatter = rcp(new ScatterFunctional<J, T>(u, i, "KS Functional"));
  auto op = scatter->evaluatedFields()[0];
  fm->registerEvaluator<J>(qoi);
  fm->registerEvaluator<J>(scatter);
  fm->requireField<J>(*op);
}

void require_qoi_pnorm(RCP<const ParameterList> p, RCP<Field> u,
    RCP<Indexer> i, FieldManager fm) {
  using T = goal::Traits;
  using J = goal::Traits::Jacobian;
  p->validateParameters(*get_valid_p_qoi_params(), 0);
  auto qoi_n = p->get<std::string>("scalar name");
  auto qoi_p = p->get<double>("p");
  auto qoi_m = p->get<double>("m");
  auto qoi = rcp(new QoIPNorm<J, T>(u, qoi_n, qoi_p, qoi_m));
  auto scatter = rcp(new ScatterFunctional<J, T>(u, i, "P-Norm Functional"));
  auto op = scatter->evaluatedFields()[0];
  fm->registerEvaluator<J>(qoi);
  fm->registerEvaluator<J>(scatter);
  fm->requireField<J>(*op);
}

void require_qoi_scalar_point(RCP<Field> u, RCP<Indexer> i,
    std::string const& set, FieldManager fm) {
  using T = goal::Traits;
  using J = goal::Traits::Jacobian;
  auto qoi = rcp(new QoIScalarPoint<J, T>(u, i, set));
  auto scatter =
    rcp(new ScatterFunctional<J, T>(u, i, "Point-Wise " + u->get_name()));
  auto op = scatter->evaluatedFields()[0];
  fm->registerEvaluator<J>(qoi);
  fm->registerEvaluator<J>(scatter);
  fm->requireField<J>(*op);
}

void require_qoi_vector_point(RCP<Field> u, RCP<Indexer> i, int c,
    std::string const& set, FieldManager fm) {
  using T = goal::Traits;
  using J = goal::Traits::Jacobian;
  auto qoi = rcp(new QoIVectorPoint<J, T>(u, i, set));
  auto scatter =
    rcp(new ScatterFunctional<J, T>(u, i, "Point-Wise " + u->get_name()));
  auto op = scatter->evaluatedFields()[0];
  fm->registerEvaluator<J>(qoi);
  fm->registerEvaluator<J>(scatter);
  fm->requireField<J>(*op);
  (void)c;
}

void require_qoi(RCP<const ParameterList> p, RCP<Field> u,
    RCP<Indexer> i, FieldManager fm) {
  auto n = p->get<std::string>("name");
  if (n == "ks") require_qoi_ks(p, u, i, fm);
  else if (n == "p norm") require_qoi_pnorm(p, u, i, fm);
  else fail("unknown qoi: %s", n.c_str());
}

void require_error(RCP<Field> u, RCP<Field> e, RCP<Indexer> i,
    FieldManager fm) {
  using T = goal::Traits;
  using R = goal::Traits::Residual;
  auto type = u->get_value_type();
  RCP<PHX::Evaluator<goal::Traits> > error;
  if (type == SCALAR)
    error = rcp(new ScalarError<R, T>(u, e, i));
  else if (type == VECTOR)
    error = rcp(new VectorError<R, T>(u, e, i));
  else
    fail("cannot require error");
  auto op = error->evaluatedFields()[0];
  fm->registerEvaluator<R>(error);
  fm->requireField<R>(*op);
}

template void register_dof<goal::Traits::Residual>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);
template void register_dof<goal::Traits::Jacobian>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);
template void register_dual<goal::Traits::Residual>(RCP<Field> z, RCP<Field> z_fine, FieldManager fm);
template void require_primal_scatter<goal::Traits::Residual>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);
template void require_primal_scatter<goal::Traits::Jacobian>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);
template void require_adjoint_scatter<goal::Traits::Residual>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);
template void require_adjoint_scatter<goal::Traits::Jacobian>(RCP<Field> u, RCP<Indexer> i, FieldManager fm);
template void require_dbc<goal::Traits::Residual>(RCP<const ParameterList> p, RCP<Indexer> i, bool condense, bool adj, FieldManager fm);
template void require_dbc<goal::Traits::Jacobian>(RCP<const ParameterList> p, RCP<Indexer> i, bool condense, bool adj, FieldManager fm);

}  /* namespace goal */
