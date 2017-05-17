#ifndef ELAST_EV_EW_TRACTION_BCS_HPP
#define ELAST_EV_EW_TRACTION_BCS_HPP

/** \file elast_ev_ew_traction_bcs.hpp */

#include <Phalanx_Evaluator_WithBaseImpl.hpp>
#include <Phalanx_Evaluator_Derived.hpp>

#include "goal_dimension.hpp"

/** \cond */
namespace goal {
class Indexer;
class Discretization;
class SolutionInfo;
}
/** \endcond */

namespace elast {

/** \cond */
using Teuchos::RCP;
using Teuchos::ParameterList;
/** \endcond */

template <typename EVALT, typename TRAITS>
class TractionEWBCs : public PHX::EvaluatorWithBaseImpl<TRAITS>,
                    public PHX::EvaluatorDerived<EVALT, TRAITS> {
 public:
  /** \cond */
  typedef typename TRAITS::SetupData SetupData;
  typedef typename TRAITS::PreEvalData PreEvalData;
  typedef typename TRAITS::PostEvalData PostEvalData;
  typedef typename TRAITS::EvalData EvalData;
  typedef typename EVALT::ScalarT ScalarT;
  using Elem = goal::Elem;
  using Node = goal::Node;
  using Dim = goal::Dim;
  using IP = goal::IP;
  using Dummy = goal::Dummy;
  /** \endcond */

  /** \brief Construct the traction BCs evaluator.
    * \param i The relevant \ref goal::Indexer structure.
    * \param p The traction boundary condition parameterlist, containing
    * entries of type Teuchos::Array<std::string> of the form:
    *  - [ side set name, t_x, t_y, t_z ]. */
  TractionEWBCs(RCP<goal::Indexer> i, RCP<const ParameterList> p, 
      RCP<goal::Field> z, RCP<goal::Field> zfine);

  /** \brief Finalize the field manager registration. */
  void postRegistrationSetup(SetupData d, PHX::FieldManager<TRAITS>& fm);

  /** \brief Grab the linear algebra data structures.
    * \param info The PreEvalData structure (\ref goal::SolutionInfo). */
  void preEvaluate(PreEvalData info);

  /** \brief Fill in the local multidimensional arrays. */
  void evaluateFields(EvalData workset);

 private:
  void validate_params();
  void apply_bc(EvalData workset, Teuchos::Array<std::string> const& a);

  RCP<const ParameterList> params;
  RCP<goal::Indexer> indexer;
  RCP<goal::Discretization> disc;
  RCP<goal::SolutionInfo> info;
  RCP<goal::Field> z;
  RCP<goal::Field> zfine;  
  int num_dims;
};

}  /* namespace elast */

#endif
