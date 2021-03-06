#include <PCU.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <Teuchos_ParameterList.hpp>

#include "goal_discretization.hpp"
#include "goal_output.hpp"

namespace goal {

using Teuchos::rcp;

static RCP<ParameterList> get_valid_params() {
  using Teuchos::Array;
  Array<std::string> dummy(0);
  auto p = rcp(new ParameterList);
  p->set<std::string>("out file", "");
  p->set<int>("interval", 1);
  p->set<bool>("turn off", false);
  p->set<Array<std::string> >("interpolate", dummy);
  return p;
}

static void validate_params(RCP<const ParameterList> p) {
  assert(p->isType<std::string>("out file"));
  p->validateParameters(*get_valid_params(), 0);
}

static void write_initial_pvd(std::string const& name, int& pos) {
  if (!PCU_Comm_Self()) {
    auto pvd = name + ".pvd";
    std::fstream pvdf;
    pvdf.open(pvd.c_str(), std::ios::out);
    pvdf << "<VTKFile type=\"Collection\" version=\"0.1\">" << std::endl;
    pvdf << "  <Collection>" << std::endl;
    pos = pvdf.tellp();
    pvdf << "  </Collection>" << std::endl;
    pvdf << "</VTKFile>" << std::endl;
    pvdf.close();
  }
}

Output::Output(RCP<const ParameterList> p, RCP<Discretization> d) {
  params = p;
  disc = d;
  validate_params(params);
  turn_off = false;
  interval = 1;
  pos = 0;
  index = 0;
  name = params->get<std::string>("out file");
  if (params->isParameter("turn off")) turn_off = params->get<bool>("turn off");
  if (params->isParameter("interval")) interval = params->get<int>("interval");
  if (params->isParameter("interpolate"))
    fields = params->get<Teuchos::Array<std::string> >("interpolate");
  if (!turn_off) write_initial_pvd(name, pos);
}

static void update_pvd(
    std::string const& name, std::string const& vtu, int& pos, const double t) {
  if (!PCU_Comm_Self()) {
    std::string pvd = name + ".pvd";
    std::fstream pvdf;
    pvdf.open(pvd.c_str(), std::ios::out | std::ios::in);
    pvdf.seekp(pos);
    pvdf << "    <DataSet timestep=\"" << t << "\" group=\"\" ";
    pvdf << "part=\"0\" file=\"" << vtu << "/" << vtu;
    pvdf << ".pvtu\"/>" << std::endl;
    pos = pvdf.tellp();
    pvdf << "  </Collection>" << std::endl;
    pvdf << "</VTKFile>" << std::endl;
    pvdf.close();
  }
}

static void interpolate(
    apf::Mesh* m, Teuchos::Array<std::string> const& names) {
  for (int i = 0; i < names.size(); ++i) {
    auto name = names[i].c_str();
    apf::Field* f = m->findField(name);
    assert(f);
    auto iname = names[i] + "_interp";
    auto type = apf::getValueType(f);
    apf::Field* g = apf::createFieldOn(m, iname.c_str(), type);
    apf::projectField(g, f);
  }
}

static void destroy(
    apf::Mesh* m, Teuchos::Array<std::string> const& names) {
  for (int i = 0; i < names.size(); ++i) {
    auto name = names[i] + "_interp";
    apf::Field* g = m->findField(name.c_str());
    assert(g);
    apf::destroyField(g);
  }
}

void Output::write_vtk(const double t) {
  auto m = disc->get_apf_mesh();
  std::ostringstream oss;
  oss << name << "_" << index;
  std::string vtu = oss.str();
  update_pvd(name, vtu, pos, t);
  interpolate(m, fields);
  apf::writeVtkFiles(vtu.c_str(), m);
  destroy(m, fields);
  ++index;
}

void Output::write(const double t) {
  if (turn_off) return;
  static int my_out_interval = 0;
  if (my_out_interval++ % interval) return;
  write_vtk(t);
}

}  /* namespace goal */
