#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/shared.h>

#include <vector>
#include <map>
#include "minimizers/scaling_functions.cpp"
#include "minimizers/G_wI_LM.h"

using namespace boost::python;

namespace xscale{
namespace boost_python { namespace {

  void
  xscale_module() {
    using namespace boost::python;
    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;
    typedef xscale::xscale_G_wI_LM wt;
    def("Load_object", &Load_object);
    def("read_objects", &read_objects);
    def("get_f_g", &get_f_g);
    def("calculate_residual_jacobian_wI_G", &calculate_residual_jacobian_wI_G);
    typedef return_internal_reference<> rir;
    scitbx::af::boost_python::shared_wrapper<Miller, rir>::wrap("Miller");
    scitbx::af::boost_python::shared_wrapper<Frame, rir>::wrap("Frame");

    class_<wt,
           bases<scitbx::lstbx::normal_equations::non_linear_ls<double> > >("xscale_G_wI_LM", init<int>())
             .def("build_up_normal_matrix",&wt::build_up_normal_matrix)
             .def("set_cpp_data",&wt::set_cpp_data)
//             .def("build_up_sparse_normal_matrix", &wt::build_up_sparse_normal_matrix)
    ;
  }

}
}} // namespace xfel::boost_python::<anonymous>

BOOST_PYTHON_MODULE(xscale_ext)
{
  xscale::boost_python::xscale_module();

}
