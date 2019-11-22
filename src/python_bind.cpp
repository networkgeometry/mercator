#include "../include/embeddingS1.hpp"
#include <pybind11/pybind11.h>


void embed(std::string edgelist_filename, std::string rootname_output, std::string already_inferred_parameters_filename, bool fast_mode, bool screen_mode, bool post_kappa, bool quiet_mode, bool validation_mode, bool clean_mode, int seed, double beta)
{
  // Initialize graph object.
  embeddingS1_t the_graph;

  // Sets the edgelist filename.
  the_graph.EDGELIST_FILENAME = edgelist_filename;

  // Sets the output rootname.
  if(rootname_output.length() != 0)
  {
    the_graph.CUSTOM_OUTPUT_ROOTNAME_MODE = true;
    the_graph.ROOTNAME_OUTPUT = rootname_output;
  }

  // Activates the refine mode.
  if(already_inferred_parameters_filename.length() != 0)
  {
    the_graph.REFINE_MODE = true;
    the_graph.ALREADY_INFERRED_PARAMETERS_FILENAME = already_inferred_parameters_filename;
  }

  // Activates the fast mode.
  if(fast_mode)
  {
    the_graph.MAXIMIZATION_MODE = false;
  }

  // Activates the validation mode.
  if(validation_mode)
  {
    the_graph.VALIDATION_MODE = true;
    the_graph.CHARACTERIZATION_MODE = true;
  }

  // Deactivates the post-processing of kappas.
  if(!post_kappa)
  {
    the_graph.KAPPA_POST_INFERENCE_MODE = false;
  }

  // Activates the quiet mode.
  if(quiet_mode)
  {
    the_graph.QUIET_MODE = true;
  }

  // Activates the clean mode.
  if(clean_mode)
  {
    the_graph.CLEAN_RAW_OUTPUT_MODE = true;
  }

  // Activates the verbose mode.
  if(screen_mode)
  {
    the_graph.VERBOSE_MODE = true;
  }

  // Sets a custom seed, if required.
  if(seed != -1)
  {
    the_graph.CUSTOM_SEED = true;
    the_graph.SEED = seed;
  }

  // Sets a custom value of beta, if required.
  if(beta != -1)
  {
    the_graph.CUSTOM_BETA = true;
    the_graph.beta = beta;
  }

  // Performs the embedding.
  the_graph.embed();
}

namespace py = pybind11;

PYBIND11_MODULE(mercator, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: python_example

        .. autosummary::
           :toctree: _generate

           embed
    )pbdoc";

    m.def("embed", &embed, "",
          py::arg("edgelist_filename"),
          py::arg("output_name") = "",
          py::arg("inf_coord") = "",
          py::arg("fast_mode") = false,
          py::arg("screen_mode") = false,
          py::arg("post_kappa") = true,
          py::arg("quiet_mode") = false,
          py::arg("validation_mode") = false,
          py::arg("clean_mode") = false,
          py::arg("seed") = -1,
          py::arg("beta") = -1);


#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}
