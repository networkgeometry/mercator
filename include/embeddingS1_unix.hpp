/*
 *
 * Provides the functions related to UNIX operating system.
 *
 * Author:  Antoine Allard
 * WWW:     antoineallard.info
 * Date:    September 2017
 *
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 */

#ifndef EMBEDDINGS1_UNIX_HPP_INCLUDED
#define EMBEDDINGS1_UNIX_HPP_INCLUDED

// Standard Template Library
#include <cstdlib>
#include <iostream>
#include <string>
// Operating System
#include <unistd.h>
// embeddingS1
#include "embeddingS1.hpp"





// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// Prints the information on the way the command line UNIX program should be used.
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void print_usage()
{
  std::cout << std::endl;
  std::cout << "NAME"                                                                                  << std::endl;
  std::cout << "\tMercator: Inference of high-quality embeddings of complex networks into the"         << std::endl;
  std::cout << "\t          hyperbolic disk"                                                           << std::endl;
  std::cout << std::endl;
  std::cout << "SYNOPSIS"                                                                              << std::endl;
  std::cout << "\tmercator [options] <edgelist_filename>"                                              << std::endl;
  std::cout << std::endl;
  std::cout << "INPUT"                                                                                 << std::endl;
  std::cout << "\tThe structure of the graph is provided by a text file containing it edgelist. Each"  << std::endl;
  std::cout << "\tline in the file corresponds to an edge in the graph (i.e., [VERTEX1] [VERTEX2])."   << std::endl;
  std::cout << "\t  - The name of the vertices need not be integers (they are stored as std::string)." << std::endl;
  std::cout << "\t  - Directed graphs will be converted to undirected."                                << std::endl;
  std::cout << "\t  - Multiple edges, self-loops and weights will be ignored."                         << std::endl;
  std::cout << "\t  - Lines starting with '# ' are ignored (i.e., comments)."                          << std::endl;
  // std::cout << std::endl;
  // std::cout << "\tOuput: a file containing the inferred coordinates (kappa, theta) for each vertex." << std::endl;
  std::cout << std::endl;
}





// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// Prints the information about the options of the command line UNIX program should be used.
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void print_help()
{
  std::cout << std::endl;
  std::cout << "The following options are available:"                                                 << std::endl;
  std::cout << "\t-a             Screen mode. Program outputs details about its progress on screen"   << std::endl;
  std::cout << "\t               (through std::clog) instead of in a log file. Useful to gather all"  << std::endl;
  std::cout << "\t               output in a single file if mercator is a subroutine of a script."    << std::endl;
  std::cout << "\t-b [VALUE]     Specify the value for beta to be used for the embedding. By "        << std::endl;
  std::cout << "\t               default the program infers the value of beta based on the average"   << std::endl;
  std::cout << "\t               local clustering coefficient of the original edgelist."              << std::endl;
  std::cout << "\t-c             Clean mode. Writes the inferred coordinates in clean file without"   << std::endl;
  std::cout << "\t               any lines starting by # to facilitate their use in custom computer"  << std::endl;
  std::cout << "\t               programs."                                                           << std::endl;
  std::cout << "\t-f             Fast mode. Does not infer the positions based on likelihood"         << std::endl;
  std::cout << "\t               maximization, rather uses only the EigenMap method."                 << std::endl;
  std::cout << "\t-g             Forces the condition beta > 1, the so called geometric phase,"       << std::endl;
  std::cout << "\t               ignoring the quasi- and non-geometric phases at beta < 1."           << std::endl;
  std::cout << "\t               The default does take these regions into account."                   << std::endl;
  std::cout << "\t-k             No post-processing of the values of kappa based on the inferred"     << std::endl;
  std::cout << "\t               angular positions (theta) resulting in every vertices with the same" << std::endl;
  std::cout << "\t               degree ending at the same radial position in the hyperbolic disk."   << std::endl;
  // std::cout << "\t-h             Print this message on screen and exit."                              << std::endl;
  std::cout << "\t-m             Uses an analytic approximation for mu that holds for N >> 1."        << std::endl;
  std::cout << "\t               The default calculates mu numerically, taking finite size effects "  << std::endl;
  std::cout << "\t               into account."                                                       << std::endl; 
  std::cout << "\t-o [ROOTNAME]  Specify the rootname used for all output files. Default: uses the"   << std::endl;
  std::cout << "\t               rootname of the edgelist file as (i.e., rootname.edge)."             << std::endl;
  std::cout << "\t-r [FILENAME]  Refine mode. Reads the inferred positions from a previous run of"    << std::endl;
  std::cout << "\t               this program (file *.inf_coord) and refines the inferred positions." << std::endl;
  std::cout << "\t-q             Quiet mode. Program does not output information about the network"   << std::endl;
  std::cout << "\t               and the embedding procedure."                                        << std::endl;
  std::cout << "\t-s [SEED]      Program uses a custom seed for the random number generator."         << std::endl;
  std::cout << "\t               Default: EPOCH."                                                     << std::endl;
  std::cout << "\t-v             Validation mode. Validates and characterizes the inferred random"    << std::endl;
  std::cout << "\t               network ensemble."                                                   << std::endl;
  std::cout << std::endl;
  // std::cout << "Please see file embeddingS1.hpp for more options (will require recompilation)."       << std::endl;
  // std::cout << std::endl;
}





// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
// Parses the options (for UNIX-like command line use) and returns the filename of the edgelist.
// =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
void parse_options(int argc , char *argv[], embeddingS1_t &the_graph)
{
  // Shows the options if no argument is given.
  if(argc == 1)
  {
    print_usage();
    print_help();
    std::exit(0);
  }

  // <edgelist_filename>
  the_graph.EDGELIST_FILENAME = argv[argc - 1];

  // Parsing options.
  int opt;
  while ((opt = getopt(argc,argv,"ab:gmcfko:r:qs:v")) != -1)
  {
    switch(opt)
    {
      case 'a':
        the_graph.VERBOSE_MODE = true;
        break;

      case 'b':
        the_graph.CUSTOM_BETA = true;
        the_graph.beta = std::stod(optarg);
        break;

      case 'c':
        the_graph.CLEAN_RAW_OUTPUT_MODE = true;
        break;

      case 'f':
        the_graph.MAXIMIZATION_MODE = false;
        break;

      // case 'h':
      //   print_usage();
      //   print_help();
      //   std::exit(0);

      case 'g':
        the_graph.ALL_BETA_MODE = false;
        break;

      case 'm':
        the_graph.NUMERIC_MU_MODE = false;
        break;

      case 'k':
        the_graph.KAPPA_POST_INFERENCE_MODE = false;
        break;

      case 'o':
        the_graph.CUSTOM_OUTPUT_ROOTNAME_MODE = true;
        the_graph.ROOTNAME_OUTPUT = optarg;
        break;

      case 'r':
        the_graph.REFINE_MODE = true;
        the_graph.ALREADY_INFERRED_PARAMETERS_FILENAME = optarg;
        break;

      case 'q':
        the_graph.QUIET_MODE = true;
        break;

      case 's':
        the_graph.CUSTOM_SEED = true;
        the_graph.SEED = std::stoi(optarg);
        break;

      case 'v':
        the_graph.VALIDATION_MODE = true;
        the_graph.CHARACTERIZATION_MODE = true;
        break;

      default:
        print_usage();
        print_help();
        std::exit(0);
    }
  }
}

#endif // EMBEDDINGS1_UNIX_HPP_INCLUDED
