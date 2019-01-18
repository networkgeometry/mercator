# Mercator
Inference of high-quality embeddings of complex networks into the hyperbolic disk



## Table of content

1. Installation
    * [Command line executable](#command-line-executable)
    * [Python module](#python-module)
2. Usage
    * [Input files](#format-of-input-files)
    * [Running the code](#running-the-code)
    * [Ouput files](#output-files)
    * [Options](#options)



## Installation

Two interfaces of the program are available: a command line executable and a Python module.


### Command line executable

Requirements
* A C++11 (or newer) compliant compiler
* The header `unistd.h`.

```
# Unix (Linux / MAC OS)
g++ -O3 -std=c++11 -I include/ src/embeddingS1_unix.cpp -o mercator
```


### Python module

Requirements
* A C++11 (or newer) compliant compiler
* python 3.x
* [pybind11] (will be installed automatically)

```
# Unix (Linux / MAC OS)
pip3 install setup.py
```


## Usage


### Format of input files

The structure of the network to be embedded is passed to the program via a file containing its edgelist (one link per line). The edgelist file consists in a simple text file with the following convention
```
# lines beginning with "#" are ignored (comments).
# note that nodes' name must be separated by at least one white space.
# there may be white space at the beginning of a line.
[name of node1]  [name of node2]  [remaining information will be ignored]
[name of node2]  [name of node3]  [remaining information will be ignored]
[name of node4]  [name of node5]  [remaining information will be ignored]
# comments can be inserted between links
[name of node5]  [name of node6]  [remaining information will be ignored]
[name of node7]  [name of node6]  [remaining information will be ignored]
...
```
Note that the nodes' name will be imported as `std::string` and can therefore be virtually anything as long as they do not include white spaces (i.e., there is not need for the nodes to be identified by contiguous integers).

**IMPORTANT**: this class only considers **simple undirected** networks **without self-loops**. Any multiple edges (e.g., if the graph is originally directed) or self-loops will be ignored.

**IMPORTANT**: in the actual version of the code, the network must have **only one component**.


### Running the code

Running `mercator` is quite straightforward
```
# Command line
./mercator <edgelist_filename>

# Python module
mercator.embed(<edgelist_filename>)
```


### Output files

The program outputs several files during and after the embedding procedure

* `*.inf_coord`: Contains information about the embedding procedure, the inferred parameters and the inferred positions.
* `*.inf_log`: Contains detailled information about the embedding procedure.


### Options

Several options are provided to adjust the embedding procedure to specific needs and can be combined

* [Custom output filename](#custom-output-filename)
* [Custom value for beta](#custom-value-for-beta)
* [Custom value for the seed of the random number generator](#custom-value-for-the-seed-of-the-random-number-generator)
* [Fast mode](#fast-mode)
* [Post-processing of the inferred values of the radial positions](#post-processing-of-the-inferred-values-of-the-radial-positions)
* [Quiet mode](#quiet-mode)
* [Refine mode](#refine-mode)
* [Screen mode](#screen-mode)
* [Validation mode](#validation-mode)

#### Custom output filename

All generated files are named `<output_rootname>.<extension>` and a custom `<output_rootname>` can be provided. If none is provided, the `<output_rootname>` is extracted from the `<edgelist_filename>` by removing its extension, otherwise the full `<edgelist_filename>` is used as `<output_rootname>` if `<edgelist_filename>`does not have any extension.

```
# Command line
./mercator -o <custom_output_rootname> <edgelist_filename>

# Python module
mercator.embed(<edgelist_filename>, output_name=<custom_output_rootname>)
```

#### Custom value for beta

A custom value for the parameter `beta` can be provided. By default a value for beta is inferred to reproduce the average local clustering coefficient of the original edgelist.

```
# Command line
./mercator -b <beta_value> <edgelist_filename>

# Python module
mercator.embed(<edgelist_filename>, beta=<beta_value>)
```

#### Custom value for the seed of the random number generator

A custom seed for the random number generator can be provided (useful when several embeddings are performed in parallel). By default, `EPOCH` is used.

```
# Command line
./mercator -s <seed_value> <edgelist_filename>

# Python module
mercator.embed(<edgelist_filename>, seed=<seed_value>)
```

#### Fast mode

Skip the likelihood maximization step (i.e., only infers the positions using the EigenMap methods). Default is **`false`**.

```
# Command line
./mercator -f <edgelist_filename>

# Python module
mercator.embed(<edgelist_filename>, fast_mode=True)
```

#### Post-processing of the inferred values of the radial positions

The inferred radial positions are updated based on the inferred angular positions. When deactivated, nodes with the same degree have the same radial position in the hyperbolic disk. Default is **`true`**.

```
# Command line
./mercator -k <edgelist_filename>

# Python module
mercator.embed(<edgelist_filename>, post_kappa=False)
```

#### Quiet mode

The program does not output details about the progress of the embedding procedure (i.e., the file `*.inf_log` is not generated). This mode supersedes the *verbose* mode (i.e., no output on screen). Default is **`false`**.

```
# Command line
./mercator -q <edgelist_filename>

# Python module
mercator.embed(<edgelist_filename>, quiet_mode=True)
```

#### Refine mode

When a file containing the previously inferred coordinates is provided (`*.inf_coord`), the program uses the already inferred positions and parameters as a starting point and perform another round of the likelihood maximization step to refine the inferred positions. The use of a different output filename is recommended to keep track of the changes. Default is **`false`**.

```
# Command line
./mercator -r <inferred_coordinates_filename> <edgelist_filename>

# Python module
mercator.embed(<edgelist_filename>, inf_coord=<inferred_coordinates_filename>)
```

#### Screen mode

The program outputs details about the progress of the embedding procedure on screen instead of generating the file `*.inf_log`. Default is **`false`**.

```
# Command line
./mercator -a <edgelist_filename>

# Python module
mercator.embed(<edgelist_filename>, screen_mode=True)
```

#### Validation mode

Validates and characterizes the inferred random network ensemble. This is done by generating a large number of networks based on the inferred parameters and positions. The following files are generated

* `*.inf_inf_pconn`: Compares the inferred probability of connection with the theoretical one based on the inferred parameters.
* `*.inf_theta_density`: Density the angular distribution
* `*.inf_vprop`: Contains the following topological properties for every nodes in the inferred ensemble and in the original network:
    * degree
    * sum of the degree of the neighbors
    * average degree of the neighbors
    * number of triangles
    * local clustering coefficient
* `*.inf_vstat`: Contains the following statistics of the inferred ensemble.
    * degree distribution
    * spectrum of the average degree of neighbors
    * clustering spectrum
* `*.obs_vstat`: Contains the same statistics as above but for the original network.

Default is **`false`**. A python script `plot_validation_of_embedding.py` is provided to visualize the results of the validation mode.

```
# Command line
./mercator -v <edgelist_filename>

# Python module
mercator.embed(<edgelist_filename>, validation_mode=True)

# To plot the various characterization quantities
python3 scripts/plot_validation_of_embedding.py <rootname of the edgelist (everything before the last ., if any)> <custom rootname (if any, otherwise leave blank)>
```





[pybind11]:  https://github.com/pybind/pybind11
