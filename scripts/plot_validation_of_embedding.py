# -​*- coding: utf-8 -*​-
# @author: Antoine Allard <antoine.allard.1@gmail.com>

# Packages
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.stats as stt
import random as rnd
# import warnings
# warnings.filterwarnings("ignore")
import numpy as np
import sys
import os

# Global parameters for the figures.
plt.rcParams["text.usetex"] = True
plt.rcParams["font.size"] = 20
# plt.rcParams["font.family"] = "serif"
# plt.rcParams["font.serif"] = "Charter"
plt.rcParams["xtick.major.width"] = 1
plt.rcParams["xtick.major.size"] = 8
plt.rcParams["ytick.major.width"] = 1
plt.rcParams["ytick.major.size"] = 8





# ==================================================================================================
# Probability of connection (using .inf_pconn file)
# ==================================================================================================
def plot_connection_probability():
  global id_fig
  # Creates the subfigure.
  id_fig += 1
  ax = fig.add_subplot(nb_rows, nb_columns, id_fig)
  # Loads the data.
  pconn = np.loadtxt(pconn_inferred)
  # Plots the data.
  ax.plot(pconn[:, 0], pconn[:, 2],
          color="#555555", linewidth=1, linestyle="--", marker="None",
          markersize=5, markeredgecolor="k", markeredgewidth=1,
          markerfacecolor="kxkcd:light teal",
          label=r"inferred")
  ax.plot(pconn[:, 0], pconn[:, 1],
          color="#555555", linewidth=1, linestyle="None", marker="o",
          markersize=5, markeredgecolor="k", markeredgewidth=1,
          markerfacecolor="k",
          label=r"original")
  # Organizes the plot.
  ax.set_xscale("log")
  ax.set_yscale("log")
  xlim_left   = 10**(np.floor(np.log10(min([x for x,y in zip(pconn[:, 0], pconn[:, 1]) if x > 0 and y > 0]))))
  xlim_right  = 10**(np.ceil (np.log10(max([x for x,y in zip(pconn[:, 0], pconn[:, 1]) if x > 0 and y > 0]))))
  ylim_top    = 2
  ylim_bottom = 10**(np.floor(np.log10(min([y for y in pconn[:, 1] if y > 0]))))
  ax.set_xlim(left=xlim_left, right=xlim_right)
  ax.set_ylim(bottom=ylim_bottom, top=ylim_top)
  ax.set_xlabel(r"rescaled distance $(\chi)$")
  ax.set_ylabel(r"probability of connection")
  # Legend.
  lines1, labels1 = ax.get_legend_handles_labels()
  ax.legend(lines1, labels1,
            shadow=False, fancybox=False, prop={"size": 20},
            frameon=False, numpoints=1,
            loc="lower left", ncol=1)





# ==================================================================================================
# Probability of connection (using raw data)
# ==================================================================================================
def plot_connection_probability_using_raw_data():
  global id_fig
  # Creates the subfigure.
  id_fig += 1
  ax = fig.add_subplot(nb_rows, nb_columns, id_fig)
  # Builds the adjacency list
  el = np.loadtxt(ename_original, dtype="str")
  adjacency_list = {}
  for i in range(el.shape[0]):
    v1 = el[i, 0]
    v2 = el[i, 1]
    if v1 not in adjacency_list:
      adjacency_list[v1] = []
    if v2 not in adjacency_list:
      adjacency_list[v2] = []
    adjacency_list[v1].append(v2)
    adjacency_list[v2].append(v1)
  # Gets the values of beta, mu and R.
  with open(coord_inferred) as f:
    for i in range(8):
      sline=next(f).strip().split()
    beta=float(next(f).strip().split()[3])
    mu=float(next(f).strip().split()[3])
    R=float(next(f).strip().split()[3])
  # Builds the connection probability.
  is_edge = []
  distance = []
  name = {}
  names = np.loadtxt(coord_inferred, usecols=[0], dtype="str")
  for i in range(len(names)):
    name[names[i]] = i
  kappa = np.loadtxt(coord_inferred, usecols=[1])
  theta = np.loadtxt(coord_inferred, usecols=[2])
  for v1 in adjacency_list:
    for v2 in adjacency_list:
      if v1 != v2:
        # Marks whether the two vertices are connected.
        if v2 in adjacency_list[v1]:
          is_edge.append(1)
        else:
          is_edge.append(0)
        # Computes the distances between the two vertices.
        i1 = int(name[v1])
        i2 = int(name[v2])
        da = np.pi - np.fabs( np.pi - np.fabs( theta[i1] - theta[i2] ) )
        distance.append( (R * da) / (mu * kappa[i1] * kappa[i2]) )
  x_min = 10**(np.floor(np.log10(min([i for i in distance if i > 0]))))
  x_max = 10**(np.ceil(np.log10(max(distance))))
  bins = np.logspace(np.log10(x_min), np.log10(x_max), 100)
  x_his = np.histogram(distance, bins, weights=distance, density=False)[0]
  y_his = np.histogram(distance, bins, weights=is_edge,  density=False)[0]
  n_his = np.histogram(distance, bins,                   density=False)[0]
  x_his = [xx / nn for xx, nn in zip(x_his, n_his) if nn > 0]
  y_his = [yy / nn for yy, nn in zip(y_his, n_his) if nn > 0]
  # Plots the data.
  ax.plot(bins, 1 / (1 + bins**beta),
          color="#555555", linewidth=1, linestyle="--", marker="None",
          markersize=5, markeredgecolor="k", markeredgewidth=1,
          markerfacecolor="k",
          label=r"inferred")
  ax.plot(x_his, y_his,
          color="#555555", linewidth=1, linestyle="None", marker="o",
          markersize=5, markeredgecolor="k", markeredgewidth=1,
          markerfacecolor="k",
          label=r"original")
  # Organizes the plot.
  ax.set_xscale("log")
  ax.set_yscale("log")
  xlim_left   = 10**(np.floor(np.log10(min([x for x,y in zip(x_his, y_his) if y > 0]))))
  xlim_right  = 10**(np.ceil (np.log10(max([x for x,y in zip(x_his, y_his) if y > 0]))))
  ylim_top    = 2
  ylim_bottom = 10**(np.floor(np.log10(min([y for y in y_his if y > 0]))))
  ax.set_xlim(left=xlim_left, right=xlim_right)
  ax.set_ylim(bottom=ylim_bottom, top=ylim_top)
  ax.set_xlabel(r"rescaled distance $(\chi)$")
  ax.set_ylabel(r"probability of connection")
  # Legend.
  lines1, labels1 = ax.get_legend_handles_labels()
  ax.legend(lines1, labels1,
            shadow=False, fancybox=False, prop={"size": 20},
            frameon=False, numpoints=1,
            loc="lower left", ncol=1)





# ==================================================================================================
# Embedding in the hyperbolic disk (plots edges by group to accelerate the process)
# ==================================================================================================
def plot_hyperbolic_disk(nb_seg=100):
  global id_fig
  # Creates the subfigure.
  id_fig += 1
  ax = fig.add_subplot(nb_rows, nb_columns, id_fig, projection='polar')
  # Prepares the data prior to plot.
  names  = np.loadtxt(coord_inferred, usecols=[0], dtype="str")
  name2id = {}
  for i in range(len(names)):
    name2id[names[i]] = int(i)
  theta  = np.loadtxt(coord_inferred, usecols=[2])
  radius = np.loadtxt(coord_inferred, usecols=[3])
  edges  = np.loadtxt(ename_original, dtype="str")
  # Plots the edges by group of seg_size (to accelerate run time).
  nb_edges = edges.shape[0]
  rnd_id = np.arange(nb_edges)
  rnd.shuffle(rnd_id)
  seg_size = np.ceil(nb_edges / nb_seg)
  for e_i in np.arange(0, nb_edges, seg_size):
    e_f = (e_i + seg_size) if (e_i + seg_size) < nb_edges else nb_edges
    t_seg = []
    r_seg = []
    for i in np.arange(e_i, e_f, 1):
      i1 = name2id[edges[rnd_id[int(i)], 0]]
      i2 = name2id[edges[rnd_id[int(i)], 1]]
      t_seg.extend([theta[i1], theta[i2], np.nan])
      r_seg.extend([radius[i1], radius[i2], np.nan])
    ax.plot(t_seg[:-1], r_seg[:-1],
            color="#555555", linewidth=0.5, linestyle="-", marker="None",
            markersize=4, markeredgecolor="None", markeredgewidth=None,
            markerfacecolor="b", alpha=0.1)
  # Plots the vertices.
  ax.plot(theta, radius,
              color="#333333", linewidth=0.5, linestyle="None", marker="o",
              markersize=3, markeredgecolor="None", markeredgewidth=None,
              markerfacecolor="b", alpha=1.0)
  # Removes labels.
  ax.set_rticks([])
  ax.set_xticklabels("")





# ==================================================================================================
# Embedding in the hyperbolic disk (plots every individual edges to maximizes the additive transparency, slower)
# ==================================================================================================
def plot_hyperbolic_disk_slow():
  global id_fig
  # Creates the subfigure.
  id_fig += 1
  ax = fig.add_subplot(nb_rows, nb_columns, id_fig, projection='polar')
  # Prepares the data prior to plot.
  names  = np.loadtxt(coord_inferred, usecols=[0], dtype="str")
  name2id = {}
  for i in range(len(names)):
    name2id[names[i]] = int(i)
  theta  = np.loadtxt(coord_inferred, usecols=[2])
  radius = np.loadtxt(coord_inferred, usecols=[3])
  edges  = np.loadtxt(ename_original, dtype="str")
  # Plots the edges individually (to maximize the additive transparency).
  for edge in edges:
    i1 = name2id[edge[0]]
    i2 = name2id[edge[1]]
    ax.plot([theta[i1], theta[i2]], [radius[i1], radius[i2]],
             color="#555555", linewidth=0.5, linestyle="-", marker="None",
             markersize=4, markeredgecolor="None", markeredgewidth=None,
             markerfacecolor="b", alpha=0.1, zorder=-1)
  # Plots the vertices.
  ax.plot(theta, radius,
              color="#333333", linewidth=0.5, linestyle="None", marker="o",
              markersize=3, markeredgecolor="None", markeredgewidth=None,
              markerfacecolor="b", alpha=1.0)
  # Removes labels.
  ax.set_rticks([])
  ax.set_xticklabels("")





# ==================================================================================================
# Theta density (using raw data)
# ==================================================================================================
def plot_theta_density():
  global id_fig
  # Creates the subfigure.
  id_fig += 1
  ax = fig.add_subplot(nb_rows, nb_columns, id_fig)
  # Prepares the data prior to plot.
  theta = np.loadtxt(theta_inferred)
  # Plots the data.
  ax.plot(theta[:, 0], theta[:, 2], c="#555555", marker="None", ls="--", lw=2, label=r"uniform")
  ax.bar(theta[:, 0], theta[:, 1], align="center", width=(2 * np.pi / len(theta)))
  # Organizes the plot.
  ax.set_ylim(bottom=0, top=(2 / (len(theta[:, 0]) - 1)))
  ax.set_xlabel(r"inferred angle $(\theta)$")
  ax.set_ylabel(r"density")
  # Legend.
  lines1, labels1 = ax.get_legend_handles_labels()
  ax.legend(lines1, labels1,
            shadow=False, fancybox=False, prop={"size": 28},
            # shadow=False, fancybox=False, prop={"size": 20},
            frameon=False, numpoints=1,
            loc="upper left", ncol=1)





# ==================================================================================================
# Theta density (using raw data)
# ==================================================================================================
def plot_theta_density_using_raw_data():
  global id_fig
  # Creates the subfigure.
  id_fig += 1
  ax = fig.add_subplot(nb_rows, nb_columns, id_fig)
  # Prepares the data prior to plot.
  theta = np.loadtxt(coord_inferred, usecols=[2])
  # Plots the data.
  bins = np.linspace(0.0, 2 * np.pi, 26 )
  # bins = np.linspace(0.0, 2 * np.pi, np.min([int(len(theta) / 10), 25]) )
  x_his = np.histogram(theta, bins, weights=theta, density=False)[0]
  n_his = np.histogram(theta, bins, density=False)[0]
  bins  = [bb for bb, nn in zip(bins[:-1], n_his) if nn > 0]
  x_his = [xx / nn for xx, nn in zip(x_his, n_his) if nn > 0]
  n_his = [nn / float(np.sum(n_his)) for nn in n_his if nn > 0]
  uniform = 1 / (26 - 1)
  # uniform = 1 / (np.min([int(len(theta) / 10), 25]) - 1)
  ax.plot([0, 2 * np.pi], uniform * np.ones(2), c="#555555", marker="None", ls="--", lw=2, label=r"uniform")
  ax.bar(bins, n_his, align="edge", width=(2 * np.pi / np.min([int(len(theta) / 10), 25]) ))
  # Organizes the plot.
  ax.set_ylim(bottom=0, top=2 / (len(bins) - 1))
  ax.set_xlabel(r"inferred angle $(\theta)$")
  ax.set_ylabel(r"density")
  # Legend.
  lines1, labels1 = ax.get_legend_handles_labels()
  ax.legend(lines1, labels1,
            shadow=False, fancybox=False, prop={"size": 20},
            frameon=False, numpoints=1,
            loc="upper left", ncol=1)





# ==================================================================================================
# Inferred theta vs original theta
# ==================================================================================================
def plot_inferred_theta_vs_original_theta():
  global id_fig
  # Creates the subfigure.
  id_fig += 1
  ax = fig.add_subplot(nb_rows, nb_columns, id_fig)
  # Loads the list of names.
  names = np.loadtxt(coord_inferred, usecols=[0], dtype='str')
  # Sets the name2num dictionary and the number of vertices.
  name2num = {}
  nb_vertices = 0
  for n in names:
    name2num[n] = nb_vertices
    nb_vertices += 1
  # Extracts the inferred values of theta.
  theta_inferred = np.zeros(nb_vertices)
  names = np.loadtxt(coord_inferred, usecols=[0], dtype='str')
  values = np.loadtxt(coord_inferred, usecols=[2])
  for i in range(len(names)):
      theta_inferred[name2num[names[i]]] = values[i]
  # Extracts the original values of theta.
  theta_original = np.zeros(nb_vertices)
  names = np.loadtxt(coord_original, usecols=[0], dtype='str')
  values = np.loadtxt(coord_original, usecols=[2])
  for i in range(len(names)):
    if names[i] in name2num:
      theta_original[name2num[names[i]]] = values[i]
  # Plots the values of both thetas
  ax.plot(np.arange(0, 2 * np.pi), np.arange(0, 2 * np.pi), c="#555555", marker="None", ls="--", lw=2)
  ax.plot(theta_original, theta_inferred, c="b", marker="o", ms=3, ls="None", lw=2, mec="None")
  # # Prints quantification of the agreement.
  # ax.annotate(r"$\rho = $ {:0.2f}".format(stt.pearsonr(theta_original, theta_inferred)[0]),
  #             xy=(0.025, 0.975),  xycoords='axes fraction',
  #             horizontalalignment='left', verticalalignment='top')
  # Bounds.
  ax.set_xbound(0, 2 * np.pi)
  ax.set_ybound(0, 2 * np.pi)
  # Axes label.
  ax.set_xlabel(r"original angle $(\theta)$")
  ax.set_ylabel(r"inferred angle $(\theta)$")





# ==================================================================================================
# Inferred vertices property obtained by simulations.
# ==================================================================================================
def plot_comparison_vprop_using_cols(cols):
  global id_fig
  # Creates the subfigure.
  id_fig += 1
  ax = fig.add_subplot(nb_rows, nb_columns, id_fig)
  # Loads the list of names.
  data = np.loadtxt(vprop_inferred, usecols=cols)
  # Plots the values of both thetas
  min_bound = 2 * np.min( [ np.min([v for v in data[:, 0] if v > 0]), np.min([v for v in data[:, 0] if v > 0]) ] ) / 3
  max_bound = 3 * np.max( [np.max(data[:, 0]), np.max(data[:, 1])] ) / 2
  ax.plot([min_bound, max_bound], [min_bound, max_bound], c="#555555", marker="None", ls="--", lw=2)
  ax.errorbar(data[:, 0], data[:, 1], yerr=(2 * data[:, 2]), alpha=0.2, #fmt="None",
              c="b", ls="None", lw=1,
              marker="o", ms=4, mfc="b", mec="None",
              ecolor="#AAAAAA", elinewidth=1, capsize=4, capthick=1)
  # ax.plot(data[:, 0], data[:, 1], c="b", marker="o", ms=4, ls="None", lw=2, mec="None", alpha=0.2)
  # Prints quantification of the agreement.
  ax.annotate(r"$\rho = {:0.3f}$".format(stt.pearsonr(data[:, 0], data[:, 1])[0]),
              xy=(0.025, 0.95),  xycoords='axes fraction',
              horizontalalignment='left', verticalalignment='center')
  # Computes the chi-squared normalized score.
  chi_squared = 0;
  nb_vertices = len(data[:, 0])
  for i in range(nb_vertices):
    if data[i, 2] > 0:
      chi_squared += ( (data[i, 0] - data[i, 1]) / data[i, 2] )**2
  ax.annotate(r"$\chi^2 / N = {:0.3f}$".format(chi_squared / nb_vertices),
  # ax.annotate(r"$\chi^2 = $ {:0.2f}".format( (chi_squared - nb_vertices) / np.sqrt(2 * nb_vertices) ),
  # ax.annotate(r"$\chi^2 = $ {:0.2f}".format( (chi_squared - nb_vertices) ),
              xy=(0.025, 0.875),  xycoords='axes fraction',
              horizontalalignment='left', verticalalignment='center')
  # Fraction of vertices for which the property lies outside 1 sigma of the ensemble.
  number_of_outsiders = 0
  for i in range(nb_vertices):
    if np.fabs(data[i, 0] - data[i, 1]) > 1.5 * data[i, 2]:
      number_of_outsiders += 1
  ax.annotate(r"$\zeta = {:0.3f}$".format(number_of_outsiders / nb_vertices),
              xy=(0.025, 0.80),  xycoords='axes fraction',
              horizontalalignment='left', verticalalignment='center')
  # Axes scales.
  ax.set_xscale("log")
  ax.set_yscale("log")
  # Bounds.
  ax.set_xbound(min_bound, max_bound)
  ax.set_ybound(min_bound, max_bound)
  # Returns the axis object for individual adjustments.
  return ax





# ==================================================================================================
# Degree distribution
# ==================================================================================================
def plot_comparison_vstat_using_cols(cols_ori, cols_inf, pos="lower right"):
  global id_fig
  # Creates the subfigure.
  id_fig += 1
  ax = fig.add_subplot(nb_rows, nb_columns, id_fig)
  # Loads and plots the original degree distribution.
  data = np.loadtxt(vstat_inferred, usecols=cols_inf)
  ax.plot(data[:, 0], data[:, 1], c="#222222", marker="None", ls="-", lw=1, ms=1, markerfacecolor="r", label=r"inferred")
  ax.fill_between(data[:, 0], data[:, 1]-(2 * data[:, 2]), data[:, 1]+(2 * data[:, 2]), color="#CCCCCC")
  # ax.plot(data[:, 0], data[:, 1]-(2 * data[:, 2]), c="k")
  # ax.plot(data[:, 0], data[:, 1]+(2 * data[:, 2]), c="k")
  # Loads and plots the original degree distribution.
  data = np.loadtxt(vstat_original, usecols=cols_ori)
  ax.plot(data[:, 0], data[:, 1], c="b", marker="s", ls="-", lw=1, markerfacecolor="b", label=r"original")
  # Legend.
  lines1, labels1 = ax.get_legend_handles_labels()
  ax.legend(lines1, labels1,
            shadow=False, fancybox=False, prop={"size": 20},
            frameon=False, numpoints=1,
            loc=pos, ncol=1)
  return ax





# ==================================================================================================
# Inferred vertices properties obtained by simulations.
# ==================================================================================================
def plot_comparison_vprop():
  # Degrees.
  ax = plot_comparison_vprop_using_cols([1, 2, 3])
  ax.set_xlabel(r"original degree $(k)$")
  ax.set_ylabel(r"inferred degree $(k)$")
  # Sum of the degree of neighbors.
  ax = plot_comparison_vprop_using_cols([4, 5, 6])
  ax.set_xlabel(r"original sum degree of neighbors")
  ax.set_ylabel(r"inferred sum degree of neighbors")
  # Average degree of neighbors.
  ax = plot_comparison_vprop_using_cols([7, 8, 9])
  ax.set_xlabel(r"original average degree of neighbors $(k_\mathrm{nn})$")
  ax.set_ylabel(r"inferred average degree of neighbors $(k_\mathrm{nn})$")
  # Sum of triangles to which each vertex participates.
  ax = plot_comparison_vprop_using_cols([10, 11, 12])
  ax.set_xlabel(r"original number of triangles")
  ax.set_ylabel(r"inferred number of triangles")
  # Clustering.
  ax = plot_comparison_vprop_using_cols([13, 14, 15])
  ax.set_xlabel(r"original clustering $(c)$")
  ax.set_ylabel(r"inferred clustering $(c)$")
  # Degree distribution.
  ax = plot_comparison_vstat_using_cols([0, 1], [0, 1, 2], "upper right")
  ax.set_xscale("log")
  ax.set_yscale("log")
  ax.set_xlabel(r"degree $(k)$")
  ax.set_xlabel(r"degree")
  # ax.set_ylabel(r"degree distribution $[p(k)]$"))
  ax.set_ylabel(r"degree distribution")
  # Sum degree of neighbors spectrum.
  ax = plot_comparison_vstat_using_cols([0, 2], [0, 3, 4], "lower right")
  ax.set_xlabel(r"degree")
  # ax.set_xlabel(r"degree $(k)$"))
  ax.set_ylabel(r"sum degree of neighbors")
  # Average degree of neighbors spectrum.
  ax = plot_comparison_vstat_using_cols([0, 3], [0, 5, 6], "lower right")
  ax.set_xlabel(r"degree")
  # ax.set_xlabel(r"degree $(k)$"))
  # ax.set_ylabel(r"average degree of neighbors $[\bar{k}_\mathrm{nn}(k)]$"))
  ax.set_ylabel(r"average degree of neighbors")
  # Number of triangles.
  ax = plot_comparison_vstat_using_cols([0, 4], [0, 7, 8], "lower right")
  # ax.set_xlabel(r"degree $(k)$"))
  ax.set_xlabel(r"degree")
  ax.set_ylabel(r"number of triangles")
  # Clustering spectrum.
  ax = plot_comparison_vstat_using_cols([0, 5], [0, 9, 10], "upper right")
  # ax.set_xlabel(r"degree $(k)$")
  ax.set_xlabel(r"degree")
  # ax.set_ylabel(r"clustering coefficient $[\bar{c}(k)]$"))
  ax.set_ylabel(r"clustering coefficient")





# ==================================================================================================
# ==================================================================================================
if __name__ == "__main__":

  # Sets the rootnames of the files.
  if len(sys.argv) == 2:
    name1 = sys.argv[1]
    name2 = sys.argv[1]
  elif len(sys.argv) == 3:
    name1 = sys.argv[1]
    name2 = sys.argv[2]
  else:
    print("WARNING: Bad arguments. Please refer to the script for details. Exiting...")
    exit()

  # Sets the filenames.
  ename_original = name1 + ".edge"              # edgelist
  coord_original = name1 + ".gen_coord"         # original coordinates
  vstat_original = name2 + ".obs_vstat"         # original properties
  coord_inferred = name2 + ".inf_coord"         # inferred coordinates
  vprop_inferred = name2 + ".inf_vprop"         # inferred properties
  vstat_inferred = name2 + ".inf_vstat"         # inferred properties
  pconn_inferred = name2 + ".inf_pconn"         # inferred probability of connection
  theta_inferred = name2 + ".inf_theta_density" # inferred angular density
  fname_outvalid = name2 + "_validation.pdf"    # figure's filename
  figure_title   = ('%r'%os.path.basename(name2)).replace('_', '\_').strip("'")

  # Sets the figure objects.
  id_fig, nb_rows, nb_columns = 0, 2, 2
  if os.path.isfile(vprop_inferred):
    id_fig, nb_rows, nb_columns = 0, 5, 3
    # id_fig, nb_rows, nb_columns = 0, 2, 3
  fig = plt.figure(figsize=(8 * nb_columns, 6 * nb_rows))

  # Plots the probablity of connection.
  if os.path.isfile(pconn_inferred):
    plot_connection_probability()
  else:
    plot_connection_probability_using_raw_data()

  # Plots the graph in the hyperbolic disk.
  plot_hyperbolic_disk()
  # plot_hyperbolic_disk_slow()

  # Plots the angular density.
  if os.path.isfile(theta_inferred):
    plot_theta_density()
  else:
    plot_theta_density_using_raw_data()

  # Plots
  if os.path.isfile(vprop_inferred):
    plot_comparison_vprop()


  # Plots the comparison between the inferred and original angles.
  if os.path.isfile(coord_original):
    plot_inferred_theta_vs_original_theta()

  # Adds figure title.
  plt.suptitle(figure_title, y=0.99, va="top", ha="center")

  # Save to file.
  plt.tight_layout(rect=[-0.012, -0.012, 1.012, 0.988])
  fig.savefig(fname_outvalid)
