#!python

import numpy as np
import numba
import multiprocessing as mp











@numba.njit(cache=True, nogil=True)
def calculate_scores(
    queries,
    spec_indptr,
    low_limits,
    high_limits,
    peptide_count, #len(peptide_sequences),
    peptide_pointers,
):
    peptide_results = np.zeros(spec_indptr[-1], np.int64)
    score_results = np.zeros(spec_indptr[-1], np.float64)
    for spectrum_index in queries:
        spectrum_start = spec_indptr[spectrum_index]
        spectrum_end = spec_indptr[spectrum_index + 1]
        candidates = np.zeros(peptide_count, np.int64)
        for ion_index in range(spectrum_start, spectrum_end):
            peptide_low = low_limits[ion_index]
            peptide_high = high_limits[ion_index]
            if peptide_low == peptide_high:
                continue
            peptides = peptide_pointers[peptide_low: peptide_high]
            candidates[peptides] += 1
        for ion_index in range(spectrum_start, spectrum_end):
            peptide_low = low_limits[ion_index]
            peptide_high = high_limits[ion_index]
            if peptide_low == peptide_high:
                continue
            peptides = peptide_pointers[peptide_low: peptide_high]
            local_candicates = candidates[peptides]
            frequencies = np.bincount(local_candicates)
            frequencies = np.cumsum(frequencies[:0:-1])[::-1]
            for regression_index, value in enumerate(frequencies):
                if value == 1:
                    break
            else:
                continue
            if regression_index < 2:
                continue
            max_count = (len(frequencies) - 1)
            regression_constant = np.log(frequencies[0])
            regression_slope = (np.log(frequencies[regression_index]) - regression_constant) / regression_index
            score = regression_constant + regression_slope * max_count
            if score < 0:
                score_results[ion_index] = score
                peptide = peptides[local_candicates == max_count + 1][0]
                peptide_results[ion_index] = peptide
    return score_results, peptide_results





















# for f in /media/proteomics2/RAW/Synapt/HDMSe/Covid19_eSwab_Dilution/*blanco*; do apex3d -pRawDirName "$f" -outputDirName /home/sander/data/covid_eswab_hdmse/ -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 -bEnableCentroids 0; done



































#
#
# for i in range(60, 70):
#     in_name = f"/home/sander/Documents/Proteomics/data/ecoli/28Oct2016_0{i}.inet.csv"
#     out_name = f"/home/sander/Documents/Proteomics/data/ecoli_test/28Oct2016_0{i}.inet.csv"
#     df = pd.read_csv(in_name)
#     row_filter = (df["PRECURSOR_RT"] > 50) & (df["PRECURSOR_RT"] < 60)
#     # row_filter |= (df["PRECURSOR_MZ"] > 500) & (df["PRECURSOR_MZ"] < 600)
#     # row_filter |= (df["FRAGMENT_MZ"] > 500) & (df["FRAGMENT_MZ"] < 600)
#     # row_filter |= (df["FRAGMENT_LOGINT"] > 10) & (df["FRAGMENT_LOGINT"] < 20)
#     df[row_filter].to_csv(out_name, index=False)































#
#
# inet=inets[0]
# evi=evis[0]
#
# coords = inet.get_ion_coordinates(inet.dimensions)
# repro = evi.get_aligned_nodes_from_group()
# edges = inet.get_edges(data_as_index=True)
# positive_edge_evidence = evi.get_edge_mask_from_group()
# negative_edge_evidence = evi.get_edge_mask_from_group(
#     positive=False
# )
# print(inet.dimensions)
#


#
#
# node_selection = repro>=25
# node_selection &= (coords[3] > 10.0) & (coords[3] < 10.1)
# node_selection &= (coords[2] > 50) & (coords[2] < 60)
#
#
# selected_neighbors = edges[node_selection].T.tocsr()[node_selection]
# a, b = selected_neighbors.nonzero()
# positive_counts = positive_edge_evidence[
#     selected_neighbors.data
# ]
# negative_counts = negative_edge_evidence[
#     selected_neighbors.data
# ]
# edge_selection=...


#
# # edge_selection = (positive_counts >= self.min_positive_threshold)
# # edge_selection &= (positive_counts <= self.max_positive_threshold)
# # edge_selection &= (negative_counts >= self.min_negative_threshold)
# # edge_selection &= (negative_counts <= self.max_negative_threshold)
# a = a[edge_selection]
# b = b[edge_selection]
# start_edges = list(zip(coords[1][node_selection][a], coords[3][node_selection][a], coords[2][node_selection][a]))
# end_edges = list(zip(coords[1][node_selection][b], coords[3][node_selection][b], coords[2][node_selection][b]))
# colors = positive_counts[edge_selection] - negative_counts[edge_selection]
# plot_edges = np.array(list(zip(start_edges, end_edges)))
# color_order = np.argsort(colors)
# colors = colors[color_order]
# plot_edges = plot_edges[color_order]






# %matplotlib notebook
#
# import matplotlib.animation as animation
# from IPython.display import HTML
# from mpl_toolkits.mplot3d import axes3d
# from mpl_toolkits.mplot3d import art3d
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
#
# # fig = plt.figure(1, figsize=(13, 9))
# fig = plt.figure()
#
# ax = fig.add_subplot(111, projection='3d')
#
#
# def init_plot():
# #     edge_collection = ax.add_collection(
# #         matplotlib.collections.LineCollection([], [], [])
# #     )
#     ax.scatter(
#         coords[1][node_selection],
#         coords[3][node_selection],
#         coords[2][node_selection],
# #         marker=".",
#         c=coords[0][node_selection],
#         cmap="RdYlGn",
#         zorder=99999,
#         s=50
#     )
#     ax.set_xlabel(inet.dimensions[1])
#     ax.set_ylabel(inet.dimensions[3])
#     ax.set_zlabel(inet.dimensions[2])
#     edge_collection = art3d.Line3DCollection(plot_edges, alpha=0.2)
#     ax.add_collection(edge_collection)
#     ax.view_init(elev=15., azim=0)
#     edge_color_mapper = matplotlib.cm.ScalarMappable(
#         norm=matplotlib.colors.Normalize(
#             # vmin=-self.evidence.evidence_count,
#             vmin=-evi.evidence_count - 1,
#             # vmin=0,
#             vmax=evi.evidence_count
#         ),
#         cmap="RdYlGn"
#     )
#     edge_collection.set_color(edge_color_mapper.to_rgba(colors))
# #     ax.set_facecolor('black')
#     return fig,
#
# init_plot()
#
# def animate(i):
#     ax.view_init(elev=15., azim=3.6*i)
#     return fig,
#
# # Animate
# ani = animation.FuncAnimation(fig, animate, init_func=init_plot,
#                                frames=100, interval=100, blit=True)
#
# # # HTML(ani.to_html5_video())
# ani.save("/home/sander/test.mp4", dpi=300)
# print("done saving")
