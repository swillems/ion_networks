#!python

import numpy as np
import numba
import multiprocessing as mp










def get_mzs_and_frequencies(mzs, ppm):
    new_mzs = np.sort(mzs)
    lower = np.searchsorted(new_mzs, new_mzs * (1 - ppm / 10**6), "left")
    upper = np.searchsorted(new_mzs, new_mzs * (1 + ppm / 10**6), "right")
    return new_mzs, upper - lower


@numba.njit
def find_peaks(intensities, mzs, max_distance=0.9):
    peaks = np.zeros(int(mzs[-1]), np.int64)
    current_max_mz = 0
    current_max_int = 0
    current_max_index = 0
    for index, (intensity, mz) in enumerate(zip(intensities, mzs)):
        if mz > current_max_mz + max_distance:
            peaks[int(current_max_mz)] = current_max_index
            current_max_mz = mz
            current_max_int = intensity
            current_max_index = index
        elif intensity > current_max_int:
            current_max_mz = mz
            current_max_int = intensity
            current_max_index = index
    return peaks


















@numba.njit
def determine_isotopic_rt_difference(
    mzs,
    rts,
    isotopic_distance,
    ppm,
    best_hit_only=False
):
    mz_order = np.argsort(mzs)
    mzs_in_mz_order = mzs[mz_order]
    rts_in_mz_order = rts[mz_order]
    if isotopic_distance > 0:
        lower_limits = np.searchsorted(
            mzs_in_mz_order,
            (mzs_in_mz_order + isotopic_distance) / (1 + ppm * 10**-6),
            "left"
        )
    else:
        lower_limits = np.arange(len(mzs)) + 1
    upper_limits = np.searchsorted(
        mzs_in_mz_order,
        (mzs_in_mz_order + isotopic_distance) * (1 + ppm * 10**-6),
        "right"
    )
    if not best_hit_only:
        first_rts = np.repeat(rts_in_mz_order, upper_limits - lower_limits)
        second_rts = np.concatenate(
            [
                rts_in_mz_order[l: h] for l, h in zip(
                    lower_limits,
                    upper_limits
                )
            ]
        )
    else:
        isotopes_found = upper_limits > lower_limits
        first_rts = rts_in_mz_order[isotopes_found]
        lower_limits = lower_limits[isotopes_found]
        upper_limits = upper_limits[isotopes_found]
        second_rts = np.array(
            [
                np.min(rts_in_mz_order[l: h]) for l, h in zip(
                    lower_limits,
                    upper_limits
                )
            ]
        )
    isotope_rt_diffs = np.abs(second_rts - first_rts)
    isotope_rt_diffs, isotope_count = np.unique(
        isotope_rt_diffs,
        return_counts=True
    )
    isotope_count = np.cumsum(isotope_count)
    return isotope_rt_diffs, isotope_count


def matrix_multiplication(a, b):
    aa_indices, aa_values = a.nonzero()
    bb_values, bb_indices = b.nonzero()
    return matrix_multiplication_numba_wrapper(
        aa_indices,
        bb_indices,
        aa_values,
        bb_values,
        np.argsort(aa_values),
        np.argsort(bb_values),
    )


@numba.njit(
    numba.i4[:, :](
        numba.i4[:],
        numba.i4[:],
        numba.i4[:],
        numba.i4[:],
        numba.i8[:],
        numba.i8[:]
    )
)
def matrix_multiplication_numba_wrapper(
    aa_indices,
    bb_indices,
    aa_values,
    bb_values,
    aa_order,
    bb_order
):
    aa_indices = aa_indices[aa_order]
    aa_values = aa_values[aa_order]
    bb_indices = bb_indices[bb_order]
    bb_values = bb_values[bb_order]
    low_limits = np.searchsorted(
        aa_values,
        bb_values,
        "left"
    )
    high_limits = np.searchsorted(
        aa_values,
        bb_values,
        "right"
    )
    diffs = high_limits - low_limits
    ends = np.cumsum(diffs)
    result = np.empty((2, ends[-1]), np.int32)
    result[1, :] = np.repeat(
        bb_indices,
        high_limits - low_limits
    )
    for l, h, e, d in zip(low_limits, high_limits, ends, diffs):
        result[0, e - d: e] = aa_indices[l: h]
    return result


def parallelizedGenerator(
    process_count,
    function,
    iterable,
    *args,
    **kwargs
):
    in_queue = mp.Queue()
    out_queue = mp.Queue()

    def _parallel_function():
        while True:
            element = in_queue.get()
            if element is None:
                out_queue.put(None)
                return
            result = function(
                element,
                *args,
                **kwargs
            )
            out_queue.put(result)

    for element in iterable:
        in_queue.put(element)
    processes = []
    for process_index in range(process_count):
        process = mp.Process(target=_parallel_function,)
        in_queue.put(None)
        process.start()
        processes.append(process)
    running_processes = process_count
    while running_processes != 0:
        result = out_queue.get()
        if result is None:
            running_processes -= 1
        else:
            yield result
    in_queue.close()
    out_queue.close()
    in_queue.join_thread()
    out_queue.join_thread()
    while processes:
        process = processes.pop()
        process.join()


import functools
import time
import multiprocessing as mp


def generate_in_parallel(cpu_count=0):
    max_cpu_count = mp.cpu_count()
    if cpu_count > max_cpu_count:
        cpu_count = max_cpu_count
    elif cpu_count <= 0:
        if cpu_count + max_cpu_count <= 0:
            cpu_count = 1
        else:
            cpu_count += max_cpu_count
    def decorator_repeat(func):
        @functools.wraps(func)
        def wrapper_generate_in_parallel(*args, **kwargs):
            in_queue = mp.Queue()
            out_queue = mp.Queue()
            def _parallel_function():
                while True:
                    element = in_queue.get()
                    if element is None:
                        out_queue.put(None)
                        return
                    result = func(
                        element,
                        *args[1:],
                        **kwargs
                    )
                    out_queue.put(result)
            for element in args[0]:
                in_queue.put(element)
            processes = []
            for process_index in range(cpu_count):
                process = mp.Process(
                    target=_parallel_function,
                )
                in_queue.put(None)
                process.start()
                processes.append(process)
            running_processes = cpu_count
            while running_processes != 0:
                result = out_queue.get()
                if result is None:
                    running_processes -= 1
                else:
                    yield result
            in_queue.close()
            out_queue.close()
            in_queue.join_thread()
            out_queue.join_thread()
            while processes:
                process = processes.pop()
                process.join()
        return wrapper_generate_in_parallel
    return decorator_repeat


























def calculate(
    inet,
    to_select_per_sample,
    ppm = 10,
    calibration_prior_estimates = {
        "PRECURSOR_RT": 0.5,
        "PRECURSOR_MZ": 5,
    },
    calibration_mzs = {
        "target": 113.084064,
        "decoy": 2*113.084064,
    },
    min_window = 100,
    usable = (0.1, 0.5),
):
    calibration_diffs = {
        "target": {},
        "decoy": {},
    }
    random_diffs = {
        "target": {},
        "decoy": {},
    }
    estimations = {}
    fints, fmzs = inet.get_ion_coordinates(
        ["FRAGMENT_LOGINT", "FRAGMENT_MZ"]
    )
    result_coordinates = {}
    precursor_arrays = {
        dimension: inet.get_ion_coordinates(dimension) for dimension in inet.precursor_dimensions
    }
    c = np.argpartition(fints, - to_select_per_sample)[-to_select_per_sample:]
    c = c[np.argsort(fmzs[c])]
    mzs = fmzs[c]
    for selection, mz_distance in calibration_mzs.items():
        if mz_distance > 0:
            lower_limits = np.searchsorted(
                mzs,
                (mzs + mz_distance) / (1 + ppm * 10**-6),
                "left"
            )
        else:
            lower_limits = np.arange(len(mzs)) + 1
        upper_limits = np.searchsorted(
            mzs,
            (mzs + mz_distance) * (1 + ppm * 10**-6),
            "right"
        )
        indptr = np.zeros(inet.node_count + 1, np.int64)
        indptr[c + 1] = upper_limits - lower_limits
        indptr = np.cumsum(indptr)
        order = np.argsort(c)
        indices = np.concatenate(
            [
                c[low: high] for low, high in zip(lower_limits[order], upper_limits[order])
            ]
        )
        for dimension, prior_estimate in calibration_prior_estimates.items():
            calibration_diffs[selection][dimension] = np.abs(
                np.repeat(
                    precursor_arrays[dimension],
                    np.diff(indptr)
                ) - precursor_arrays[dimension][indices]
            )
            random_diffs[selection][dimension] = np.sum(
                calibration_diffs[selection][dimension] > prior_estimate
            )
    for dimension in precursor_arrays:
        x = np.concatenate(
            [
                calibration_diffs["target"][dimension],
                calibration_diffs["decoy"][dimension],
            ]
        )
        y = np.repeat(
            [
                1,
                -random_diffs["target"][dimension] / random_diffs["decoy"][dimension]
            ],
            [
                calibration_diffs["target"][dimension].shape[0],
                calibration_diffs["decoy"][dimension].shape[0]
            ]
        )
        o = np.argsort(x)
        x = x[o]
        y = np.cumsum(y[o])
        result_coordinates[dimension] = (x, y)
    for dimension in precursor_arrays:
        x=np.sort(calibration_diffs["target"][dimension])
        use_slice = slice(
            int(x.shape[0]*usable[0]),
            int(x.shape[0]*usable[1])
        )
        # offset = np.max(
        #     x[min_window:]-x[:-min_window]
        # )
        offset = np.max(
            x[use_slice][min_window:]-x[use_slice][:-min_window]
        )
        # offset=0.001
#         print(offset)
        l = np.searchsorted(x, x - offset)
        r = np.searchsorted(x, x)
        y = r-l
        x, y = x[l>0],y[l>0]
        # plt.plot(x,np.cumsum(y)/np.arange(y.shape[0]))
#         plt.plot(x,y)
        ransac = sklearn.linear_model.RANSACRegressor()
        ransac.fit(
            x[use_slice].reshape(-1, 1),
            y[use_slice].reshape(-1, 1),
        )
#         plt.plot(
#             x[[0,-1]],
#             ransac.predict(x[[0,-1]].reshape(-1, 1))
#         )
#         plt.plot(
#             x[use_slice],
#             y[use_slice]
#         )
        min_index = np.argmin(
            ransac.predict(x.reshape(-1, 1)).flatten() < y
        )
#         plt.title(f"{dimension}: {x[min_index]}")
        estimations[dimension] = x[min_index]
    return result_coordinates, calibration_diffs, estimations





















# for f in /media/proteomics2/RAW/Synapt/HDMSe/Covid19_eSwab_Dilution/*blanco*; do apex3d -pRawDirName "$f" -outputDirName /home/sander/data/covid_eswab_hdmse/ -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 -bEnableCentroids 0; done



























#!python

# builtin
import os
# external
import numpy as np
import matplotlib.pyplot as plt
import PySimpleGUI as sg
import matplotlib
matplotlib.use('TkAgg')
# local
import ms_utils


BROWSER_IMAGES_PATH = os.path.join(
    ms_utils.LIB_PATH,
    "browser_images"
)
LOADING_ANIMATION_FILE_NAME = os.path.join(
    BROWSER_IMAGES_PATH,
    "load_animation.gif"
)


class Browser(object):
    # TODO: Docstring

    def __init__(self, evidence):
        # TODO: Docstring
        self.evidence = evidence
        self.ion_network = self.evidence.ion_network
        sg.PopupAnimated(
            LOADING_ANIMATION_FILE_NAME,
            message='Loading',
            # time_between_frames=10
            # non_blocking=True
        )
        self.create_overview_window()
        self.create_plot_window()
        while True:
            if self.update_window(self.overview_window) is not None:
                break
            sg.PopupAnimated(None)
            # if self.update_window(self.plot_window) is not None:
            #     break
        self.overview_window.close()
        self.plot_window.close()

    def create_overview_window(self):
        # TODO: Docstring
        layout = []
        self.node_dimensions = {}
        for dimension in self.ion_network.dimensions:
            coordinates = self.ion_network.get_ion_coordinates(dimension)
            min_coordinate = float(np.min(coordinates))
            max_coordinate = float(np.max(coordinates))
            self.node_dimensions[dimension] = [min_coordinate, max_coordinate]
            row = [
                sg.Text(f"{dimension}", size=(25, 1)),
                sg.InputText(
                    min_coordinate,
                    key=f"min_{dimension}",
                    size=(10, 1),
                    # enable_events=True
                ),
                sg.InputText(
                    max_coordinate,
                    key=f"max_{dimension}",
                    size=(10, 1),
                    # enable_events=True
                ),
            ]
            layout.append(row)
        max_node_count = float(self.evidence.evidence_count)
        self.node_dimensions["node_evidence"] = [max_node_count, max_node_count]
        self.node_threshold = 2 * len(self.node_dimensions)
        self.nodes = np.repeat(self.node_threshold, self.ion_network.node_count)
        self.edges = self.ion_network.get_edges(data_as_index=True)
        self.positive_edge_evidence = self.evidence.get_edge_mask_from_group()
        self.negative_edge_evidence = self.evidence.get_edge_mask_from_group(
            positive=False
        )
        node_evidence = self.evidence.get_aligned_nodes_from_group()
        self.nodes -= node_evidence < max_node_count
        layout.append(
            [
                sg.Text('NODE EVIDENCE', size=(25, 1)),
                sg.InputText(
                    max_node_count,
                    key=f"min_node_count",
                    size=(10, 1),
                    # enable_events=True
                ),
                sg.InputText(
                    max_node_count,
                    key=f"max_node_count",
                    size=(10, 1),
                    # enable_events=True
                ),
            ]
        )
        self.x_axis = "FRAGMENT_MZ"
        layout.append(
            [
                sg.Text('x-axis', size=(25, 1)),
                sg.Combo(
                    self.ion_network.dimensions,# + ["node_count"],
                    size=(21, 1),
                    default_value=self.x_axis,
                    key="x_axis",
                    # enable_events=True
                )
            ]
        )
        self.y_axis = "FRAGMENT_LOGINT"
        layout.append(
            [
                sg.Text('y-axis', size=(25, 1)),
                sg.Combo(
                    self.ion_network.dimensions,# + ["node_count"],
                    size=(21, 1),
                    default_value=self.y_axis,
                    key="y_axis",
                    # enable_events=True
                )
            ]
        )
        self.show_edges = False
        self.c_map = "RdYlGn"
        layout.append(
            [
                sg.Checkbox(
                    'Show edges',
                    default=self.show_edges,
                    key="show_edges"
                ),
                sg.Combo(
                    sorted(matplotlib.cm.__dict__['datad']),
                    size=(21, 1),
                    default_value=self.c_map,
                    key="c_map",
                    # enable_events=True
                )
            ]
        )
        self.bg_color = "white"
        layout.append(
            [
                sg.Combo(
                    sorted(matplotlib.colors.__dict__['CSS4_COLORS']),
                    size=(21, 1),
                    default_value=self.bg_color,
                    key="bg_color",
                    # enable_events=True
                )
            ]
        )
        self.min_positive_threshold = max_node_count
        self.max_positive_threshold = max_node_count
        layout.append(
            [
                sg.Text('POSITIVE EDGE EVIDENCE', size=(25, 1)),
                sg.InputText(
                    self.min_positive_threshold,
                    key=f"min_positive_edge_count",
                    size=(10, 1),
                    # enable_events=True
                ),
                sg.InputText(
                    self.max_positive_threshold,
                    key=f"max_positive_edge_count",
                    size=(10, 1),
                    # enable_events=True
                ),
            ]
        )
        self.min_negative_threshold = 0
        self.max_negative_threshold = 0
        layout.append(
            [
                sg.Text('NEGATIVE EDGE EVIDENCE', size=(25, 1)),
                sg.InputText(
                    self.min_negative_threshold,
                    key=f"min_negative_edge_count",
                    size=(10, 1),
                    # enable_events=True
                ),
                sg.InputText(
                    self.max_negative_threshold,
                    key=f"max_negative_edge_count",
                    size=(10, 1),
                    # enable_events=True
                ),
            ]
        )
        layout.append(
            [
                sg.Button('Refresh'),
                sg.Button("Save"),
                sg.InputText(
                    "",
                    key=f"file_name",
                    size=(10, 1),
                    # enable_events=True
                ),
                sg.InputText(
                    "1",
                    key="dpi",
                    size=(10, 1),
                    # enable_events=True
                )
            ]
        )
        self.overview_window = sg.Window('Settings', layout)

    def create_plot_window(self):
        # TODO: Docstring
        # self.fig = plt.figure(1, figsize=(13, 9))
        # a4 1200ppi 9921 x 14032 (40% is cropped of)
        self.fig = plt.figure(
            1,
            figsize=(
                (14032 / 0.6) / 2400,
                (9921 / 0.6) / 2400)
        )
        # self.fig = plt.figure(1, figsize=(29.7 / 2.5, 21 / 2.5))
        self.aggregate_ax = self.fig.add_subplot(111)
        self.aggregate_ax.set_facecolor('white')
        figure_x, figure_y, figure_w, figure_h = self.fig.bbox.bounds
        layout = [
            [sg.Canvas(size=(figure_w, figure_h), key='canvas')]
        ]
        self.plot_window = sg.Window('Plot', layout, finalize=True)
        mpl_backend = matplotlib.backends.backend_tkagg
        self.figure_canvas_agg = mpl_backend.FigureCanvasTkAgg(
            self.fig,
            self.plot_window['canvas'].TKCanvas
        )
        self.figure_canvas_agg.get_tk_widget().pack(
            side='top',
            fill='both',
            expand=1
        )
        self.update_plot()

    def update_plot(self):
        # TODO: Docstring
        nodes = np.flatnonzero(self.nodes == self.node_threshold)
        x_coordinates, y_coordinates = self.ion_network.get_ion_coordinates(
            [self.x_axis, self.y_axis],
            nodes
        )
        if not hasattr(self, "node_scatter"):
            self.node_scatter = self.aggregate_ax.scatter(
                x_coordinates,
                y_coordinates,
                marker=".",
                c="lightgrey"
            )
        else:
            self.node_scatter.set_offsets(np.c_[x_coordinates, y_coordinates])
        self.aggregate_ax.set_xlabel(self.x_axis)
        self.aggregate_ax.set_ylabel(self.y_axis)
        self.aggregate_ax.set_xlim(
            [np.min(x_coordinates), np.max(x_coordinates)]
        )
        self.aggregate_ax.set_ylim(
            [np.min(y_coordinates), np.max(y_coordinates)]
        )
        if not hasattr(self, "edge_collection"):
            self.edge_collection = self.aggregate_ax.add_collection(
                matplotlib.collections.LineCollection([], [])
            )
        else:
            if self.show_edges:
                selected_neighbors = self.edges[nodes].T.tocsr()[nodes]
                a, b = selected_neighbors.nonzero()
                positive_counts = self.positive_edge_evidence[
                    selected_neighbors.data
                ]
                negative_counts = self.negative_edge_evidence[
                    selected_neighbors.data
                ]
                selection = (positive_counts >= self.min_positive_threshold)
                selection &= (positive_counts <= self.max_positive_threshold)
                selection &= (negative_counts >= self.min_negative_threshold)
                selection &= (negative_counts <= self.max_negative_threshold)
                a = a[selection]
                b = b[selection]
                start_edges = list(zip(x_coordinates[a], y_coordinates[a]))
                end_edges = list(zip(x_coordinates[b], y_coordinates[b]))
                colors = positive_counts[selection] - negative_counts[selection]
                edges = np.array(list(zip(start_edges, end_edges)))
                color_order = np.argsort(colors)
                self.colors = colors[color_order]
                edges = edges[color_order]
            else:
                edges = []
                self.colors = []
            self.edge_collection.set_segments(edges)
            self.update_edge_colors()
        self.figure_canvas_agg.draw()
        self.figure_canvas_agg.flush_events()

    def update_edge_colors(self):
        edge_color_mapper = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(
                # vmin=-self.evidence.evidence_count,
                vmin=-self.evidence.evidence_count - 1,
                # vmin=0,
                vmax=self.evidence.evidence_count
            ),
            cmap=self.c_map
        )
        self.edge_collection.set_color(edge_color_mapper.to_rgba(self.colors))
        self.figure_canvas_agg.draw()
        self.figure_canvas_agg.flush_events()

    def update_window(self, window):
        # TODO: Docstring
        event, values = window.read(timeout=100)
        if event == sg.TIMEOUT_KEY:
            return
        if event in (None, 'Exit'):
            return "exit"
        if (event == "Save") and (values["file_name"] != ""):
            plt.savefig(values["file_name"], dpi=2 * int(values["dpi"]))
        if window.Title == "Settings":
            # print(event, values)
            self.perform_overview_action(event, values)
        if window.Title == "Plot":
            self.perform_plot_action(event, values)

    def perform_overview_action(self, event, values):
        # TODO: Docstring
        update = False
        if self.show_edges != values["show_edges"]:
            self.show_edges = values["show_edges"]
            update = True
        if self.x_axis != values["x_axis"]:
            self.x_axis = values["x_axis"]
            update = True
        if self.y_axis != values["y_axis"]:
            self.y_axis = values["y_axis"]
            update = True
        if self.min_positive_threshold != values["min_positive_edge_count"]:
            self.min_positive_threshold = float(
                values["min_positive_edge_count"]
            )
            update = True
        if self.max_positive_threshold != values["max_positive_edge_count"]:
            self.max_positive_threshold = float(
                values["max_positive_edge_count"]
            )
            update = True
        if self.min_negative_threshold != values["min_negative_edge_count"]:
            self.min_negative_threshold = float(
                values["min_negative_edge_count"]
            )
            update = True
        if self.max_negative_threshold != values["max_negative_edge_count"]:
            self.max_negative_threshold = float(
                values["max_negative_edge_count"]
            )
            update = True
        for key, (low, high) in list(self.node_dimensions.items()):
            if key == "node_evidence":
                if low != float(values["min_node_count"]):
                    update = True
                    node_count = self.evidence.get_aligned_nodes_from_group()
                    self.nodes -= node_count >= low
                    self.nodes += node_count >= float(values["min_node_count"])
                    self.node_dimensions["node_evidence"][0] = float(
                        values["min_node_count"]
                    )
                if high != float(values["max_node_count"]):
                    update = True
                    node_count = self.evidence.get_aligned_nodes_from_group()
                    self.nodes -= node_count <= high
                    self.nodes += node_count <= float(values["max_node_count"])
                    self.node_dimensions["node_evidence"][1] = float(
                        values["max_node_count"]
                    )
            else:
                if low != float(values[f"min_{key}"]):
                    update = True
                    coordinates = self.ion_network.get_ion_coordinates(key)
                    self.nodes -= coordinates >= low
                    self.nodes += coordinates >= float(values[f"min_{key}"])
                    self.node_dimensions[key][0] = float(values[f"min_{key}"])
                if high != float(values[f"max_{key}"]):
                    update = True
                    coordinates = self.ion_network.get_ion_coordinates(key)
                    self.nodes -= coordinates <= high
                    self.nodes += coordinates <= float(values[f"max_{key}"])
                    self.node_dimensions[key][1] = float(values[f"max_{key}"])
        if update:
            # print(np.bincount(self.nodes == self.node_threshold))
            self.update_plot()
        if self.c_map != values["c_map"]:
            self.c_map = values["c_map"]
            self.update_edge_colors()
        if self.bg_color != values["bg_color"]:
            self.bg_color = values["bg_color"]
            self.aggregate_ax.set_facecolor(self.bg_color)
            self.figure_canvas_agg.draw()
            self.figure_canvas_agg.flush_events()

    def perform_plot_action(self, event, values):
        # TODO: Docstring
        print("plot", event, values)

#
# def updateSelectedNodes(self, gui, event):
#     self.log.printMessage("Updating selected nodes")
#     rt = event.xdata
#     dt = event.ydata
#     rt_low, rt_high, dt_low, dt_high = gui.getVisibleBoundaries()
#     rt_distance = (
#         self.anchors["RT"][self.visible_nodes] - rt
#     ) / (rt_high - rt_low)
#     dt_distance = (
#         self.anchors["DT"][self.visible_nodes] - dt
#     ) / (dt_high - dt_low)
#     distance = np.sqrt(rt_distance**2 + dt_distance**2)
#     if event.key != "control":
#         self.unselectAllVisible()
#     selected_index = np.argmin(distance)
#     self.selected_nodes[selected_index] = not self.selected_nodes[selected_index]
#     self.updateLabelSelection(gui)
#
# def onclick(event):
#     if (event is None) or (not event.dblclick):
#         return
#     self.dataset.updateSelectedNodes(self, event)
#     self.dataset.plotSelectedNodes(self)
#     self.dataset.plotIons(self)
#     self.refreshAggregateCanvas()
#     self.refreshIonCanvas()
#     self.aggregate_canvas.mpl_disconnect(self.click_connection)
#     self.click_connection = self.aggregate_canvas.mpl_connect(
#         'button_press_event',
#         onclick
#     )
#
# self.click_connection = self.aggregate_canvas.mpl_connect(
#     'button_press_event',
#     onclick
# )
