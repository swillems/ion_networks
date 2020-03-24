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
import utils


BROWSER_IMAGES_PATH = os.path.join(
    utils.LIB_PATH,
    "browser_images"
)
LOADING_ANIMATION_FILE_NAME = os.path.join(
    BROWSER_IMAGES_PATH,
    "load_animation.gif"
)


class Browser(object):
    # TODO: Docstring

    def __init__(self, ion_network, evidence, logger):
        # TODO: Docstring
        self.ion_network = ion_network
        self.evidence = evidence
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
        max_node_count = float(self.evidence.network_count)
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
                vmin=-self.evidence.network_count,
                # vmin=0,
                vmax=self.evidence.network_count
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

    def perform_plot_action(self, event, values):
        # TODO: Docstring
        print("plot", event, values)
