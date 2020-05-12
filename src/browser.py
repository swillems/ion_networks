#!python

# builtin
import contextlib
# external
import numpy as np
import matplotlib.pyplot as plt
import PySimpleGUI as sg
import matplotlib
# local
import ms_run_files


matplotlib.use('TkAgg')
# plt.rcParams['toolbar'] = 'toolmanager'

MATPLOTLIB_COLORS_CONTINUOUS = sorted(matplotlib.cm.__dict__['datad'])
MATPLOTLIB_COLORS_FIXED = sorted(matplotlib.colors.__dict__['CSS4_COLORS'])


class Browser(object):
    # TODO: Docstring

    def __init__(
        self,
        start=False,
        widget_size=20,
    ):
        # TODO: Docstring
        self.widget_size = widget_size
        self.window = {}
        self.evaluate_window = {}
        self.init_main_window()
        self.window["Main"] = sg.Window("Main", self.window["Main"])
        self.active_window_name = "Main"
        self.figs = {
            "network": self.init_figure("network"),
            "evidence": self.init_figure("evidence"),
        }
        plt.show(block=False)
        if start:
            self.run()
            self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, type=None, value=None, traceback=None):
        for window in list(self.window.values()):
            if not isinstance(window, list):
                window.close()
        for fig in list(self.figs.values()):
            plt.close(fig)

    def run(self):
        # TODO: Docstring
        while self.active_window_name is not None:
            window = self.window[self.active_window_name]
            event, values = window.read(timeout=10)
            if event == sg.TIMEOUT_KEY:
                continue
            elif event is None:
                self.active_window_name = None
            elif event == "Return to main menu":
                self.swap_active_window("Main")
            else:
                self.evaluate_window[self.active_window_name](event, values)

    def init_main_window(self):
        self.window["Main"] = [
            [
                sg.Text("Ion-network", size=(self.widget_size, 1)),
                sg.Button(
                    'None',
                    size=(self.widget_size, 1),
                    key='Select ion-network'
                ),
            ],
            [
                sg.Text("Nodes", size=(self.widget_size, 1)),
                sg.Button(
                    'Settings',
                    size=(self.widget_size, 1),
                    key="Node settings"
                )
            ],
            [
                sg.Text("Edges", size=(self.widget_size, 1)),
                sg.Button(
                    'Settings',
                    size=(self.widget_size, 1),
                    key="Edge settings"
                )
            ],
            [
                sg.Text("Compare", size=(self.widget_size, 1)),
                sg.Button(
                    'Compare',
                    size=(self.widget_size, 1),
                    key="Compare settings"
                )
            ],
            [
                sg.Text("Misc", size=(self.widget_size, 1)),
                sg.Button(
                    'Misc',
                    size=(self.widget_size, 1),
                    key="Misc settings"
                )
            ],
            # Compare with other samples
            # Other settings (bg color, save, export)
        ]
        self.evaluate_window["Main"] = self.evaluate_main_window

    def init_node_settings_window(self):
        try:
            del self.nodes
        except AttributeError:
            pass
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
                ),
                sg.InputText(
                    max_coordinate,
                    key=f"max_{dimension}",
                    size=(10, 1),
                ),
            ]
            layout.append(row)
        max_node_count = float(self.evidence.evidence_count)
        self.node_dimensions["node_evidence"] = [max_node_count, max_node_count]
        self.node_threshold = 2 * len(self.node_dimensions)
        self.nodes = np.repeat(self.node_threshold, self.ion_network.node_count)
        node_evidence = self.evidence.get_aligned_nodes_from_group()
        self.nodes -= node_evidence < max_node_count
        layout.append(
            [
                sg.Text('NODE EVIDENCE', size=(25, 1)),
                sg.InputText(
                    max_node_count,
                    key=f"min_node_count",
                    size=(10, 1),
                ),
                sg.InputText(
                    max_node_count,
                    key=f"max_node_count",
                    size=(10, 1),
                ),
            ]
        )
        self.x_axis = "PRECURSOR_RT"
        layout.append(
            [
                sg.Text('X-AXIS', size=(25, 1)),
                sg.Combo(
                    self.ion_network.dimensions,
                    size=(21, 1),
                    default_value=self.x_axis,
                    key="x_axis",
                )
            ]
        )
        self.y_axis = "FRAGMENT_MZ"
        layout.append(
            [
                sg.Text('Y-AXIS', size=(25, 1)),
                sg.Combo(
                    self.ion_network.dimensions,
                    size=(21, 1),
                    default_value=self.y_axis,
                    key="y_axis",
                )
            ]
        )
        self.node_color = "FRAGMENT_LOGINT"
        # self.node_color = "lightgrey"
        self.node_color_c_map = "RdYlGn"
        self.node_color_normalize = True
        layout.append(
            [
                sg.Text('NODE COLOR', size=(25, 1)),
                sg.Combo(
                    self.ion_network.dimensions + ['NODE EVIDENCE'] + MATPLOTLIB_COLORS_FIXED,
                    size=(21, 1),
                    default_value=self.node_color,
                    key="node_color",
                ),
                sg.Combo(
                    MATPLOTLIB_COLORS_CONTINUOUS,
                    size=(21, 1),
                    default_value=self.node_color_c_map,
                    key="node_color_c_map",
                ),
                sg.Checkbox(
                    'Relative',
                    default=self.node_color_normalize,
                    key="node_color_normalize"
                ),
            ]
        )
        layout.append(self.add_main_menu_and_continue_buttons_to_layout())
        self.window["Node settings"] = sg.Window('Node settings', layout)
        self.evaluate_window["Node settings"] = self.evaluate_node_settings_window

    def init_edge_settings_window(self):
        # TODO: Docstring
        try:
            del self.edges
            del self.positive_edge_evidence
            del self.negative_edge_evidence
        except AttributeError:
            pass
        layout = []
        self.edges = self.ion_network.get_edges(data_as_index=True)
        self.positive_edge_evidence = self.evidence.get_edge_mask_from_group()
        self.negative_edge_evidence = self.evidence.get_edge_mask_from_group(
            positive=False
        )
        self.show_edges = False
        layout.append(
            [
                sg.Checkbox(
                    'SHOW EDGES',
                    default=self.show_edges,
                    key="show_edges"
                ),
            ]
        )
        max_node_count = float(self.evidence.evidence_count)
        self.min_positive_threshold = max_node_count
        self.max_positive_threshold = max_node_count
        layout.append(
            [
                sg.Text('POSITIVE EDGE EVIDENCE', size=(25, 1)),
                sg.InputText(
                    self.min_positive_threshold,
                    key=f"min_positive_edge_count",
                    size=(10, 1),
                ),
                sg.InputText(
                    self.max_positive_threshold,
                    key=f"max_positive_edge_count",
                    size=(10, 1),
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
                ),
                sg.InputText(
                    self.max_negative_threshold,
                    key=f"max_negative_edge_count",
                    size=(10, 1),
                ),
            ]
        )
        self.min_summed_threshold = max_node_count
        self.max_summed_threshold = max_node_count
        layout.append(
            [
                sg.Text('SUMMED EDGE EVIDENCE', size=(25, 1)),
                sg.InputText(
                    self.min_summed_threshold,
                    key=f"min_summed_edge_count",
                    size=(10, 1),
                ),
                sg.InputText(
                    self.max_summed_threshold,
                    key=f"max_summed_edge_count",
                    size=(10, 1),
                ),
            ]
        )
        self.edge_color = "SUMMED"
        # self.edge_color = "lightgrey"
        self.edge_color_c_map = "RdYlGn"
        self.edge_color_normalize = True
        layout.append(
            [
                sg.Text('EDGE COLOR', size=(25, 1)),
                sg.Combo(
                    ['POSITIVE', 'NEGATIVE', 'SUMMED'] + MATPLOTLIB_COLORS_FIXED,
                    size=(21, 1),
                    default_value=self.edge_color,
                    key="edge_color",
                ),
                sg.Combo(
                    MATPLOTLIB_COLORS_CONTINUOUS,
                    size=(21, 1),
                    default_value=self.edge_color_c_map,
                    key="edge_color_c_map",
                ),
                sg.Checkbox(
                    'Relative',
                    default=self.edge_color_normalize,
                    key="edge_color_normalize"
                ),
            ]
        )
        layout.append(self.add_main_menu_and_continue_buttons_to_layout())
        self.window["Edge settings"] = sg.Window('Edge settings', layout)
        self.evaluate_window["Edge settings"] = self.evaluate_edge_settings_window

    def init_compare_settings_window(self):
        layout = []
        self.evidence_axis = "FRAGMENT_LOGINT"
        layout.append(
            [
                sg.Text('X-AXIS', size=(25, 1)),
                sg.Combo(
                    self.ion_network.dimensions,
                    size=(21, 1),
                    default_value=self.evidence_axis,
                    key="evidence_axis",
                )
            ]
        )
        layout.append(self.add_main_menu_and_continue_buttons_to_layout())
        self.window["Compare settings"] = sg.Window('Compare settings', layout)
        self.evaluate_window["Compare settings"] = self.evaluate_compare_settings_window

    def init_misc_settings_window(self):
        layout = []
        self.network_bg_color = "white"
        layout.append(
            [
                sg.Text('NETWORK BACKGROUND COLOR', size=(25, 1)),
                sg.Combo(
                    MATPLOTLIB_COLORS_FIXED,
                    size=(21, 1),
                    default_value=self.network_bg_color,
                    key="network_bg_color",
                    # enable_events=True
                )
            ]
        )
        self.evidence_bg_color = "white"
        layout.append(
            [
                sg.Text('EVIDENCE BACKGROUND COLOR', size=(25, 1)),
                sg.Combo(
                    MATPLOTLIB_COLORS_FIXED,
                    size=(21, 1),
                    default_value=self.evidence_bg_color,
                    key="evidence_bg_color",
                    # enable_events=True
                )
            ]
        )
        layout.append(
            [
                sg.Button("Save network"),
                sg.Button("Save evidence")
            ]
        )
        self.figs["network"].axes[0].set_facecolor(self.network_bg_color)
        self.figs["evidence"].axes[0].set_facecolor(self.evidence_bg_color)
        layout.append(self.add_main_menu_and_continue_buttons_to_layout())
        self.window["Misc settings"] = sg.Window('Misc settings', layout)
        self.evaluate_window["Misc settings"] = self.evaluate_misc_settings_window

    def evaluate_main_window(self, event, values):
        # TODO: Docstring
        if event == "Select ion-network":
            self.init_network()
        else:
            if event not in self.window:
                sg.popup_error('Please select an ion-network first.')
                return
        if event == 'Node settings':
            self.swap_active_window("Node settings")
        if event == 'Edge settings':
            self.swap_active_window("Edge settings")
        if event == 'Compare settings':
            self.swap_active_window("Compare settings")
        if event == 'Misc settings':
            self.swap_active_window("Misc settings")

    def evaluate_misc_settings_window(
        self,
        event,
        values,
        update_node_selection=False,
        update_node_colors=False,
    ):
        if self.network_bg_color != values["network_bg_color"]:
            self.network_bg_color = values["network_bg_color"]
            self.figs["network"].axes[0].set_facecolor(self.network_bg_color)
            self.figs["network"].canvas.draw()
            self.figs["network"].canvas.flush_events()
        if self.evidence_bg_color != values["evidence_bg_color"]:
            self.evidence_bg_color = values["evidence_bg_color"]
            self.figs["evidence"].axes[0].set_facecolor(self.evidence_bg_color)
            self.figs["evidence"].canvas.draw()
            self.figs["evidence"].canvas.flush_events()
        if event == "Save network":
            file_name = sg.popup_get_file("Save network", save_as=True)
            if file_name is None:
                return
            with loading_window():
                self.figs["network"].savefig(file_name)
        if event == "Save evidence":
            file_name = sg.popup_get_file("Save evidence", save_as=True)
            if file_name is None:
                return
            with loading_window():
                self.figs["evidence"].savefig(file_name)

    def evaluate_compare_settings_window(
        self,
        event,
        values,
        update_node_selection=False,
        update_node_colors=False,
    ):
        if self.evidence_axis != values["evidence_axis"]:
            self.evidence_axis = values["evidence_axis"]
            self.figure_recenter(
                self.figs["evidence"],
                reset_axis="x",
                flush=False
            )
            self.evidence_figure_update_axis_selection()

    def evaluate_node_settings_window(
        self,
        event,
        values,
        update_node_selection=False,
        update_node_colors=False,
    ):
        for key, (low, high) in list(self.node_dimensions.items()):
            if key == "node_evidence":
                if low != float(values["min_node_count"]):
                    update_node_selection = True
                    node_count = self.evidence.get_aligned_nodes_from_group()
                    self.nodes -= node_count >= low
                    self.nodes += node_count >= float(values["min_node_count"])
                    self.node_dimensions["node_evidence"][0] = float(
                        values["min_node_count"]
                    )
                if high != float(values["max_node_count"]):
                    update_node_selection = True
                    node_count = self.evidence.get_aligned_nodes_from_group()
                    self.nodes -= node_count <= high
                    self.nodes += node_count <= float(values["max_node_count"])
                    self.node_dimensions["node_evidence"][1] = float(
                        values["max_node_count"]
                    )
            else:
                if low != float(values[f"min_{key}"]):
                    update_node_selection = True
                    coordinates = self.ion_network.get_ion_coordinates(key)
                    self.nodes -= coordinates >= low
                    self.nodes += coordinates >= float(values[f"min_{key}"])
                    self.node_dimensions[key][0] = float(values[f"min_{key}"])
                if high != float(values[f"max_{key}"]):
                    update_node_selection = True
                    coordinates = self.ion_network.get_ion_coordinates(key)
                    self.nodes -= coordinates <= high
                    self.nodes += coordinates <= float(values[f"max_{key}"])
                    self.node_dimensions[key][1] = float(values[f"max_{key}"])
        if self.x_axis != values["x_axis"]:
            self.x_axis = values["x_axis"]
            self.figure_recenter(
                self.figs["network"],
                reset_axis="x",
                flush=False
            )
            update_node_selection = True
        if self.y_axis != values["y_axis"]:
            self.y_axis = values["y_axis"]
            self.figure_recenter(
                self.figs["network"],
                reset_axis="y",
                flush=False
            )
            update_node_selection = True
        if update_node_selection:
            (
                nodes,
                x_coordinates,
                y_coordinates
            ) = self.network_figure_update_node_selection(flush=False)
            selected_edges = self.network_figure_update_edge_selection(
                nodes,
                x_coordinates,
                y_coordinates,
                flush=False
            )
            self.network_figure_update_edge_colors(selected_edges)
            update_node_colors = True
        if self.node_color != values["node_color"]:
            self.node_color = values["node_color"]
            update_node_colors = True
        if self.node_color_c_map != values["node_color_c_map"]:
            self.node_color_c_map = values["node_color_c_map"]
            update_node_colors = True
        if self.node_color_normalize != values["node_color_normalize"]:
            self.node_color_normalize = values["node_color_normalize"]
            update_node_colors |= (
                self.node_color in self.ion_network.dimensions + ['NODE EVIDENCE']
            )
        if update_node_colors:
            try:
                self.network_figure_update_node_colors(nodes)
            except NameError:
                self.network_figure_update_node_colors()

    def evaluate_edge_settings_window(
        self,
        event,
        values,
        update_edge_selection=False,
        update_edge_colors=False,
    ):
        # TODO: Docstring
        if self.show_edges != values["show_edges"]:
            self.show_edges = values["show_edges"]
            update_edge_selection = True
        if self.min_positive_threshold != values["min_positive_edge_count"]:
            self.min_positive_threshold = float(
                values["min_positive_edge_count"]
            )
            update_edge_selection = True
        if self.max_positive_threshold != values["max_positive_edge_count"]:
            self.max_positive_threshold = float(
                values["max_positive_edge_count"]
            )
            update_edge_selection = True
        if self.min_negative_threshold != values["min_negative_edge_count"]:
            self.min_negative_threshold = float(
                values["min_negative_edge_count"]
            )
            update_edge_selection = True
        if self.max_negative_threshold != values["max_negative_edge_count"]:
            self.max_negative_threshold = float(
                values["max_negative_edge_count"]
            )
            update_edge_selection = True
        if self.min_summed_threshold != values["min_summed_edge_count"]:
            self.min_summed_threshold = float(
                values["min_summed_edge_count"]
            )
            update_edge_selection = True
        if self.max_summed_threshold != values["max_summed_edge_count"]:
            self.max_summed_threshold = float(
                values["max_summed_edge_count"]
            )
            update_edge_selection = True
        if update_edge_selection:
            selected_edges = self.network_figure_update_edge_selection(flush=False)
            update_edge_colors = True
        if self.edge_color != values["edge_color"]:
            self.edge_color = values["edge_color"]
            update_edge_colors = True
        if self.edge_color_c_map != values["edge_color_c_map"]:
            self.edge_color_c_map = values["edge_color_c_map"]
            update_edge_colors = True
        if self.edge_color_normalize != values["edge_color_normalize"]:
            self.edge_color_normalize = values["edge_color_normalize"]
            update_edge_colors |= (
                self.edge_color in ['POSITIVE', 'NEGATIVE', 'SUMMED']
            )
        if update_edge_colors:
            self.network_figure_update_edge_colors(selected_edges)

    def network_figure_update_node_colors(
        self,
        nodes=None,
        reorder=True,
        flush=True,
    ):
        # TODO: Docstring
        if nodes is None:
            nodes = np.flatnonzero(self.nodes == self.node_threshold)
        if self.node_color in self.ion_network.dimensions + ['NODE EVIDENCE']:
            if self.node_color in self.ion_network.dimensions:
                inds = self.ion_network.get_ion_coordinates(self.node_color)
            elif self.node_color == 'NODE EVIDENCE':
                inds = self.evidence.get_aligned_nodes_from_group()
            color_inds = inds[nodes]
            if reorder:
                offsets = np.array(self.node_scatter.get_offsets())
                x_coordinates = offsets[:, 0]
                y_coordinates = offsets[:, 1]
                order = np.argsort(color_inds)
                color_inds = color_inds[order]
                x_coordinates = x_coordinates[order]
                y_coordinates = y_coordinates[order]
                self.node_scatter.set_offsets(
                    np.c_[x_coordinates, y_coordinates]
                )
            if not self.node_color_normalize:
                vmin = np.min(inds)
                vmax = np.max(inds)
            else:
                vmin = np.min(color_inds)
                vmax = np.max(color_inds)
            color_mapper = matplotlib.cm.ScalarMappable(
                norm=matplotlib.colors.Normalize(
                    vmin=vmin,
                    vmax=vmax,
                ),
                cmap=self.node_color_c_map
            )
            colors = color_mapper.to_rgba(color_inds)
        else:
            colors = [matplotlib.colors.to_rgba(self.node_color)] * len(nodes)
        self.node_scatter.set_facecolor(colors)
        if flush:
            self.figs["network"].canvas.draw()
            self.figs["network"].canvas.flush_events()

    def network_figure_update_edge_colors(
        self,
        selected_edges=None,
        reorder=True,
        flush=True,
    ):
        # TODO: Docstring
        if self.show_edges:
            if selected_edges is None:
                nodes = np.flatnonzero(self.nodes == self.node_threshold)
                selected_neighbors = self.edges[nodes].T.tocsr()[nodes]
                positive_counts = self.positive_edge_evidence[
                    selected_neighbors.data
                ]
                negative_counts = self.negative_edge_evidence[
                    selected_neighbors.data
                ]
                summed_counts = positive_counts - negative_counts
                selection = (positive_counts >= self.min_positive_threshold)
                selection &= (positive_counts <= self.max_positive_threshold)
                selection &= (negative_counts >= self.min_negative_threshold)
                selection &= (negative_counts <= self.max_negative_threshold)
                selection &= (summed_counts >= self.min_summed_threshold)
                selection &= (summed_counts <= self.max_summed_threshold)
                selected_edges = selected_neighbors.data[selection]
            if self.edge_color in ['POSITIVE', 'NEGATIVE', 'SUMMED']:
                if self.edge_color == 'POSITIVE':
                    inds = self.positive_edge_evidence
                if self.edge_color == 'NEGATIVE':
                    inds = self.negative_edge_evidence
                elif self.edge_color == 'SUMMED':
                    inds = self.positive_edge_evidence - self.negative_edge_evidence
                color_inds = inds[selected_edges]
                if reorder:
                    offsets = np.array(self.edge_collection.get_segments())
                    order = np.argsort(color_inds)
                    color_inds = color_inds[order]
                    offsets = offsets[order]
                    self.edge_collection.set_segments(offsets)
                if not self.edge_color_normalize:
                    vmin = np.min(inds)
                    vmax = np.max(inds)
                else:
                    vmin = np.min(color_inds)
                    vmax = np.max(color_inds)
                color_mapper = matplotlib.cm.ScalarMappable(
                    norm=matplotlib.colors.Normalize(
                        vmin=vmin,
                        vmax=vmax,
                    ),
                    cmap=self.edge_color_c_map
                )
                colors = color_mapper.to_rgba(color_inds)
            else:
                colors = [
                    matplotlib.colors.to_rgba(self.edge_color)
                ] * len(selected_edges)
        else:
            colors = []
        self.edge_collection.set_color(colors)
        if flush:
            self.figs["network"].canvas.draw()
            self.figs["network"].canvas.flush_events()

    def network_figure_update_edge_selection(
        self,
        nodes=None,
        x_coordinates=None,
        y_coordinates=None,
        flush=True
    ):
        # TODO: Docstring
        if not hasattr(self, "edge_collection"):
            self.edge_collection = self.figs["network"].axes[0].add_collection(
                matplotlib.collections.LineCollection([], [])
            )
        if self.show_edges:
            if nodes is None:
                nodes = np.flatnonzero(self.nodes == self.node_threshold)
            if (x_coordinates is None) or (y_coordinates is None):
                x_coordinates, y_coordinates = self.ion_network.get_ion_coordinates(
                    [self.x_axis, self.y_axis],
                    nodes
                )
            selected_neighbors = self.edges[nodes].T.tocsr()[nodes]
            a, b = selected_neighbors.nonzero()
            positive_counts = self.positive_edge_evidence[
                selected_neighbors.data
            ]
            negative_counts = self.negative_edge_evidence[
                selected_neighbors.data
            ]
            summed_counts = positive_counts - negative_counts
            selection = (positive_counts >= self.min_positive_threshold)
            selection &= (positive_counts <= self.max_positive_threshold)
            selection &= (negative_counts >= self.min_negative_threshold)
            selection &= (negative_counts <= self.max_negative_threshold)
            selection &= (summed_counts >= self.min_summed_threshold)
            selection &= (summed_counts <= self.max_summed_threshold)
            a = a[selection]
            b = b[selection]
            start_edges = list(zip(x_coordinates[a], y_coordinates[a]))
            end_edges = list(zip(x_coordinates[b], y_coordinates[b]))
            edges = np.array(list(zip(start_edges, end_edges)))
            selected_edges = selected_neighbors.data[selection]
        else:
            edges = []
            selected_edges = []
        self.edge_collection.set_segments(edges)
        if flush:
            self.figs["network"].canvas.draw()
            self.figs["network"].canvas.flush_events()
        return selected_edges

    def network_figure_update_node_selection(self, flush=True):
        # TODO: Docstring
        nodes = np.flatnonzero(self.nodes == self.node_threshold)
        x_coordinates, y_coordinates = self.ion_network.get_ion_coordinates(
            [self.x_axis, self.y_axis],
            nodes
        )
        if not hasattr(self, "node_scatter"):
            self.node_scatter = self.figs["network"].axes[0].scatter(
                x_coordinates,
                y_coordinates,
                marker=".",
            )
        else:
            self.node_scatter.set_offsets(np.c_[x_coordinates, y_coordinates])
        if flush:
            self.figs["network"].canvas.draw()
            self.figs["network"].canvas.flush_events()
        return nodes, x_coordinates, y_coordinates

    def evidence_figure_update_axis_selection(self, flush=True):
        # TODO: Docstring
        node_mask = np.zeros(self.ion_network.node_count, np.int64)
        nodes = np.flatnonzero(self.nodes == self.node_threshold)
        node_mask[nodes] = np.arange(len(nodes))
        alignments = np.zeros((len(nodes), self.evidence.evidence_count + 1))
        alignments[:, 0] = self.ion_network.get_ion_coordinates(
            self.evidence_axis,
            nodes
        )
        for i, other_run in enumerate(self.evidence.network_keys):
            evidence = self.evidence.get_other_run(other_run)
            alignment = self.evidence.get_alignment(
                other=evidence
            )
            self_ions = node_mask[alignment[:, 0]]
            other_ions = alignment[:, 1]
            alignments[self_ions, i + 1] = evidence.ion_network.get_ion_coordinates(
                self.evidence_axis,
                other_ions
            )
        if not hasattr(self, "evidence_plot"):
            self.evidence_plot = self.figs["evidence"].axes[0].plot(
                alignments.T
            )
            self.figs["evidence"].axes[0].set_xticks(
                list(range(self.evidence.evidence_count + 1))
            )
            self.figs["evidence"].axes[0].set_xticklabels(
                [self.ion_network.run_name] + self.evidence.network_keys,
                rotation=45,
                ha="right"
            )
        else:
            for i in self.evidence_plot:
                i.remove()
            del self.evidence_plot[:]
            self.evidence_plot = self.figs["evidence"].axes[0].plot(
                alignments.T
            )
        if flush:
            self.figs["evidence"].canvas.draw()
            self.figs["evidence"].canvas.flush_events()

    def init_figure(self, title):
        fig = plt.figure()
        fig.add_subplot(111)
        fig.canvas.toolbar.pack_forget()
        # for tool in list(fig.canvas.manager.toolmanager.tools):
        #     fig.canvas.manager.toolmanager.remove_tool(tool)
        fig.canvas.manager.window.protocol(
            "WM_DELETE_WINDOW",
            lambda: self.swap_active_window(None)
        )
        fig.canvas.set_window_title(title)
        fig.canvas.mpl_connect(
            'scroll_event',
            lambda event: self.figure_zoom(
                fig,
                x=event.xdata,
                y=event.ydata,
                zoom=event.button,
            )
        )
        fig.canvas.mpl_connect(
            'button_press_event',
            lambda event: self.mouse_event(fig, event)
        )
        return fig

    def _reset_axis(self, fig, axis):
        try:
            ax = fig.axes[0]
            # TODO what if it is evidence fig?
            if "x" in axis.lower():
                ax.set_xlabel(self.x_axis)
                ax.set_xlim(
                    [
                        self.node_dimensions[self.x_axis][0],
                        self.node_dimensions[self.x_axis][1],
                    ]
                )
            if "y" in axis.lower():
                ax.set_ylabel(self.y_axis)
                ax.set_ylim(
                    [
                        self.node_dimensions[self.y_axis][0],
                        self.node_dimensions[self.y_axis][1],
                    ]
                )
        except AttributeError:
            pass

    def figure_zoom(self, fig, x=None, y=None, zoom=None, flush=True):
        ax = fig.axes[0]
        if x is not None:
            x_min, x_max = ax.get_xlim()
            x_center = (x_min + x_max) / 2
            x_range = x_max - x_min
            x_offset = 0.5 - np.abs(x - x_center) / x_range
            if zoom == "down":
                x_scale = 1 + x_offset
            elif zoom == "up":
                x_scale = 1 - x_offset
            new_x_lim = [
                x_center - x_range / 2 * x_scale,
                x_center + x_range / 2 * x_scale
            ]
            ax.set_xlim(new_x_lim)
        if x is not None:
            y_min, y_max = ax.get_ylim()
            y_center = (y_min + y_max) / 2
            y_range = y_max - y_min
            y_offset = 0.5 - np.abs(y - y_center) / y_range
            if zoom == "down":
                y_scale = 1 + y_offset
            elif zoom == "up":
                y_scale = 1 - y_offset
            new_y_lim = [
                y_center - y_range / 2 * y_scale,
                y_center + y_range / 2 * y_scale
            ]
            ax.set_ylim(new_y_lim)
        if flush:
            fig.canvas.draw()
            fig.canvas.flush_events()

    def mouse_event(self, fig, event):
        if event.button == matplotlib.backend_bases.MouseButton.MIDDLE:
            reset_axis = "xy" if event.dblclick else None
            self.figure_recenter(
                fig,
                x=event.xdata,
                y=event.ydata,
                reset_axis=reset_axis
            )

    def figure_recenter(self, fig, x=None, y=None, reset_axis=None, flush=True):
        if reset_axis is not None:
            self._reset_axis(fig, reset_axis)
        else:
            ax = fig.axes[0]
            if x is not None:
                x_min, x_max = ax.get_xlim()
                x_range = x_max - x_min
                new_x_lim = [
                    x - x_range / 2,
                    x + x_range / 2
                ]
                ax.set_xlim(new_x_lim)
            if y is not None:
                y_min, y_max = ax.get_ylim()
                y_range = y_max - y_min
                new_y_lim = [
                    y - y_range / 2,
                    y + y_range / 2
                ]
                ax.set_ylim(new_y_lim)
        if flush:
            fig.canvas.draw()
            fig.canvas.flush_events()

    def add_main_menu_and_continue_buttons_to_layout(
        self,
        main_menu_button=True,
        continue_button=True
    ):
        # TODO: Docstring
        row = []
        if main_menu_button:
            row.append(
                sg.Button("Return to main menu", size=(self.widget_size, 1))
            )
        if continue_button:
            row.append(sg.Button("Submit", size=(self.widget_size, 1)))
        return row

    def swap_active_window(self, new_window_name=""):
        # TODO: Docstring, implement
        if new_window_name is None:
            self.active_window_name = None
        if self.active_window_name is None:
            return
        if self.active_window_name != "":
            self.window[self.active_window_name].Hide()
        if new_window_name != "":
            if isinstance(self.window[new_window_name], list):
                self.window[new_window_name] = sg.Window(
                    new_window_name,
                    self.window[new_window_name]
                )
            self.window[new_window_name].UnHide()
        self.active_window_name = new_window_name

    def init_network(self):
        file_name = sg.popup_get_file(
            'Please select an ion-network',
            file_types=(('Ion-network', '*.inet.hdf'),)
        )
        if file_name is None:
            return
        try:
            ion_network = ms_run_files.Network(file_name)
        except (OSError, ValueError):
            sg.popup_error('This is not a valid ion_network')
            return
        try:
            evidence = ms_run_files.Evidence(file_name)
        except (OSError, ValueError):
            sg.popup_error('This ion_network has no valid evidence')
            return
        if hasattr(self, "ion_network") and (self.ion_network == ion_network):
            pass
        else:
            self.window[self.active_window_name].Hide()
            opts = 5
            for i in range(opts):
                sg.one_line_progress_meter(
                    'Loading ion-network and evidence',
                    i,
                    opts - 1,
                    'load_window',
                    'Please wait...',
                    orientation="h",
                )
                if i == 0:
                    self.ion_network = ion_network
                    self.evidence = evidence
                    self.init_node_settings_window()
                    self.init_misc_settings_window()
                if i == 1:
                    self.init_edge_settings_window()
                if i == 2:
                    self.init_compare_settings_window()
                if i == 3:
                    self.figure_recenter(
                        self.figs["network"],
                        reset_axis="xy"
                    )
                    (
                        nodes,
                        x_coordinates,
                        y_coordinates
                    ) = self.network_figure_update_node_selection(flush=False)
                    self.network_figure_update_node_colors(nodes, flush=False)
                    selected_edges = self.network_figure_update_edge_selection(
                        nodes=nodes,
                        x_coordinates=x_coordinates,
                        y_coordinates=y_coordinates,
                        flush=False
                    )
                    self.network_figure_update_edge_colors(selected_edges)
                    self.window["Main"]["Select ion-network"](
                        ion_network.run_name
                    )
            self.window[self.active_window_name].UnHide()


@contextlib.contextmanager
def loading_window():
    # TODO: Docstring
    sg.one_line_progress_meter(
        'Please wait',
        0,
        1,
        'load_window',
        'Please wait...',
        orientation="h",
    )
    yield None
    sg.one_line_progress_meter(
        'Please wait',
        1,
        1,
        'load_window',
        'Please wait...',
        orientation="h",
    )
