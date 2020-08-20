#!python

# builtin
import contextlib
import os
import warnings
# external
import numpy as np
import numexpr as ne
import matplotlib.pyplot as plt
import PySimpleGUI as sg
import matplotlib
import pandas as pd
# local
try:
    from . import ms_run_files
    from . import ms_utils
except (ImportError, ModuleNotFoundError):
    import ms_run_files
    me_utils


matplotlib.use('TkAgg')
plt.rcParams['toolbar'] = 'toolmanager'

BASE_PATH = ms_utils.BASE_PATH
LIB_PATH = ms_utils.LIB_PATH
DEFAULT_BROWSER_PATH = os.path.join(LIB_PATH, "browser_images")
DEFAULT_BROWSER_IMAGES = {
    "pointer": "pointer_25x25.png",
}
NODE_LABEL_DECIMALS = 2

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
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.widget_size = widget_size
            self.window = {}
            self.evaluate_window = {}
            self.init_main_window()
            self.window["Main"] = sg.Window(
                "Main", self.window["Main"], finalize=True
            )
            self.active_window_name = "Main"
            self.figs = {
                "network": self.init_figure("network"),
                "evidence": self.init_figure("evidence"),
            }
            self.figs["network"].canvas.manager.toolmanager.add_tool(
                'node_select',
                PointerTool,
                callback_function=self.evidence_figure_update_axis_selection
            )
            self.figs["network"].canvas.manager.toolbar.add_tool(
                'node_select',
                'navigation',
                1
            )
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
        ]
        self.evaluate_window["Main"] = self.evaluate_main_window

    def init_node_settings_window(self):
        try:
            del self.node_mask
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
        max_node_count = float(self.evidence.run_count)
        self.node_dimensions["node_evidence"] = [max_node_count, max_node_count]
        self.node_threshold = 2 * len(self.node_dimensions)
        self.node_mask = np.repeat(self.node_threshold, self.ion_network.node_count)
        node_evidence = self.evidence.get_nodes()
        self.node_mask -= node_evidence < max_node_count
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
        self.node_labels = ""
        layout.append(
            [
                sg.Text('NODE LABELS', size=(25, 1)),
                sg.Combo(
                    self.ion_network.dimensions + self.ion_network.ion_comments + ['NODE EVIDENCE', ""],
                    size=(21, 1),
                    default_value=self.node_labels,
                    key="node_labels",
                )
            ]
        )
        self.node_color = "FRAGMENT_LOGINT"
        # self.node_color = "lightgrey"
        self.node_color_c_map = "RdYlGn"
        self.node_color_normalize = False
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
        self.window["Node settings"] = sg.Window(
            'Node settings', layout, finalize=True
        )
        self.window["Node settings"].hide()
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
        self.edges = self.ion_network.get_edges(
            return_as_scipy_csr=True,
            return_pointers=True
        )
        self.positive_edge_evidence = self.evidence.get_edges()
        self.negative_edge_evidence = self.evidence.get_edges(
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
        max_node_count = float(self.evidence.run_count)
        self.edge_formula = "p - n"
        self.min_edge_threshold = max_node_count
        self.max_edge_threshold = max_node_count
        layout.append(
            [
                sg.Text('EDGE EVIDENCE', size=(25, 1)),
                sg.InputText(
                    self.edge_formula,
                    key=f"edge_formula",
                    size=(10, 1),
                ),
                sg.InputText(
                    self.min_edge_threshold,
                    key=f"min_edge_threshold",
                    size=(10, 1),
                ),
                sg.InputText(
                    self.max_edge_threshold,
                    key=f"max_edge_threshold",
                    size=(10, 1),
                ),
            ]
        )
        self.edge_color = "EDGE EVIDENCE"
        # self.edge_color = "lightgrey"
        self.edge_color_c_map = "RdYlGn"
        self.edge_color_normalize = False
        layout.append(
            [
                sg.Text('EDGE COLOR', size=(25, 1)),
                sg.Combo(
                    # ['POSITIVE', 'NEGATIVE', 'SUMMED'] + MATPLOTLIB_COLORS_FIXED,
                    ['EDGE EVIDENCE'] + MATPLOTLIB_COLORS_FIXED,
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
        self.window["Edge settings"] = sg.Window(
            'Edge settings', layout, finalize=True
        )
        self.window["Edge settings"].hide()
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
        self.figs["evidence"].axes[0].set_xticks(
            list(range(self.evidence.run_count + 1))
        )
        self.figs["evidence"].axes[0].set_xticklabels(
            [self.ion_network.run_name] + self.evidence.runs,
            rotation=45,
            ha="right"
        )
        self.figs["evidence"].axes[0].set_xlabel("Run")
        self.figs["evidence"].axes[0].set_xlim(
            [
                0,
                len(self.evidence.runs),
            ]
        )
        flush_figure(self.figs["evidence"], True)
        self.window["Compare settings"] = sg.Window(
            'Compare settings', layout, finalize=True
        )
        self.window["Compare settings"].hide()
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
        self.nodes_on_top = True
        layout.append(
            [
                sg.Checkbox(
                    'NODES ON TOP',
                    default=self.nodes_on_top,
                    key="nodes_on_top"
                ),
            ]
        )
        layout.append(
            [
                sg.Button("Save network plot"),
                sg.Button("Save evidence plot")
            ]
        )
        layout.append(
            [
                # sg.Button("Export network data"),
                sg.Button("Export evidence data")
            ]
        )
        self.figs["network"].axes[0].set_facecolor(self.network_bg_color)
        self.figs["evidence"].axes[0].set_facecolor(self.evidence_bg_color)
        layout.append(self.add_main_menu_and_continue_buttons_to_layout())
        self.window["Misc settings"] = sg.Window(
            'Misc settings', layout, finalize=True
        )
        self.window["Misc settings"].hide()
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
            flush_figure(self.figs["network"], True)
        if self.evidence_bg_color != values["evidence_bg_color"]:
            self.evidence_bg_color = values["evidence_bg_color"]
            self.figs["evidence"].axes[0].set_facecolor(self.evidence_bg_color)
            flush_figure(self.figs["evidence"], True)
        if self.nodes_on_top != values["nodes_on_top"]:
            # TODO
            e = self.edge_collection.get_zorder()
            n = self.node_scatter.get_zorder()
            self.edge_collection.set_zorder(n)
            self.node_scatter.set_zorder(e)
            flush_figure(self.figs["network"], True)
            self.nodes_on_top = values["nodes_on_top"]
        if event == "Save network plot":
            file_name = sg.popup_get_file("Save network plot", save_as=True)
            if file_name is None:
                return
            with loading_window():
                self.figs["network"].savefig(file_name)
        if event == "Save evidence plot":
            file_name = sg.popup_get_file("Save evidence plot", save_as=True)
            if file_name is None:
                return
            with loading_window():
                self.figs["evidence"].savefig(file_name)
        if event == "Export evidence data":
            if hasattr(self, "evidence_plot"):
                file_name = sg.popup_get_file("Save evidence data", save_as=True)
                if file_name is None:
                    return
                with loading_window():
                    data = np.array(
                        [
                            lines[0].get_data()[1] for lines in self.evidence_plot
                        ]
                    )
                    data = pd.DataFrame(
                        data,
                        columns=[self.ion_network.run_name] + self.evidence.runs
                    )
                    data.to_csv(file_name, index=False)

    def evaluate_compare_settings_window(
        self,
        event,
        values,
        update_node_selection=False,
        update_node_colors=False,
    ):
        # if self.evidence_axis != values["evidence_axis"]:
        if True:
            self.evidence_axis = values["evidence_axis"]
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
                    node_count = self.evidence.get_nodes()
                    self.node_mask -= node_count >= low
                    self.node_mask += node_count >= float(values["min_node_count"])
                    self.node_dimensions["node_evidence"][0] = float(
                        values["min_node_count"]
                    )
                if high != float(values["max_node_count"]):
                    update_node_selection = True
                    node_count = self.evidence.get_nodes()
                    self.node_mask -= node_count <= high
                    self.node_mask += node_count <= float(values["max_node_count"])
                    self.node_dimensions["node_evidence"][1] = float(
                        values["max_node_count"]
                    )
            else:
                if low != float(values[f"min_{key}"]):
                    update_node_selection = True
                    coordinates = self.ion_network.get_ion_coordinates(key)
                    self.node_mask -= coordinates >= low
                    self.node_mask += coordinates >= float(values[f"min_{key}"])
                    self.node_dimensions[key][0] = float(values[f"min_{key}"])
                if high != float(values[f"max_{key}"]):
                    update_node_selection = True
                    coordinates = self.ion_network.get_ion_coordinates(key)
                    self.node_mask -= coordinates <= high
                    self.node_mask += coordinates <= float(values[f"max_{key}"])
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
        if self.node_labels != values["node_labels"]:
            self.node_labels = values["node_labels"]
            try:
                self.network_figure_update_node_labels(nodes)
            except NameError:
                self.network_figure_update_node_labels()

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
        if self.min_edge_threshold != values["min_edge_threshold"]:
            self.min_edge_threshold = float(
                values["min_edge_threshold"]
            )
            update_edge_selection = True
        if self.max_edge_threshold != values["max_edge_threshold"]:
            self.max_edge_threshold = float(
                values["max_edge_threshold"]
            )
            update_edge_selection = True
        if self.edge_formula != values["edge_formula"]:
            self.edge_formula = values["edge_formula"]
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
        reorder=False,
        flush=True,
    ):
        # TODO: Docstring
        if nodes is None:
            nodes = self.get_filtered_nodes()
        if len(nodes) == 0:
            flush_figure(self.figs["network"], flush)
            return
        if self.node_color in self.ion_network.dimensions + ['NODE EVIDENCE']:
            if self.node_color in self.ion_network.dimensions:
                inds = self.ion_network.get_ion_coordinates(self.node_color)
            elif self.node_color == 'NODE EVIDENCE':
                inds = self.evidence.get_nodes()
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
        flush_figure(self.figs["network"], flush)
        self.reset_selected_nodes()

    def network_figure_update_node_labels(
        self,
        nodes=None,
        reorder=False,
        flush=True
    ):
        if nodes is None:
            nodes = self.get_filtered_nodes()
        if len(nodes) == 0:
            flush_figure(self.figs["network"], flush)
            return
        # TODO:
        for child in self.figs["network"].axes[0].get_children():
            if isinstance(
                child, matplotlib.text.Annotation
            ):
                child.remove()
        if self.node_labels in self.ion_network.dimensions:
            labels = np.round(
                self.ion_network.get_ion_coordinates(
                    self.node_labels,
                    nodes
                ),
                NODE_LABEL_DECIMALS
            )
        elif self.node_labels in self.ion_network.ion_comments:
            labels = self.ion_network.get_ion_comments(
                self.node_labels,
                nodes
            )
        elif self.node_labels == 'NODE EVIDENCE':
            labels = self.evidence.get_nodes()[nodes]
        else:
            flush_figure(self.figs["network"], flush)
            return
        x_coordinates, y_coordinates = self.ion_network.get_ion_coordinates(
            [self.x_axis, self.y_axis],
            nodes
        )
        for x, y, label in zip(
            x_coordinates,
            y_coordinates,
            labels
        ):
            self.figs["network"].axes[0].annotate(label, (x, y))
        flush_figure(self.figs["network"], flush)

    def network_figure_update_edge_colors(
        self,
        selected_edges=None,
        reorder=True,
        flush=True,
    ):
        # TODO: Docstring
        if self.show_edges:
            if selected_edges is None:
                nodes = self.get_filtered_nodes()
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
            if len(selected_edges) == 0:
                flush_figure(self.figs["network"], flush)
                return
            if self.edge_color in ['EDGE EVIDENCE']:
                inds = ne.evaluate(
                    self.edge_formula,
                    local_dict={
                        "p": self.positive_edge_evidence,
                        "n": self.negative_edge_evidence
                    },
                    global_dict={},
                )
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
        flush_figure(self.figs["network"], flush)

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
                nodes = self.get_filtered_nodes()
            if (x_coordinates is None) or (y_coordinates is None):
                x_coordinates, y_coordinates = self.ion_network.get_ion_coordinates(
                    [self.x_axis, self.y_axis],
                    nodes
                )
            selected_neighbors = self.edges[nodes].T.tocsr()[nodes]
            a, b = selected_neighbors.nonzero()
            try:
                values = ne.evaluate(
                    self.edge_formula,
                    local_dict={
                        "p": self.positive_edge_evidence[
                            selected_neighbors.data
                        ],
                        "n": self.negative_edge_evidence[
                            selected_neighbors.data
                        ]
                    },
                    global_dict={},
                )
                selection = (values >= self.min_edge_threshold)
                selection &= (values <= self.max_edge_threshold)
                a = a[selection]
                b = b[selection]
                start_edges = list(zip(x_coordinates[a], y_coordinates[a]))
                end_edges = list(zip(x_coordinates[b], y_coordinates[b]))
                edges = np.array(list(zip(start_edges, end_edges)))
                selected_edges = selected_neighbors.data[selection]
            except NotImplementedError:
                sg.popup_error(
                    f"The edge formula '{self.edge_formula}' is not valid"
                )
                edges = []
                selected_edges = []
            except KeyError:
                sg.popup_error(
                    f"The edge formula '{self.edge_formula}' is not valid"
                )
                edges = []
                selected_edges = []
        else:
            edges = []
            selected_edges = []
        self.edge_collection.set_segments(edges)
        flush_figure(self.figs["network"], flush)
        return selected_edges

    def network_figure_update_node_selection(self, flush=True):
        # TODO: Docstring
        nodes = self.get_filtered_nodes()
        x_coordinates, y_coordinates = self.ion_network.get_ion_coordinates(
            [self.x_axis, self.y_axis],
            nodes
        )
        if not hasattr(self, "node_scatter"):
            self.node_scatter = self.figs["network"].axes[0].scatter(
                x_coordinates,
                y_coordinates,
                marker=".",
                picker=True,
                zorder=4
            )
        else:
            self.node_scatter.set_offsets(np.c_[x_coordinates, y_coordinates])
        return nodes, x_coordinates, y_coordinates

    def evidence_figure_update_axis_selection(self, flush=True):
        # TODO: Docstring
        self.figs["evidence"].axes[0].set_ylabel(self.evidence_axis)
        node_mask = np.zeros(self.ion_network.node_count, np.int64)
        nodes = self.get_selected_nodes()
        if len(nodes) == 0:
            flush_figure(self.figs["evidence"], flush)
            return
        node_mask[nodes] = np.arange(len(nodes))
        alignments = np.zeros((len(nodes), self.evidence.run_count + 1))
        alignments[:, 0] = self.ion_network.get_ion_coordinates(
            self.evidence_axis,
            nodes
        )
        for i, other_evidence in enumerate(self.evidence):
            self_alignment = self.evidence.get_nodes(other_evidence)
            other_alignment = other_evidence.get_nodes(self.evidence)
            selected = np.isin(self_alignment, nodes)
            self_ions = node_mask[self_alignment[selected]]
            other_ions = other_alignment[selected]
            alignments[self_ions, i + 1] = other_evidence.ion_network.get_ion_coordinates(
                self.evidence_axis,
                other_ions
            )
            if self.evidence_axis == "FRAGMENT_LOGINT":
                alignments[self_ions, i + 1] += self.evidence.read_attr(
                    "median_intensity_correction",
                    parent_group_name=f"runs/{other_evidence.run_name}/nodes"
                )
        if hasattr(self, "evidence_plot"):
            for evidence_plot in self.evidence_plot:
                for i in evidence_plot:
                    i.remove()
                del evidence_plot[:]
            del self.evidence_plot[:]
        color_mapper = matplotlib.cm.ScalarMappable(
            norm=matplotlib.colors.Normalize(
                vmin=0,
                vmax=20,
            ),
            cmap="tab20"
        )
        colors = color_mapper.to_rgba(nodes % 20)
        self.evidence_plot = []
        for i, c in zip(alignments, colors):
            evidence_plot = self.figs["evidence"].axes[0].plot(
                i,
                c=c
            )
            self.evidence_plot.append(evidence_plot)

        self.figs["evidence"].axes[0].set_ylim(
            [
                np.min(alignments),
                np.max(alignments),
            ]
        )
        flush_figure(self.figs["evidence"], flush)

    def init_figure(self, title):
        fig = plt.figure()
        fig.add_subplot(111)
        fig.canvas.manager.toolmanager.remove_tool('back')
        fig.canvas.manager.toolmanager.remove_tool('forward')
        fig.canvas.manager.toolmanager.remove_tool('home')
        try:
            fig.canvas.manager.toolmanager.remove_tool('copy')
        except KeyError:
            pass
        try:
            fig.canvas.manager.toolmanager.remove_tool('allnav')
        except KeyError:
            pass
        try:
            fig.canvas.manager.toolmanager.remove_tool('nav')
        except KeyError:
            pass
        fig.canvas.manager.window.protocol(
            "WM_DELETE_WINDOW",
            lambda: self.swap_active_window(None)
        )
        fig.canvas.set_window_title(title)
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
        flush_figure(fig, flush)

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
        flush_figure(fig, flush)

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

    def get_filtered_nodes(self):
        return np.flatnonzero(self.node_mask == self.node_threshold)

    def get_selected_nodes(self):
        nodes = self.get_filtered_nodes()
        node_selection = self.figs[
            "network"
        ].canvas.manager.toolmanager._tools['node_select']._current_selection
        return nodes[sorted(node_selection)]

    def reset_selected_nodes(self):
        self.figs["network"].canvas.manager.toolmanager._tools['node_select'].hard_reset()

    def init_network(self):
        file_name = sg.popup_get_file(
            'Please select an ion-network',
            file_types=(('Ion-network', '*.inet.hdf'),)
        )
        if file_name is None:
            return
        try:
            ion_network = ms_run_files.HDF_Network_File(file_name)
        except (OSError, ValueError):
            sg.popup_error('This is not a valid ion_network')
            return
        try:
            evidence = ms_run_files.HDF_Evidence_File(file_name)
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


def flush_figure(figure, flush):
    if flush:
        figure.canvas.draw()
        figure.canvas.flush_events()


def bind(instance, func, as_name=None):
    """
    Bind the function *func* to *instance*, with either provided name *as_name*
    or the existing name of *func*. The provided *func* should accept the
    instance as the first argument, i.e. "self".
    """
    if as_name is None:
        as_name = func.__name__
    bound_method = func.__get__(instance, instance.__class__)
    setattr(instance, as_name, bound_method)
    return bound_method


class PointerTool(matplotlib.backend_tools.ToolToggleBase):
    description = 'Select/deselected nodes with left/right mouseclick, single/double middle mouseclick to select none/all, scroll to undo/redo'
    image = os.path.join(DEFAULT_BROWSER_PATH, DEFAULT_BROWSER_IMAGES["pointer"])
    radio_group = 'default'
    default_keymap = 's'

    def __init__(self, *args, callback_function=None):
        super().__init__(*args)
        self._pick_connection = None
        self._zoom_connection = None
        self._button_connection = None
        self._modifiable = False
        self._stack = []
        self._stack_pointer = 0
        self._current_selection = set()
        self._scatter_to_update = False
        self.callback_function = callback_function

    def enable(self, event):
        self._pick_connection = self.figure.canvas.mpl_connect(
            'pick_event',
            self.pick_event
        )
        self._zoom_connection = self.figure.canvas.mpl_connect(
            'scroll_event',
            self.scroll_event
        )
        self._button_connection = self.figure.canvas.mpl_connect(
            'button_press_event',
            self.button_event
        )
        self._modifiable = True

    def disable(self, event):
        self._modifiable = False
        self.figure.canvas.mpl_disconnect(self._button_connection)
        self.figure.canvas.mpl_disconnect(self._zoom_connection)
        self.figure.canvas.mpl_disconnect(self._pick_connection)

    def pick_event(self, event):
        if self._modifiable:
            self._modifiable = False
            for index in event.ind:
                if event.mouseevent.button == 1:
                    self.update_stack(index + 1)
                if event.mouseevent.button == 3:
                    self.update_stack(-index - 1)
            self._modifiable = True
            self.update_scatter()

    def scroll_event(self, event):
        if self._modifiable:
            self._modifiable = False
            if event.button == "up":
                self.update_stack_pointer(direction="up")
            elif event.button == "down":
                self.update_stack_pointer(direction="down")
            self._modifiable = True
            self.update_scatter()

    def button_event(self, event):
        if self._modifiable:
            self._modifiable = False
            if event.button == 2:
                if event.dblclick:
                    self.set_all(mode="select")
                else:
                    self.set_all(mode="delete")
            self._modifiable = True
            self.update_scatter()

    def hard_reset(self):
        self._stack = []
        self._stack_pointer = 0
        self._current_selection = set()
        self._scatter_to_update = True
        self.update_scatter()

    def set_all(self, mode="delete"):
        if mode == "delete":
            for index in list(self._current_selection):
                self.update_stack(-(index - 1))
        elif mode == "select":
            try:
                scatter_size = self.figure.axes[0].collections[0]._offsets.shape[0]
            except IndexError:
                return
            for index in range(1, scatter_size + 1):
                self.update_stack(index)
        self.update_scatter()

    def update_stack(self, index):
        if (index > 0) and (index - 1 in self._current_selection):
            return
        if (index < 0) and (-index - 1 not in self._current_selection):
            return
        self._stack = self._stack[:self._stack_pointer]
        self._stack.append(index)
        self.update_stack_pointer("up")

    def update_stack_pointer(self, direction="up"):
        if direction == "up":
            if self._stack_pointer < len(self._stack):
                new_element = self._stack[self._stack_pointer]
                self.update_current_selection(new_element)
                self._stack_pointer += 1
        elif direction == "down":
            if self._stack_pointer > 0:
                self._stack_pointer -= 1
                new_element = self._stack[self._stack_pointer]
                self.update_current_selection(-new_element)

    def update_current_selection(self, index):
        if index > 0:
            self._current_selection.add(index - 1)
        elif index < 0:
            self._current_selection.remove(-index - 1)
        self._scatter_to_update = True

    def update_scatter(self):
        if self._scatter_to_update:
            scatter = self.figure.axes[0].collections[0]
            edge_colors = scatter.get_facecolors().copy()
            try:
                edge_colors[list(self._current_selection)] = [0, 0, 0, 1]
            except IndexError:
                edge_colors = np.repeat(edge_colors, scatter._offsets.shape[0], axis=0)
                edge_colors[list(self._current_selection)] = [0, 0, 0, 1]
            scatter.set_edgecolors(edge_colors)
            flush_figure(self.figure, True)
            if self.callback_function is not None:
                self.callback_function()
            self._scatter_to_update = False
