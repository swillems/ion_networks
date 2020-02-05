#!python

import numpy as np
import matplotlib.pyplot as plt
import PySimpleGUI as sg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.ticker import NullFormatter
import matplotlib
matplotlib.use('TkAgg')


class GUI(object):

    def __init__(self, ion_network, evidence, logger):
        self.ion_network = ion_network
        self.evidence = evidence
        self.create_overview_window()
        self.create_plot_window()
        while True:
            if self.update_window(self.overview_window) is not None:
                break
            # if self.update_window(self.plot_window) is not None:
            #     break
        self.overview_window.close()
        self.plot_window.close()

    def create_overview_window(self):
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
                    size=(10, 1)
                ),
                sg.InputText(
                    max_coordinate,
                    key=f"max_{dimension}",
                    size=(10, 1)
                ),
            ]
            layout.append(row)
        max_node_count = float(self.evidence.network_count)
        self.node_dimensions["node_evidence"] = [max_node_count, max_node_count]
        self.node_threshold = 2 * len(self.node_dimensions)
        self.nodes = np.repeat(self.node_threshold , self.ion_network.node_count)
        node_evidence = self.evidence.get_evidence(kind="nodes", return_total=True)
        self.nodes -= node_evidence < max_node_count
        layout.append(
            [
                sg.Text('NODE EVIDENCE', size=(25, 1)),
                sg.InputText(
                    max_node_count,
                    key=f"min_node_count",
                    size=(10, 1)
                ),
                sg.InputText(
                    max_node_count,
                    key=f"max_node_count",
                    size=(10, 1)
                ),
            ]
        )
        layout.append([sg.Button('Refresh')])
        self.overview_window = sg.Window('Settings', layout)

    def create_plot_window(self):
        plt.figure(1)
        plt.subplot(111)
        self.fig = plt.gcf()      # if using Pyplot then get the figure from the plot
        figure_x, figure_y, figure_w, figure_h = self.fig.bbox.bounds
        layout = [
            [sg.Canvas(size=(figure_w, figure_h), key='canvas')]
        ]
        self.plot_window = sg.Window('Plot', layout, finalize=True)
        self.figure_canvas_agg = FigureCanvasTkAgg(self.fig, self.plot_window['canvas'].TKCanvas)
        self.figure_canvas_agg.get_tk_widget().pack(side='top', fill='both', expand=1)
        self.update_plot()

    def update_plot(self):
        x_coordinates, y_coordinates = self.ion_network.get_ion_coordinates(
            ["RT", "MZ1"],
            np.flatnonzero(self.nodes == self.node_threshold)
        )
        if not hasattr(self, "old_scatter"):
            self.old_scatter = plt.scatter(x_coordinates, y_coordinates, marker=".")
        elif not hasattr(self, "new_scatter"):
            self.new_scatter = plt.scatter(x_coordinates, y_coordinates, marker=".")
        else:
            self.old_scatter.set_offsets(self.new_scatter.get_offsets())
            self.new_scatter.set_offsets(np.c_[x_coordinates, y_coordinates])
        self.figure_canvas_agg.draw()

    def update_window(self, window):
        event, values = window.read(timeout=100)
        if event == sg.TIMEOUT_KEY:
            return
        if event in (None, 'Exit'):
            return "exit"
        if window.Title == "Settings":
            self.perform_overview_action(event, values)
        if window.Title == "Plot":
            self.perform_plot_action(event, values)

    def perform_overview_action(self, event, values):
        update = False
        for key, (low, high) in list(self.node_dimensions.items()):
            if key == "node_evidence":
                if low != float(values["min_node_count"]):
                    update = True
                    node_count = self.evidence.get_evidence(kind="nodes", return_total=True)
                    self.nodes -= node_count >= low
                    self.nodes += node_count >= float(values["min_node_count"])
                    self.node_dimensions["node_evidence"][0] = float(values["min_node_count"])
                if high != float(values["max_node_count"]):
                    update = True
                    node_count = self.evidence.get_evidence(kind="nodes", return_total=True)
                    self.nodes -= node_count <= high
                    self.nodes += node_count <= float(values["max_node_count"])
                    self.node_dimensions["node_evidence"][1] = float(values["max_node_count"])
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
            print(np.bincount(self.nodes == self.node_threshold))
            self.update_plot()

    def perform_plot_action(self, event, values):
        print("plot", event, values)
