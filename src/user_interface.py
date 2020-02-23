#!python

# builtin
import threading
# external
import PySimpleGUI as sg
# local
import interface


class GUI(object):
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
        self.init_convert_window()
        self.init_create_window()
        self.init_evidence_window()
        self.init_show_window()
        self.init_terminal_window()
        self.window["Main"].activate()
        self.active_window = "Main"
        if start:
            self.run()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        for window in list(self.window.values()):
            window.close()

    def init_main_window(self):
        # TODO: Docstring
        window = Window("Ion-networks", widget_size=self.widget_size)
        window.layout = [
            [sg.Button("Convert", size=(self.widget_size, 1))],
            [sg.Button("Create", size=(self.widget_size, 1))],
            [sg.Button("Evidence", size=(self.widget_size, 1))],
            [sg.Button("Show", size=(self.widget_size, 1))],
        ]
        self.window["Main"] = window
        self.evaluate_window["Main"] = self.evaluate_main_window

    def init_convert_window(self):
        # TODO: Docstring
        window = Window("Convert", widget_size=self.widget_size)
        window.add_input_path_to_layout()
        window.add_output_directory_to_layout()
        window.layout.append(
            [
                sg.Text('Data type', size=(self.widget_size, 1)),
                sg.Combo(
                    ['DDA', 'HDMSE', "SONAR", "SWIMDIA"],
                    default_value='HDMSE',
                    key="data_type",
                    size=(self.widget_size * 2, 1)
                    # enable_events=True
                )
            ]
        )
        window.add_parameter_file_to_layout()
        window.add_log_file_to_layout()
        window.add_continue_and_cancel_buttons_to_layout()
        self.window["Convert"] = window
        self.evaluate_window["Convert"] = self.evaluate_convert_window

    def init_create_window(self):
        # TODO: Docstring
        window = Window("Create", widget_size=self.widget_size)
        window.add_input_path_to_layout()
        window.add_output_directory_to_layout()
        window.add_parameter_file_to_layout()
        window.add_log_file_to_layout()
        window.add_continue_and_cancel_buttons_to_layout()
        self.window["Create"] = window
        self.evaluate_window["Create"] = self.evaluate_create_window

    def init_evidence_window(self):
        # TODO: Docstring
        window = Window("Evidence", widget_size=self.widget_size)
        window.add_input_path_to_layout()
        window.add_output_directory_to_layout()
        window.add_parameter_file_to_layout()
        window.add_log_file_to_layout()
        window.add_continue_and_cancel_buttons_to_layout()
        self.window["Evidence"] = window
        self.evaluate_window["Evidence"] = self.evaluate_evidence_window

    def init_terminal_window(self):
        # TODO: Docstring
        window = Window("Terminal", widget_size=self.widget_size)
        window.layout.append([sg.Output(size=(150, 50))])
        window.add_continue_and_cancel_buttons_to_layout(continue_button=False)
        self.window["Terminal"] = window

    def init_show_window(self):
        # TODO: Docstring
        window = Window("Show", widget_size=self.widget_size)
        # window.add_input_path_to_layout()
        # window.add_output_directory_to_layout()
        # window.add_parameter_file_to_layout()
        # window.add_log_file_to_layout()
        window.add_continue_and_cancel_buttons_to_layout()
        self.window["Show"] = window
        self.evaluate_window["Show"] = self.evaluate_show_window

    def run(self):
        # TODO: Docstring
        keep_running = True
        while keep_running:
            for window_name, window in list(self.window.items()):
                if not window.active:
                    continue
                event, values = window.read(timeout=10)
                if event == sg.TIMEOUT_KEY:
                    continue
                elif event is None:
                    keep_running = False
                    break
                elif event == "Return to main menu":
                    self.swap_active_window("Main")
                else:
                    self.evaluate_window[window_name](event, values)

    def evaluate_main_window(self, event, values):
        # TODO: Docstring
        self.swap_active_window(event)

    def evaluate_convert_window(self, event, values):
        # TODO: Docstring
        if event == "Continue":
            self.run_terminal_command(
                interface.convert_data_formats_to_csvs,
                [values["input_path"]],
                values["data_type"],
                values["output_directory"],
                values["parameter_file_name"],
                values["log_file_name"]
            )

    def evaluate_create_window(self, event, values):
        # TODO: Docstring
        if event == "Continue":
            self.run_terminal_command(
                interface.create_ion_networks,
                [values["input_path"]],
                values["output_directory"],
                values["parameter_file_name"],
                values["log_file_name"]
            )

    def evaluate_evidence_window(self, event, values):
        # TODO: Docstring
        if event == "Continue":
            self.run_terminal_command(
                interface.evidence_ion_networks,
                [values["input_path"]],
                values["output_directory"],
                values["parameter_file_name"],
                values["log_file_name"]
            )

    def evaluate_show_window(self, event, values):
        # TODO: Docstring
        if event == "Continue":
            self.run_terminal_command(
                interface.show_ion_network,
                [values["ion_network_file_name"]],
                values["evidence_file_name"],
                values["parameter_file_name"],
                values["log_file_name"]
            )

    def run_terminal_command(self, command, *args):
        # TODO: Docstring
        self.swap_active_window("Terminal")
        self.window["Terminal"].read(timeout=10)
        self.window["Terminal"].window["Return to main menu"].Update(
            text="Executing, please wait",
            disabled=True
        )
        thread = threading.Thread(target=command, args=args)
        thread.start()
        thread.deamon = True
        while thread.isAlive():
            event, values = self.window["Terminal"].read(timeout=10)
            if event is None:
                sg.Popup(
                    "WARNING, thread is still running in the background!",
                    title="WARNING",
                )
        self.window["Terminal"].window["Return to main menu"].Update(
            text="Return to main menu",
            disabled=False
        )

    def swap_active_window(self, new_window_name):
        # TODO: Docstring, implement
        self.window[self.active_window].deactivate()
        self.window[new_window_name].activate()
        self.active_window = new_window_name


class Window(object):
    # TODO: Docstring

    def __init__(
        self,
        name="",
        widget_size=20,
    ):
        # TODO: Docstring
        self.name = name
        self.layout = []
        self.active = False
        self.widget_size = widget_size

    def read(self, *args, **kwargs):
        # TODO: Docstring
        return self.window.read(*args, **kwargs)

    def init(self):
        # TODO: Docstring
        # TODO: default_element_size?, default_button_element_size
        self.window = sg.Window(self.name, self.layout)

    def close(self):
        if hasattr(self, "window"):
            self.window.close()

    def activate(self):
        # TODO: Docstring
        self.active = True
        if not hasattr(self, "window"):
            self.init()
        self.window.UnHide()

    def deactivate(self):
        # TODO: Docstring
        self.active = False
        self.window.Hide()

    def add_input_path_to_layout(self):
        # TODO: Docstring
        # TODO: Multiple and independent files?
        self.layout.append(
            [
                sg.Text("Input path", size=(self.widget_size, 1)),
                sg.Input(
                    key="input_path",
                    size=(self.widget_size * 2, 1)
                ),
                sg.FolderBrowse(size=(self.widget_size, 1)),
            ]
        )

    def add_output_directory_to_layout(self):
        # TODO: Docstring
        self.layout.append(
            [
                sg.Text("Output directory", size=(self.widget_size, 1)),
                sg.Input(
                    key="output_directory",
                    size=(self.widget_size * 2, 1)
                ),
                sg.FolderBrowse(size=(self.widget_size, 1)),
            ]
        )

    def add_parameter_file_to_layout(self):
        # TODO: Docstring
        self.layout.append(
            [
                sg.Text("Parameter file", size=(self.widget_size, 1)),
                sg.Input(
                    key="parameter_file_name",
                    size=(self.widget_size * 2, 1)
                ),
                sg.FileBrowse(size=(self.widget_size, 1)),
            ]
        )

    def add_log_file_to_layout(self):
        self.layout.append(
            [
                sg.Text("Log file", size=(self.widget_size, 1)),
                sg.Input(
                    key="log_file_name",
                    size=(self.widget_size * 2, 1)
                ),
                sg.FileBrowse(size=(self.widget_size, 1)),
            ]
        )

    def add_continue_and_cancel_buttons_to_layout(
        self,
        cancel_button=True,
        continue_button=True
    ):
        row = []
        if cancel_button:
            row.append(
                sg.Button("Return to main menu", size=(self.widget_size, 1))
            )
        if continue_button:
            row.append(sg.Button("Continue", size=(self.widget_size, 1)))
        self.layout.append(row)
