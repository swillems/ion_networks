#!python

# external
import PySimpleGUI as sg

# local
import interface


class GUI(object):
    # TODO: Docstring

    def __init__(
        self,
        start=True,
        widget_size=25,
    ):
        # TODO: Docstring
        self.widget_size = widget_size
        self.windows = {}
        self.init_main_window()
        if start:
            self.run()

    def init_main_window(self):
        # TODO: Docstring
        self.windows["Main"] = Window("Ion-networks")
        self.windows["Main"].layout.append(
            [sg.Button("Convert", size=(self.widget_size, 1))]
        )
        self.windows["Main"].layout.append(
            [sg.Button("Create", size=(self.widget_size, 1))]
        )
        self.windows["Main"].layout.append(
            [sg.Button("Evidence", size=(self.widget_size, 1))]
        )
        self.windows["Main"].layout.append(
            [sg.Button("Show", size=(self.widget_size, 1))]
        )
        self.windows["Main"].layout.append(
            [sg.Button("Exit", size=(self.widget_size, 1))]
        )
        self.windows["Main"].activate()

    def init_convert_window(self):
        # TODO: Docstring
        window = Window("Convert")
        window.add_input_path_to_layout()
        window.add_output_directory_to_layout()
        window.layout.append(
            [
                sg.Text('Data type', size=(self.widget_size, 1)),
                sg.Combo(
                    ['DDA', 'HDMSE', "SONAR", "SWIMDIA"],
                    default_value='HDMSE',
                    key="data_type",
                    # enable_events=True
                )
            ]
        )
        window.add_parameter_file_to_layout()
        window.add_log_file_to_layout()
        window.add_continue_and_cancel_buttons_to_layout()
        self.windows["Convert"] = window

    def init_create_window(self):
        # TODO: Docstring
        window = Window("Create")
        window.add_input_path_to_layout()
        window.add_output_directory_to_layout()
        window.add_parameter_file_to_layout()
        window.add_log_file_to_layout()
        window.add_continue_and_cancel_buttons_to_layout()
        self.windows["Create"] = window

    def init_evidence_window(self):
        # TODO: Docstring
        window = Window("Evidence")
        window.add_input_path_to_layout()
        window.add_output_directory_to_layout()
        window.add_parameter_file_to_layout()
        window.add_log_file_to_layout()
        window.add_continue_and_cancel_buttons_to_layout()
        self.windows["Evidence"] = window

    def init_show_window(self):
        # TODO: Docstring
        window = Window("Show")
        # window.add_input_path_to_layout()
        # window.add_output_directory_to_layout()
        # window.add_parameter_file_to_layout()
        # window.add_log_file_to_layout()
        window.add_continue_and_cancel_buttons_to_layout()
        self.windows["Show"] = window

    def run(self):
        # TODO: Docstring
        finished = None
        while finished is None:
            for window_name, window in list(self.windows.items()):
                if not window.active:
                    continue
                event, values = window.read(timeout=100)
                if event == sg.TIMEOUT_KEY:
                    continue
                elif window_name == "Main":
                    finished = self.evaluate_main_window(event, values)
                elif window_name == "Convert":
                    self.evaluate_convert_window(event, values)
                elif window_name == "Create":
                    self.evaluate_create_window(event, values)
                elif window_name == "Evidence":
                    self.evaluate_evidence_window(event, values)
                elif window_name == "Show":
                    self.evaluate_show_window(event, values)

    def evaluate_main_window(self, event, values):
        # TODO: Docstring
        if event in (None, 'Exit'):
            return "Exit"
        if event in ["Convert", "Create", "Evidence", "Show"]:
            self.windows["Main"].deactivate(delete=False)
        if event == "Convert":
            self.init_convert_window()
        elif event == "Create":
            self.init_create_window()
        elif event == "Evidence":
            self.init_evidence_window()
        elif event == "Show":
            self.init_show_window()
        self.windows[event].activate()

    def evaluate_convert_window(self, event, values):
        # TODO: Docstring
        if event in (None, 'Cancel'):
            self.windows["Convert"].deactivate()
            self.windows["Main"].activate()
        elif event == "Continue":
            interface.convert_data_formats_to_csvs(
                [values["input_path"]],
                values["data_type"],
                values["output_directory"],
                values["parameter_file_name"],
                values["log_file_name"]
            )

    def evaluate_create_window(self, event, values):
        # TODO: Docstring
        if event in (None, 'Cancel'):
            self.windows["Create"].deactivate()
            self.windows["Main"].activate()
        elif event == "Continue":
            interface.create_ion_networks(
                [values["input_path"]],
                values["output_directory"],
                values["parameter_file_name"],
                values["log_file_name"]
            )

    def evaluate_evidence_window(self, event, values):
        # TODO: Docstring
        if event in (None, 'Cancel'):
            self.windows["Evidence"].deactivate()
            self.windows["Main"].activate()
        elif event == "Continue":
            interface.evidence_ion_networks(
                [values["input_path"]],
                values["output_directory"],
                values["parameter_file_name"],
                values["log_file_name"]
            )

    def evaluate_show_window(self, event, values):
        # TODO: Docstring
        if event in (None, 'Cancel'):
            self.windows["Show"].deactivate()
            self.windows["Main"].activate()
        elif event == "Continue":
            interface.show_ion_network(
                [values["ion_network_file_name"]],
                values["evidence_file_name"],
                values["parameter_file_name"],
                values["log_file_name"]
            )


class Window(object):
    # TODO: Docstring

    def __init__(
        self,
        name="",
        widget_size=25,
    ):
        # TODO: Docstring
        self.name = name
        self.layout = []
        self.active = False
        self.widget_size = widget_size

    def read(self, *args, **kwargs):
        # TODO: Docstring
        return self.window.read(*args, **kwargs)

    def activate(self):
        # TODO: Docstring
        self.active = True
        if not hasattr(self, "window"):
            self.window = sg.Window(self.name, self.layout)
        # self.window.Enable()

    def deactivate(self, delete=True):
        # TODO: Docstring
        self.active = False
        # self.window.Disable()
        if delete:
            self.window.Close()

    def add_input_path_to_layout(self):
        # TODO: Multiple and independent files?
        self.layout.append(
            [
                sg.Text("Input path", size=(self.widget_size, 1)),
                sg.Input(key="input_path"),
                sg.FolderBrowse(),
            ]
        )

    def add_output_directory_to_layout(self):
        # TODO: Multiple and independent files?
        self.layout.append(
            [
                sg.Text("Output directory", size=(self.widget_size, 1)),
                sg.Input(key="output_directory"),
                sg.FolderBrowse(),
            ]
        )

    def add_parameter_file_to_layout(self):
        self.layout.append(
            [
                sg.Text("Parameter file", size=(self.widget_size, 1)),
                sg.Input(key="parameter_file_name"),
                sg.FileBrowse(),
            ]
        )

    def add_log_file_to_layout(self):
        self.layout.append(
            [
                sg.Text("Log file", size=(self.widget_size, 1)),
                sg.Input(key="log_file_name"),
                sg.FileBrowse(),
            ]
        )

    def add_continue_and_cancel_buttons_to_layout(self):
        self.layout.append(
            [
                sg.Button("Cancel", size=(self.widget_size, 1)),
                sg.Button("Continue", size=(self.widget_size, 1))
            ]
        )
