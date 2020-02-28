#!python

# builtin
import os
import threading
# external
import PySimpleGUI as sg
# local
import network
import evidence
import browser
import utils


def convert_data_formats_to_csvs(
    input_path,
    data_type,
    output_directory,
    parameter_file_name,
    log_file_name
):
    """
    Convert centroided MSMS data to a unified csv that can be read as an
    ion-network.

    Parameters
    ----------
    input_path : iterable[str]
        An iterable with file and/or folder names.
    output_directory : str or None
        If provided, all new files will be saved in this directory.
    data_type : str
        The data type of the input files. Options are:
            'dda'
            'sonar'
            'hdmse'
            'swimdia'
    parameter_file_name : str or None
        If provided, parameters will be read from this file.
    log_file_name : str or None
        If provided, all logs will be written to this file.
    """
    if parameter_file_name is None:
        parameter_file_name = ""
    if output_directory is None:
        output_directory = ""
    if log_file_name is None:
        log_file_name = ""
    parameters = utils.read_parameters_from_json_file(
        file_name=parameter_file_name
    )
    with utils.open_logger(log_file_name, parameters=parameters) as logger:
        # logger.info("Running command: ")
        # logger.info("")
        if output_directory != "":
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)
        input_file_names = utils.get_file_names_with_extension(
            input_path,
            extension=utils.DATA_TYPE_FILE_EXTENSIONS[data_type]
        )
        file_count = len(input_file_names)
        logger.info(f"Found {file_count} files to process.")
        for input_file_name in sorted(input_file_names):
            data = utils.read_data_from_file(data_type, input_file_name, logger)
            file_name_base = os.path.splitext(
                os.path.basename(input_file_name)
            )[0]
            if output_directory == "":
                output_path = os.path.dirname(input_file_name)
            else:
                output_path = output_directory
            output_file_name = os.path.join(
                output_path,
                f"{file_name_base}.inet.csv"
            )
            utils.write_data_to_csv_file(data, output_file_name, logger)


def create_ion_networks(
    input_path,
    output_directory,
    parameter_file_name,
    log_file_name
):
    """
    Create ion-networks from unified csv files.

    Parameters
    ----------
    input_path : iterable[str]
        An iterable with file and/or folder names.
    output_directory : str or None
        If provided, all new files will be saved in this directory.
    parameter_file_name : str or None
        If provided, parameters will be read from this file.
    log_file_name : str or None
        If provided, all logs will be written to this file.
    """
    if parameter_file_name is None:
        parameter_file_name = ""
    if output_directory is None:
        output_directory = ""
    if log_file_name is None:
        log_file_name = ""
    parameters = utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
        default="create"
    )
    with utils.open_logger(log_file_name, parameters=parameters) as logger:
        input_file_names = utils.get_file_names_with_extension(
            input_path,
            ".inet.csv"
        )
        file_count = len(input_file_names)
        logger.info(f"Found {file_count} .inet.csv files to process.")
        for csv_file_name in input_file_names:
            local_file_name = os.path.basename(csv_file_name)
            if output_directory == "":
                output_path = os.path.dirname(csv_file_name)
            else:
                output_path = output_directory
            ion_network_file_name = os.path.join(
                output_path,
                f"{local_file_name[:-9]}.inet.hdf"
            )
            network.Network(
                network_file_name=ion_network_file_name,
                centroided_csv_file_name=csv_file_name,
                parameters=parameters,
                logger=logger
            )


def evidence_ion_networks(
    input_path,
    output_directory,
    parameter_file_name,
    log_file_name
):
    """
    Evidence ion-networks with each other.

    Parameters
    ----------
    input_path : iterable[str]
        An iterable with file and/or folder names.
    output_directory : str or None
        If provided, all new files will be saved in this directory.
    parameter_file_name : str or None
        If provided, parameters will be read from this file.
    log_file_name : str or None
        If provided, all logs will be written to this file.
    """
    if parameter_file_name is None:
        parameter_file_name = ""
    if output_directory is None:
        output_directory = ""
    if log_file_name is None:
        log_file_name = ""
    parameters = utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
        default="evidence"
    )
    with utils.open_logger(log_file_name, parameters=parameters) as logger:
        input_file_names = utils.get_file_names_with_extension(
            input_path,
            ".inet.hdf"
        )
        file_count = len(input_file_names)
        logger.info(f"Found {file_count} .inet.hdf files to process.")
        ion_networks = [
            network.Network(file_name) for file_name in input_file_names
        ]
        evidence_files = []
        for ion_network in ion_networks:
            local_file_name = os.path.basename(ion_network.file_name)
            if output_directory == "":
                output_path = os.path.dirname(ion_network.file_name)
            else:
                output_path = output_directory
            evidence_file_name = os.path.join(
                output_path,
                f"{local_file_name[:-9]}.evidence.hdf"
            )
            evidence_files.append(
                evidence.Evidence(
                    evidence_file_name=evidence_file_name,
                    ion_network=ion_network,
                    parameters=parameters,
                    logger=logger
                )
            )
        for index, evidence_file in enumerate(evidence_files[:-1]):
            for secondary_evidence_file in evidence_files[index + 1:]:
                evidence_file.mutual_collect_evidence_from(
                    secondary_evidence_file,
                    parameters=parameters,
                )


def show_ion_network(
    ion_network_file_name,
    evidence_file_name,
    parameter_file_name,
    log_file_name
):
    # TODO: Docstring
    # TODO: Implement
    if parameter_file_name is None:
        parameter_file_name = ""
    if log_file_name is None:
        log_file_name = ""
    parameters = utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
    )
    with utils.open_logger(log_file_name, parameters=parameters) as logger:
        inet = network.Network(ion_network_file_name)
        evi = evidence.Evidence(
            evidence_file_name=evidence_file_name,
            ion_network=inet,
            parameters=parameters,
            logger=logger
        )
        browser.Browser(
            inet,
            evi,
            logger
        )


def run_ion_network_gui():
    # TODO: Docstring
    with GUI() as gui:
        gui.run()


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
        self.window["Main"] = sg.Window("Main", self.window["Main"])
        self.active_window_name = "Main"
        if start:
            self.run()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        for window in list(self.window.values()):
            if not isinstance(window, list):
                window.close()

    def init_main_window(self):
        # TODO: Docstring
        self.window["Main"] = [
            [sg.Button("Convert", size=(self.widget_size, 1))],
            [sg.Button("Create", size=(self.widget_size, 1))],
            [sg.Button("Evidence", size=(self.widget_size, 1))],
            [sg.Button("Show", size=(self.widget_size, 1))],
        ]
        self.evaluate_window["Main"] = self.evaluate_main_window

    def init_convert_window(self):
        # TODO: Docstring
        self.window["Convert"] = [
            self.add_input_path_to_layout(
                file_types=(
                    ('HDMSE', '*.csv'),
                    ('SONAR', '*.csv'),
                    ('SWIMDIA', '*.csv'),
                    ('DDA', '*.mgf'),
                )
            ),
            self.add_output_directory_to_layout(),
            [
                sg.Text('Data type', size=(self.widget_size, 1)),
                sg.Combo(
                    ['DDA', 'HDMSE', "SONAR", "SWIMDIA"],
                    default_value='HDMSE',
                    key="data_type",
                    size=(self.widget_size * 2, 1)
                    # enable_events=True
                )
            ],
            self.add_parameter_file_to_layout(),
            self.add_log_file_to_layout(),
            self.add_main_menu_and_continue_buttons_to_layout()
        ]
        self.evaluate_window["Convert"] = self.evaluate_convert_window

    def init_create_window(self):
        # TODO: Docstring
        self.window["Create"] = [
            self.add_input_path_to_layout(
                file_types=(('Ion-networks', '*.inet.csv'),)
            ),
            self.add_output_directory_to_layout(),
            self.add_parameter_file_to_layout(),
            self.add_log_file_to_layout(),
            self.add_main_menu_and_continue_buttons_to_layout(),
        ]
        self.evaluate_window["Create"] = self.evaluate_create_window

    def init_evidence_window(self):
        # TODO: Docstring
        self.window["Evidence"] = [
            self.add_input_path_to_layout(
                file_types=(('Ion-networks', '*.inet.hdf'),)
            ),
            self.add_output_directory_to_layout(),
            self.add_parameter_file_to_layout(),
            self.add_log_file_to_layout(),
            self.add_main_menu_and_continue_buttons_to_layout(),
        ]
        self.evaluate_window["Evidence"] = self.evaluate_evidence_window

    def init_show_window(self):
        # TODO: Docstring
        # TODO: Implement
        self.window["Show"] = [
            self.add_input_path_to_layout(
                file_types=(('Ion-network', '*.inet.hdf'),),
                title="Ion-network",
                key="ion_network_file_name",
                multiple=False
            ),
            self.add_input_path_to_layout(
                file_types=(('Evidence', '*.evidence.hdf'),),
                title="Evidence",
                key="evidence_file_name",
                multiple=False
            ),
            self.add_main_menu_and_continue_buttons_to_layout(),
        ]
        self.evaluate_window["Show"] = self.evaluate_show_window

    def init_terminal_window(self):
        # TODO: Docstring
        self.window["Terminal"] = [
            [sg.Output(size=(150, 50))],
            self.add_main_menu_and_continue_buttons_to_layout(
                continue_button=False
            )
        ]

    def evaluate_main_window(self, event, values):
        # TODO: Docstring
        self.swap_active_window(event)

    def evaluate_convert_window(self, event, values):
        # TODO: Docstring
        if event == "Submit":
            self.run_terminal_command(
                convert_data_formats_to_csvs,
                values["input_path"].split(";"),
                values["data_type"],
                values["output_directory"],
                values["parameter_file_name"],
                values["log_file_name"]
            )

    def evaluate_create_window(self, event, values):
        # TODO: Docstring
        if event == "Submit":
            self.run_terminal_command(
                create_ion_networks,
                values["input_path"].split(";"),
                values["output_directory"],
                values["parameter_file_name"],
                values["log_file_name"]
            )

    def evaluate_evidence_window(self, event, values):
        # TODO: Docstring
        if event == "Submit":
            self.run_terminal_command(
                evidence_ion_networks,
                values["input_path"].split(";"),
                values["output_directory"],
                values["parameter_file_name"],
                values["log_file_name"]
            )

    def evaluate_show_window(self, event, values):
        # TODO: Docstring,
        # TODO: implement
        if event == "Submit":
            self.swap_active_window("")
            show_ion_network(
                values["ion_network_file_name"],
                values["evidence_file_name"],
                "",
                ""
            )
            self.swap_active_window("Main")

    def add_input_path_to_layout(
        self,
        file_types=(('ALL Files', '*.*'),),
        title="Input path",
        key="input_path",
        multiple=True
    ):
        # TODO: Docstring
        # TODO: Multiple and independent files?
        if multiple:
            browse_button = sg.FilesBrowse(
                size=(self.widget_size, 1),
                file_types=file_types
            )
        else:
            browse_button = sg.FileBrowse(
                size=(self.widget_size, 1),
                file_types=file_types
            )
        row = [
            sg.Text(title, size=(self.widget_size, 1)),
            sg.Input(
                key=key,
                size=(self.widget_size * 2, 1)
            ),
            browse_button,
        ]
        return row

    def add_output_directory_to_layout(self):
        # TODO: Docstring
        row = [
            sg.Text("Output directory", size=(self.widget_size, 1)),
            sg.Input(
                key="output_directory",
                size=(self.widget_size * 2, 1)
            ),
            sg.FolderBrowse(size=(self.widget_size, 1)),
        ]
        return row

    def add_parameter_file_to_layout(self):
        # TODO: Docstring
        row = [
            sg.Text("Parameter file", size=(self.widget_size, 1)),
            sg.Input(
                key="parameter_file_name",
                size=(self.widget_size * 2, 1),
            ),
            sg.FileBrowse(
                size=(self.widget_size, 1),
                file_types=(('Parameter', '*.json'),)
            ),
        ]
        return row

    def add_log_file_to_layout(self):
        # TODO: Docstring
        # TODO: remove overwrite warning
        row = [
            sg.Text("Log file", size=(self.widget_size, 1)),
            sg.Input(
                key="log_file_name",
                size=(self.widget_size * 2, 1),
            ),
            sg.FileSaveAs(size=(self.widget_size, 1)),
        ]
        return row

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

    def run_terminal_command(self, command, *args):
        # TODO: Docstring
        self.swap_active_window("Terminal")
        self.window["Terminal"].read(timeout=10)
        self.window["Terminal"]["Return to main menu"].Update(
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
                    button_color=("Black", "Red")
                )
                self.active_window_name = None
        self.window["Terminal"]["Return to main menu"].Update(
            text="Return to main menu",
            disabled=False
        )

    def swap_active_window(self, new_window_name=""):
        # TODO: Docstring, implement
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
