#!python

# builtin
import os
import threading
# external
import PySimpleGUI as sg
import click
# local
try:
    from . import ms_run_files
    from . import ms_database
    from . import ms_utils
    from . import browser
except (ImportError, ModuleNotFoundError):
    import ms_run_files
    import ms_database
    import ms_utils
    import browser


def convert_data_formats_to_csvs(
    input_path,
    data_type,
    output_directory,
    parameter_file_name,
    log_file_name,
    threads=None
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
            'DDA'
            'SONAR'
            'HDMSE'
            'SWIMDIA'
            'DIAPASEF'
    parameter_file_name : str or None
        If provided, parameters will be read from this file.
    log_file_name : str or None
        If provided, all logs will be written to this file.
    """
    if parameter_file_name is None:
        parameter_file_name = ""
    if parameter_file_name != "":
        parameter_file_name = os.path.abspath(parameter_file_name)
    if output_directory is None:
        output_directory = ""
    if output_directory != "":
        output_directory = os.path.abspath(output_directory)
    parameters = ms_utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
        default="convert"
    )
    if threads is not None:
        ms_utils.set_threads(threads)
    # TODO: Proper parsing of empty log...?
    if (log_file_name is None) or (log_file_name == ""):
        log_file_name = parameters["log_file_name"]
        if log_file_name == "":
            log_file_name = output_directory
    with ms_utils.open_logger(log_file_name) as logger:
        logger.info(f"Converting to generic input csvs")
        input_file_names = ms_utils.get_file_names_with_extension(
            input_path,
            extension=ms_utils.DATA_TYPE_FILE_EXTENSIONS[data_type]
        )
        file_count = len(input_file_names)
        logger.info(
            f"{file_count} input_file{'s' if file_count != 1 else ''}"
            f": {input_file_names}"
        )
        logger.info(f"data_type: {data_type}")
        logger.info(f"output_directory: {output_directory}")
        logger.info(f"parameter_file_name: {parameter_file_name}")
        logger.info(f"max_threads: {ms_utils.MAX_THREADS}")
        logger.info("")
        if output_directory != "":
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)
        for input_file_name in sorted(input_file_names):
            data = ms_utils.read_data_from_file(
                data_type,
                input_file_name,
            )
            file_name_base = os.path.basename(input_file_name)[
                :-len(ms_utils.DATA_TYPE_FILE_EXTENSIONS[data_type])
            ]
            if output_directory == "":
                output_path = os.path.dirname(input_file_name)
            else:
                output_path = output_directory
            output_file_name = os.path.join(
                output_path,
                f"{file_name_base}.inet.csv"
            )
            ms_utils.write_data_to_csv_file(data, output_file_name)


def create_ion_networks(
    input_path,
    output_directory,
    parameter_file_name,
    log_file_name,
    threads=None
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
    if parameter_file_name != "":
        parameter_file_name = os.path.abspath(parameter_file_name)
    if output_directory is None:
        output_directory = ""
    if output_directory != "":
        output_directory = os.path.abspath(output_directory)
    parameters = ms_utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
        default="create"
    )
    if threads is not None:
        ms_utils.set_threads(threads)
    # TODO: Proper parsing of empty log...?
    if (log_file_name is None) or (log_file_name == ""):
        log_file_name = parameters["log_file_name"]
        if log_file_name == "":
            log_file_name = output_directory
    with ms_utils.open_logger(log_file_name) as logger:
        logger.info(f"Creating ion-networks")
        input_file_names = ms_utils.get_file_names_with_extension(
            input_path,
            ".inet.csv"
        )
        file_count = len(input_file_names)
        logger.info(
            f"{file_count} input_file{'s' if file_count != 1 else ''}"
            f": {input_file_names}"
        )
        logger.info(f"output_directory: {output_directory}")
        logger.info(f"parameter_file_name: {parameter_file_name}")
        logger.info(f"max_threads: {ms_utils.MAX_THREADS}")
        logger.info("")
        for centroids_file_name in input_file_names:
            local_file_name = os.path.basename(centroids_file_name)
            if output_directory == "":
                output_path = os.path.dirname(centroids_file_name)
            else:
                output_path = output_directory
            ion_network_file_name = os.path.join(
                output_path,
                f"{local_file_name[:-9]}.inet.hdf"
            )
            network = ms_run_files.HDF_Network_File(
                ion_network_file_name,
                new_file=True
            )
            network.create(
                centroids_file_name=centroids_file_name,
                parameters=parameters
            )


def evidence_ion_networks(
    input_path,
    parameter_file_name,
    log_file_name,
    threads=None
):
    """
    Evidence ion-networks with each other.

    Parameters
    ----------
    input_path : iterable[str]
        An iterable with file and/or folder names.
    parameter_file_name : str or None
        If provided, parameters will be read from this file.
    log_file_name : str or None
        If provided, all logs will be written to this file.
    """
    if parameter_file_name is None:
        parameter_file_name = ""
    if parameter_file_name != "":
        parameter_file_name = os.path.abspath(parameter_file_name)
    parameters = ms_utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
        default="evidence"
    )
    if threads is not None:
        ms_utils.set_threads(threads)
    # TODO: Proper parsing of empty log...?
    if (log_file_name is None) or (log_file_name == ""):
        log_file_name = parameters["log_file_name"]
    with ms_utils.open_logger(log_file_name) as logger:
        logger.info(f"Evidencing ion-networks")
        input_file_names = ms_utils.get_file_names_with_extension(
            input_path,
            ".inet.hdf"
        )
        file_count = len(input_file_names)
        logger.info(
            f"{file_count} input_file{'s' if file_count != 1 else ''}"
            f": {input_file_names}"
        )
        logger.info(f"parameter_file_name: {parameter_file_name}")
        logger.info(f"max_threads: {ms_utils.MAX_THREADS}")
        logger.info("")
        evidence_files = [
            ms_run_files.HDF_Evidence_File(
                file_name,
                new_file=parameters["force_overwrite"],
                is_read_only=False,
            ) for file_name in input_file_names
        ]
        for index, evidence_file in enumerate(evidence_files[:-1]):
            logger.info(f"Caching edges of {evidence_file.file_name}")
            indptr, indices, pointers = evidence_file.ion_network.get_edges(
                symmetric=True,
                return_pointers=True
            )
            for secondary_evidence_file in evidence_files[index + 1:]:
                evidence_file.align(
                    secondary_evidence_file,
                    parameters=parameters,
                    indptr=indptr,
                    indices=indices,
                    pointers=pointers,
                )
            del indptr, indices, pointers


def show_ion_network(
    parameter_file_name,
    log_file_name
):
    # TODO: Docstring
    # TODO: Implementation updates
    if parameter_file_name is None:
        parameter_file_name = ""
    if parameter_file_name != "":
        parameter_file_name = os.path.abspath(parameter_file_name)
    parameters = ms_utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
    )
    # TODO: Proper parsing of empty log...?
    if (log_file_name is None) or (log_file_name == ""):
        log_file_name = parameters["log_file_name"]
    with ms_utils.open_logger(log_file_name) as logger:
        logger.info(f"Showing ion-networks")
        logger.info("")
        with browser.Browser() as browser_object:
            browser_object.run()


def run_ion_network_gui():
    # TODO: Docstring
    with GUI() as gui:
        gui.run()


def create_database(
    input_path,
    output_directory,
    parameter_file_name,
    log_file_name,
    threads=None
):
    # TODO: Docstring
    if parameter_file_name is None:
        parameter_file_name = ""
    if parameter_file_name != "":
        parameter_file_name = os.path.abspath(parameter_file_name)
    if output_directory is None:
        output_directory = ""
    if output_directory != "":
        output_directory = os.path.abspath(output_directory)
    parameters = ms_utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
        default="database"
    )
    if threads is not None:
        ms_utils.set_threads(threads)
    # TODO: Proper parsing of empty log...?
    if (log_file_name is None) or (log_file_name == ""):
        log_file_name = parameters["log_file_name"]
        if log_file_name == "":
            log_file_name = output_directory
    # TODO: turn off ms2pip logger?
    with ms_utils.open_logger(log_file_name) as logger:
        logger.info(f"Creating database")
        input_file_names = ms_utils.get_file_names_with_extension(
            input_path,
            extension=".fasta"
        )
        file_count = len(input_file_names)
        logger.info(
            f"{file_count} input_file{'s' if file_count != 1 else ''}"
            f": {input_file_names}"
        )
        logger.info(f"output_directory: {output_directory}")
        # TODO: Refer to default parameter file if necessary
        logger.info(f"parameter_file_name: {parameter_file_name}")
        logger.info(f"max_threads: {ms_utils.MAX_THREADS}")
        # TODO: include relevant individual parameters?
        logger.info("")
        base_name = "_".join(
            [
                ".".join(
                    os.path.basename(fasta_file_name).split(".")[:-1]
                ) for fasta_file_name in input_file_names
            ]
        )
        database_file_name = os.path.join(output_directory, base_name)
        if parameters["create_targets"]:
            if parameters["create_decoys"]:
                database_file_name = f"{database_file_name}_concatenated_decoy.hdf"
            else:
                database_file_name = f"{database_file_name}.hdf"
        else:
            database_file_name = f"{database_file_name}_decoy.hdf"
        db = ms_database.HDF_Database_File(
            database_file_name,
            new_file=True,
        )
        db.create_from_fastas(input_file_names, parameters)


def annotate(
    input_path,
    database_file_name,
    mgf_format,
    parameter_file_name,
    log_file_name,
    threads=None,
    fragment_ppm=None,
    export_decoys=None,
    fdr_filter=None
):
    # TODO: Docstring
    if parameter_file_name is None:
        parameter_file_name = ""
    if parameter_file_name != "":
        parameter_file_name = os.path.abspath(parameter_file_name)
    parameters = ms_utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
        default="annotation"
    )
    if fragment_ppm is not None:
        parameters["annotation_ppm"] = fragment_ppm
    if export_decoys is not None:
        parameters["export_decoys"] = export_decoys
    if fdr_filter is not None:
        parameters["fdr_filter"] = fdr_filter
    if threads is not None:
        ms_utils.set_threads(threads)
    # TODO: Proper parsing of empty log...?
    if (log_file_name is None) or (log_file_name == ""):
        log_file_name = parameters["log_file_name"]
    with ms_utils.open_logger(log_file_name) as logger:
        if mgf_format:
            logger.info(f"Annotating mgf files")
            input_file_names = ms_utils.get_file_names_with_extension(
                input_path,
                extension=".mgf"
            )
        else:
            logger.info(f"Annotating ion-networks")
            input_file_names = ms_utils.get_file_names_with_extension(
                input_path,
                extension=".evidence.hdf"
            )
        file_count = len(input_file_names)
        logger.info(
            f"{file_count} input_file{'s' if file_count != 1 else ''}"
            f": {input_file_names}"
        )
        logger.info(f"database_file_name: {database_file_name}")
        logger.info(f"parameter_file_name: {parameter_file_name}")
        logger.info(f"max_threads: {ms_utils.MAX_THREADS}")
        logger.info(f"fragment_ppm: {parameters['annotation_ppm']}")
        logger.info(f"fdr_filter: {parameters['fdr_filter']}")
        # TODO implement fdr filtering!!!
        logger.info(f"export_decoys: {parameters['export_decoys']}")
        logger.info("")
        database = ms_database.HDF_Database_File(database_file_name)
        for file_name in input_file_names:
            if mgf_format:
                out_file_name = f"{file_name[:-4]}.annotation.csv"
                ms_utils.annotate_mgf(
                    file_name,
                    database,
                    out_file_name,
                    parameters,
                )
            else:
                out_file_name = f"{'.'.join(file_name.split('.')[:-2])}.csv"
                evidence = ms_run_files.HDF_Evidence_File(file_name)
                evidence.annotate(database, out_file_name, parameters)


def create_mgfs(
    input_path,
    output_directory,
    parameter_file_name,
    log_file_name,
    threads=None
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
    if parameter_file_name != "":
        parameter_file_name = os.path.abspath(parameter_file_name)
    if output_directory is None:
        output_directory = ""
    if output_directory != "":
        output_directory = os.path.abspath(output_directory)
    parameters = ms_utils.read_parameters_from_json_file(
        file_name=parameter_file_name,
        default="mgf"
    )
    if threads is not None:
        ms_utils.set_threads(threads)
    # TODO: Proper parsing of empty log...?
    if (log_file_name is None) or (log_file_name == ""):
        log_file_name = parameters["log_file_name"]
        if log_file_name == "":
            log_file_name = output_directory
    with ms_utils.open_logger(log_file_name) as logger:
        logger.info(f"Creating ion-networks")
        input_file_names = ms_utils.get_file_names_with_extension(
            input_path,
            ".evidence.hdf"
        )
        file_count = len(input_file_names)
        logger.info(
            f"{file_count} input_file{'s' if file_count != 1 else ''}"
            f": {input_file_names}"
        )
        logger.info(f"output_directory: {output_directory}")
        logger.info(f"parameter_file_name: {parameter_file_name}")
        logger.info(f"max_threads: {ms_utils.MAX_THREADS}")
        logger.info("")
        for input_file in input_file_names:
            evidence = ms_run_files.HDF_Evidence_File(input_file)
            evidence.create_mgf(parameters=parameters)


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
        update_message = ms_utils.verify_version()
        if update_message != "":
            sg.popup(
                update_message,
                title="Update available",
                non_blocking=True
            )
        self.init_main_window()
        self.init_convert_window()
        self.init_create_window()
        self.init_evidence_window()
        self.init_terminal_window()
        self.window["Main"] = sg.Window(
            "Main",
            self.window["Main"],
            finalize=True
        )
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
        default_data_type = "HDMSE"
        self.window["Convert"] = [
            [
                sg.Text('Data type', size=(self.widget_size, 1)),
                sg.Combo(
                    sorted(ms_utils.DATA_TYPE_FILE_EXTENSIONS),
                    default_value=default_data_type,
                    key="data_type",
                    size=(self.widget_size * 2, 1)
                    # enable_events=True
                )
            ],
            self.add_input_path_to_layout(
                file_types=(
                    (
                        key,
                        f"*{value}"
                    ) for key, value in sorted(
                        ms_utils.DATA_TYPE_FILE_EXTENSIONS.items()
                    )
                ),
                # TODO: default_value=default_data_type,
            ),
            self.add_output_directory_to_layout(),
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
            self.add_parameter_file_to_layout(),
            self.add_log_file_to_layout(),
            self.add_main_menu_and_continue_buttons_to_layout(),
        ]
        self.evaluate_window["Evidence"] = self.evaluate_evidence_window

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
        if event == "Show":
            self.swap_active_window("")
            show_ion_network(
                "",
                ""
            )
            self.swap_active_window("Main")
        else:
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
                values["parameter_file_name"],
                values["log_file_name"]
            )

    def add_input_path_to_layout(
        self,
        file_types=(('ALL Files', '*.*'),),
        title="Input path",
        key="input_path",
        multiple=True,
        default_value=None,
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
                file_types=file_types,
                # TODO: default file_type?
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

    def add_parameter_file_to_layout(self, default=""):
        # TODO: Docstring
        row = [
            sg.Text("Parameter file", size=(self.widget_size, 1)),
            sg.Input(
                key="parameter_file_name",
                size=(self.widget_size * 2, 1),
                default_text=default
            ),
            sg.FileBrowse(
                size=(self.widget_size, 1),
                file_types=(('Parameter', '*.json'),)
            ),
        ]
        return row

    def add_log_file_to_layout(self, default=""):
        # TODO: Docstring
        # TODO: remove overwrite warning
        # TODO: default log / empty log not parsed properly
        row = [
            sg.Text("Log file", size=(self.widget_size, 1)),
            sg.Input(
                key="log_file_name",
                size=(self.widget_size * 2, 1),
                default_text=default
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
                    self.window[new_window_name],
                    finalize=True
                )
            self.window[new_window_name].UnHide()
        self.active_window_name = new_window_name


class CLI(object):

    CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

    def __init__(self):
        print(ms_utils.verify_version())
        # TODO: process update message
        self.main.add_command(CLI.convert)
        self.main.add_command(CLI.create)
        self.main.add_command(CLI.evidence)
        self.main.add_command(CLI.show)
        self.main.add_command(CLI.gui)
        self.main.add_command(CLI.database)
        self.main.add_command(CLI.annotate)
        self.main.add_command(CLI.mgf)
        self.main()

    @staticmethod
    @click.group(
        context_settings=CONTEXT_SETTINGS,
        help="Analysis of LC-[...]-MSMS data with ion-networks."
    )
    def main():
        pass

    @staticmethod
    @click.command(
        "convert",
        help="Convert [input.*] files with centroided ions to unified "
            "[input.inet.csv] csv files.",
        short_help="Convert various input formats to unified input."
    )
    @click.option(
        "--input_path",
        "-i",
        help="An [input.*] file with centroided ion peaks that needs to be "
            "converted to a unified [input.inet.csv] file. "
            "Individual files can be provided, as well as folders. "
            "This flag can be set multiple times.",
        multiple=True,
        required=True,
        type=click.Path(exists=True)
    )
    @click.option(
        '--data_type',
        '-d',
        help="The data type of the [input.*] file. If this is DDA, a [*.mgf] "
        "file that was centroided with ms-convert is expected with the field "
        "RTINSECONDS as LC coordinate. "
        "For HDMSE, SONAR and SWIM-DIA, a [*_Apex3DIons.csv] generated with "
        "Waters' Apex3d is expected, typically generated as follows "
        "'Apex3D64.exe -pRawDirName sample.raw -outputDirName "
        "peak_picked_sample_folder -lockMassZ2 785.8426 "
        "-lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 "
        "-leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 "
        "-bEnableCentroids 0'. "
        "For DIAPASEF, a [*_centroided.hdf] file generated with diapasef.py "
        "(https://github.com/swillems/diapasef) is expected.",
        required=True,
        type=click.Choice(
            ['DDA', 'HDMSE', "SONAR", "SWIMDIA", "DIAPASEF"],
            case_sensitive=True
        )
    )
    @click.option(
        "--output_directory",
        "-o",
        help="For each [input.*] file, an [input.inet.csv] file is created. "
            "If no output directory is provided, each [input.inet.csv] file is "
            "placed in the same folder as its corresponding [input.*] file. "
            "This output directory can also be supplied through a "
            "[parameters.json] file. "
            "WARNING: This overrides already existing files without "
            "confirmation.",
        type=click.Path(file_okay=False),
    )
    @click.option(
        "--parameter_file",
        "-p",
        "parameter_file_name",
        help="A [parameters.json] file with optional parameters.",
        type=click.Path(exists=True, dir_okay=False),
    )
    @click.option(
        "--log_file",
        "-l",
        "log_file_name",
        help="Save the log to a [log.txt] file. "
            "By default this is written to the current directory. "
            "This log file can also be supplied through a [parameters.json] "
            "file. It can be turned off by providing an empty path (i.e. ''). "
            "If the log file already exists, the new log data is appended.",
        type=click.Path(),
    )
    def convert(
        input_path,
        output_directory,
        data_type,
        parameter_file_name,
        log_file_name
    ):
        convert_data_formats_to_csvs(
            input_path,
            data_type,
            output_directory,
            parameter_file_name,
            log_file_name
        )

    @staticmethod
    @click.command(
        "create",
        help="Create [input.inet.hdf] ion-network files from unified "
            "[input.inet.csv] files.",
        short_help="Create ion-networks from unified input."
    )
    @click.option(
        "--input_path",
        "-i",
        help="A unified [input.inet.csv] file with centroided ion peaks. "
            "Columns with headers PRECURSOR_RT, FRAGMENT_MZ and "
            "FRAGMENT_LOGINT always need to be present. "
            "All columns whose header start with # are not interpreted as ion "
            "coordinates but as comments such as e.g. prior annotations. "
            "All other columns (e.g. PRECURSOR_MZ or PRECURSOR_DT) are "
            "automatically interpreted as dimensions with ion coordinates. "
            "All PRECURSOR_* dimensions are used to create edges between ions. "
            "Individual files can be provided, as well as folders. "
            "This flag can be set multiple times.",
        multiple=True,
        required=True,
        type=click.Path(exists=True)
    )
    @click.option(
        "--output_directory",
        "-o",
        help="For each [input.inet.csv] file, an [input.inet.hdf] ion-network "
            "file is created. "
            "If no output directory is provided, each [input.inet.hdf] file is "
            "placed in the same folder as its corresponding [input.inet.csv] "
            "file. "
            "This output directory can also be supplied through a "
            "[parameters.json] file. "
            "WARNING: This overrides already existing files without "
            "confirmation.",
        type=click.Path(file_okay=False)
    )
    @click.option(
        "--parameter_file",
        "-p",
        "parameter_file_name",
        help="A [parameters.json] file with optional parameters.",
        type=click.Path(exists=True, dir_okay=False)
    )
    @click.option(
        "--log_file",
        "-l",
        "log_file_name",
        help="Save the log to a [log.txt] file. "
            "By default this is written to the current directory. "
            "This log file can also be supplied through a [parameters.json] "
            "file. It can be turned off by providing an empty path (i.e. ''). "
            "If the log file already exists, the new log data is appended.",
        type=click.Path()
    )
    @click.option(
        "--threads",
        "-t",
        "threads",
        help="The number of threads to use.",
        type=int
    )
    def create(
        input_path,
        output_directory,
        parameter_file_name,
        log_file_name,
        threads
    ):
        create_ion_networks(
            input_path,
            output_directory,
            parameter_file_name,
            log_file_name,
            threads
        )

    @staticmethod
    @click.command(
        "evidence",
        help="Collect pairwise evidence for [input.inet.hdf] ion-network files "
            "as [input.evidence.hdf] evidence files.",
        short_help="Collect evidence for ion-networks."
    )
    @click.option(
        "--input_path",
        "-i",
        help="An [input.inet.hdf] ion-network file."
            "Individual files can be provided, as well as folders. "
            "This flag can be set multiple times.",
        required=True,
        multiple=True,
        type=click.Path(exists=True)
    )
    @click.option(
        "--parameter_file",
        "-p",
        "parameter_file_name",
        help="A [parameters.json] file with optional parameters.",
        type=click.Path(exists=True, dir_okay=False)
    )
    @click.option(
        "--log_file",
        "-l",
        "log_file_name",
        help="Save the log to a [log.txt] file. "
            "By default this is written to the current directory. "
            "This log file can also be supplied through a [parameters.json] "
            "file. It can be turned off by providing an empty path (i.e. ''). "
            "If the log file already exists, the new log data is appended.",
        type=click.Path()
    )
    @click.option(
        "--threads",
        "-t",
        "threads",
        help="The number of threads to use.",
        type=int
    )
    def evidence(
        input_path,
        parameter_file_name,
        log_file_name,
        threads
    ):
        evidence_ion_networks(
            input_path,
            parameter_file_name,
            log_file_name,
            threads
        )

    @staticmethod
    # TODO: Implement
    @click.command(
        "show",
        help="Show and browse ion-networks and their evidence.",
        short_help="Show and browse ion-networks."
    )
    @click.option(
        "--parameter_file",
        "-p",
        "parameter_file_name",
        help="A [parameters.json] file with optional parameters.",
        type=click.Path(exists=True, dir_okay=False)
    )
    @click.option(
        "--log_file",
        "-l",
        "log_file_name",
        help="Save the log to a [log.txt] file. "
            "By default this is written to the current directory. "
            "This log file can also be supplied through a [parameters.json] "
            "file. It can be turned off by providing an empty path (i.e. ''). "
            "If the log file already exists, the new log data is appended.",
        type=click.Path()
    )
    def show(
        parameter_file_name,
        log_file_name
    ):
        show_ion_network(
            parameter_file_name,
            log_file_name
        )

    @staticmethod
    @click.command(
        "gui",
        help="Graphical user interface to analyse ion-networks.",
    )
    def gui():
        run_ion_network_gui()

    @staticmethod
    @click.command(
        "database",
        help="Create a [database.hdf] from fasta files.",
        short_help="Create database from fasta files."
    )
    @click.option(
        "--input_path",
        "-i",
        help="A fasta file with protein sequences. "
            "Individual files can be provided, as well as folders. "
            "This flag can be set multiple times.",
        multiple=True,
        required=True,
        type=click.Path(exists=True)
    )
    @click.option(
        "--output_directory",
        "-o",
        help="The output directory fot the database. The file name is "
            "automatically set as a concatenation of the input fasta files, "
            "potentially appedned with _concatenated_decoy. "
            "This output directory can also be supplied through a "
            "[parameters.json] file. "
            "WARNING: This overrides already existing files without "
            "confirmation.",
        type=click.Path(file_okay=False)
    )
    @click.option(
        "--parameter_file",
        "-p",
        "parameter_file_name",
        help="A [parameters.json] file with optional parameters.",
        type=click.Path(exists=True, dir_okay=False)
    )
    @click.option(
        "--log_file",
        "-l",
        "log_file_name",
        help="Save the log to a [log.txt] file. "
            "By default this is written to the current directory. "
            "This log file can also be supplied through a [parameters.json] "
            "file. It can be turned off by providing an empty path (i.e. ''). "
            "If the log file already exists, the new log data is appended.",
        type=click.Path()
    )
    @click.option(
        "--threads",
        "-t",
        "threads",
        help="The number of threads to use.",
        type=int
    )
    def database(
        input_path,
        output_directory,
        parameter_file_name,
        log_file_name,
        threads
    ):
        create_database(
            input_path,
            output_directory,
            parameter_file_name,
            log_file_name,
            threads
        )

    @staticmethod
    @click.command(
        "annotate",
        help="Annotate ion-network files.",
        short_help="Annotate ion-network files."
    )
    @click.option(
        "--input_path",
        "-i",
        help="An [input.evidence.hdf] evidence file. "
            "Individual files can be provided, as well as folders. "
            "This flag can be set multiple times.",
        required=True,
        multiple=True,
        type=click.Path(exists=True)
    )
    @click.option(
        "--database_file",
        "-d",
        "database_file_name",
        help="A [database.hdf] file.",
        required=True,
        type=click.Path(exists=True, dir_okay=False),
    )
    @click.option(
        '--mgf_format',
        '-m',
        'mgf_format',
        help="Input is in mgf format",
        is_flag=True,
        default=False,
        show_default=True,
    )
    @click.option(
        "--parameter_file",
        "-p",
        "parameter_file_name",
        help="A [parameters.json] file with optional parameters.",
        type=click.Path(exists=True, dir_okay=False)
    )
    @click.option(
        "--log_file",
        "-l",
        "log_file_name",
        help="Save the log to a [log.txt] file. "
            "By default this is written to the current directory. "
            "This log file can also be supplied through a [parameters.json] "
            "file. It can be turned off by providing an empty path (i.e. ''). "
            "If the log file already exists, the new log data is appended.",
        type=click.Path()
    )
    @click.option(
        "--threads",
        "-t",
        "threads",
        help="The number of threads to use.",
        type=int
    )
    @click.option(
        "--fragment_ppm_error",
        "-e",
        "fragment_ppm_error",
        help="The maximum allowed fragment ppm error.",
        type=float
    )
    @click.option(
        "--export_decoys",
        "-x",
        "export_decoys",
        help="Include decoy hits in the export.",
        is_flag=True,
        default=False,
        show_default=True,
    )
    @click.option(
        "--fdr_filter",
        "-f",
        "fdr_filter",
        help="Remove annotations above this fdr.",
        type=float
    )
    def annotate(
        input_path,
        database_file_name,
        mgf_format,
        parameter_file_name,
        log_file_name,
        threads,
        fragment_ppm_error,
        export_decoys,
        fdr_filter
    ):
        annotate(
            input_path,
            database_file_name,
            mgf_format,
            parameter_file_name,
            log_file_name,
            threads,
            fragment_ppm_error,
            export_decoys,
            fdr_filter
        )

    @staticmethod
    @click.command(
        "mgf",
        help="Create [input.MGF] ion-network files from evidence "
            "[input.evidence.hdf] files.",
        short_help="Create mgf files from evidence."
    )
    @click.option(
        "--input_path",
        "-i",
        help="An evidence [input.evidence.hdf] file. "
            "Individual files can be provided, as well as folders. "
            "This flag can be set multiple times.",
        multiple=True,
        required=True,
        type=click.Path(exists=True)
    )
    @click.option(
        "--output_directory",
        "-o",
        help="For each [input.evidence.hdf] file, an [input.mgf] "
            "file is created. "
            "If no output directory is provided, each [input.mgf] file is "
            "placed in the same folder as its corresponding "
            "[input.evidence.hdf] file. "
            "This output directory can also be supplied through a "
            "[parameters.json] file. "
            "WARNING: This overrides already existing files without "
            "confirmation.",
        type=click.Path(file_okay=False)
    )
    @click.option(
        "--parameter_file",
        "-p",
        "parameter_file_name",
        help="A [parameters.json] file with optional parameters.",
        type=click.Path(exists=True, dir_okay=False)
    )
    @click.option(
        "--log_file",
        "-l",
        "log_file_name",
        help="Save the log to a [log.txt] file. "
            "By default this is written to the current directory. "
            "This log file can also be supplied through a [parameters.json] "
            "file. It can be turned off by providing an empty path (i.e. ''). "
            "If the log file already exists, the new log data is appended.",
        type=click.Path()
    )
    def mgf(
        input_path,
        output_directory,
        parameter_file_name,
        log_file_name
    ):
        create_mgfs(
            input_path,
            output_directory,
            parameter_file_name,
            log_file_name
        )


# TODO: Rename "unified" and "peaks" referring to .inet.csv files?
# TODO: Define help text in separate json files?
# TODO: Show help text popups in GUI
# TODO: Database interface (CLI + GUI)
# TODO: Annotation interface (CLI + GUI)


# -pRawDirName ~/sandbox/HDMSE_test/171114_HDMSE_Mclass_K562_30min_01.raw -outputDirName ~/sandbox/HDMSE_test/ -lockMassZ2 785.8426 -lockmassToleranceAMU 0.25 -bCSVOutput 1 -writeFuncCsvFiles 0 -leThresholdCounts 1 -heThresholdCounts 1 -apexTrackSNRThreshold 1 -bEnableCentroids 0
