#!python

import json
import os


BASE_PATH = os.path.dirname(os.path.dirname(__file__))
LIB_PATH = os.path.join(BASE_PATH, "lib")
DEFAULT_PARAMETER_PATH = os.path.join(LIB_PATH, "default_parameters")
CREATE_PARAMETERS = "create_parameters.json"
ALIGN_PARAMETERS = "align_parameters.json"
EVIDENCE_PARAMETERS = "evidence_parameters.json"
SHOW_PARAMETERS = "show_parameters.json"


def read(default=None, file_name=None):
    """
    Read parameter file.

    Parameters
    ----------
    default : str
        The default parameters that should be loaded. Options are:
            "create"
            "align"
            "evidence"
            "show"
            None
    file_name : str
        The name of a .json file that contains parameters defined by the user.
        These will override the default parameters.

    Returns
    -------
    dict
        A dictionary with parameters.
    """
    if default is None:
        parameters = {}
    else:
        if default == "create":
            default_parameter_file_name = os.path.join(
                DEFAULT_PARAMETER_PATH,
                CREATE_PARAMETERS
            )
        elif default == "align":
            default_parameter_file_name = os.path.join(
                DEFAULT_PARAMETER_PATH,
                ALIGN_PARAMETERS
            )
        elif default == "evidence":
            default_parameter_file_name = os.path.join(
                DEFAULT_PARAMETER_PATH,
                EVIDENCE_PARAMETERS
            )
        elif default == "show":
            default_parameter_file_name = os.path.join(
                DEFAULT_PARAMETER_PATH,
                SHOW_PARAMETERS
            )
        with open(default_parameter_file_name, "r") as in_file:
            parameters = json.load(in_file)
    if file_name is not None:
        with open(file_name, "r") as in_file:
            user_defined_parameters = json.load(in_file)
        parameters.update(user_defined_parameters)
    return parameters
