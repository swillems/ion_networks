#!python

import os
import sys
import logging

import h5py
import scipy
import pandas as pd
import numpy as np
import multiprocessing as mp


class Network(object):

    def __init__(
        self,
        raw_file_name=None,
        network_file_name=None,
        parameters=None
    ):
        # TODO
        print("create", raw_file_name, network_file_name, parameters)
        pass

    def align(
        self,
        other,
        alignment_file_name=None,
        parameters=None
    ):
        # TODO
        print("align", self, other, alignment_file_name, parameters)
        pass

    def evidence(
        self,
        alignment_file_name,
        evidence_file_name=None,
        parameters=None
    ):
        # TODO
        print("evidence", self, alignment_file_name, evidence_file_name, parameters)
        pass


if __name__ == "__main__":
    pass
