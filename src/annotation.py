#!python

# builtin
import logging
import os
# external
import h5py
import numpy as np
import numba
# local
import evidence as evi
import database as db


class Annotation(object):
    # TODO: Docstring

    def __init__(
        self,
        annotation_file_name=None,
        evidence=None,
        database=None,
        parameters={},
        logger=logging.getLogger()
    ):
        # TODO: Docstring
        self.logger = logger
        if evidence is not None:
            self.evidence = evidence
            self.file_name = os.path.join(
                evidence.file_name_path,
                evidence.file_name_base
            ) + ".annotation.hdf"
        elif annotation_file_name is not None:
            self.file_name = annotation_file_name
            self.evidence = evi.Evidence(
                os.path.join(
                    self.file_name_path,
                    self.file_name_base
                ) + ".evidence.hdf"
            )
        else:
            raise ValueError("Annotation file needs corresponding evidence.")
        self.ion_network = evidence.ion_network
        self.evidence.annotation = self
        self.ion_network.annotation = self
        if database is not None:
            if isinstance(database, str):
                database = db.Database(database)
            self.create_annotation(database, parameters)

    def create_annotation(self, database, parameters):
        self.write_parameters(database, parameters)
        self.write_candidates(database, parameters)
        # with h5py.File(self.file_name, "w") as hdf_file:
        #     self.write_parameters(hdf_file, database, parameters)
        #     self.write_candidates(hdf_file, database, parameters)

    def write_parameters(self, database, parameters):
        # TODO: Docstring
        with h5py.File(self.file_name, "w") as hdf_file:
            self.logger.info(f"Writing parameters to {self.file_name}")
            # TODO: Include prot, pep, frag counts
            hdf_file.attrs["database_file_name"] = str(database.file_name)
            for arg, value in parameters.items():
                if isinstance(value, str):
                    hdf_file.attrs[arg] = value
                else:
                    try:
                        iter(value)
                    except TypeError:
                        hdf_file.attrs[arg] = value
                    else:
                        hdf_file.attrs[arg] = str(value)

    def write_candidates(self, database, parameters):
        # TODO: Docstring
        self.logger.info(f"Writing candidates to {self.file_name}")
        low_peptide_indices, high_peptide_indices = self.get_candidate_peptide_indices_for_nodes(
            database,
            parameters
        )
        with h5py.File(self.file_name, "r+") as hdf_file:
            candidate_group = hdf_file.create_group("node_candidates")
            candidate_group.create_dataset(
                "low_peptide_indices",
                data=low_peptide_indices,
                compression="lzf",
            )
            candidate_group.create_dataset(
                "high_peptide_indices",
                data=high_peptide_indices,
                compression="lzf",
            )
        edge_indptr, edge_indices = self.get_candidate_peptide_indices_for_edges(
            database,
            low_peptide_indices,
            high_peptide_indices
        )
        with h5py.File(self.file_name, "r+") as hdf_file:
            candidate_group = hdf_file.create_group("edge_candidates")
            candidate_group.create_dataset(
                "indptr",
                data=edge_indptr,
                compression="lzf",
            )
            candidate_group.create_dataset(
                "indices",
                data=edge_indices,
                compression="lzf",
            )

    def get_candidate_peptide_indices_for_edges(
        self,
        database,
        low_peptide_indices,
        high_peptide_indices
    ):
        self.logger.info(f"Writing edge candidates to {self.file_name}")
        database_peptides = database.get_fragment_coordinates("peptide_index")
        indptr, indices = self.ion_network.get_edges(
            indptr_and_indices=True
        )
        max_batch = np.max(
            np.diff(indptr) * (high_peptide_indices - low_peptide_indices)
        )
        # TODO: remove init and just compile upfront?
        self.__get_candidate_peptide_indices_for_edges(
            indptr[:10],
            indices,
            low_peptide_indices,
            high_peptide_indices,
            database_peptides,
            max_batch
        )
        edge_indptr, edge_indices = self.__get_candidate_peptide_indices_for_edges(
            indptr,
            indices,
            low_peptide_indices,
            high_peptide_indices,
            database_peptides,
            max_batch
        )
        return edge_indptr, edge_indices

    @staticmethod
    @numba.njit
    def __get_candidate_peptide_indices_for_edges(
        indptr,
        indices,
        low_peptide_indices,
        high_peptide_indices,
        db_peptides,
        max_batch
    ):
        # TODO: Docstring
        result_indptr = np.empty(indptr[-1], np.int64)
        result_indices = np.empty(max_batch, np.int64)
        current_index = 0
        for start, end, low, high in zip(
            indptr[:-1],
            indptr[1:],
            low_peptide_indices,
            high_peptide_indices,
        ):
            if (low == high) or (start == end):
                result_indptr[start:end] = current_index
                continue
            if (
                (end - start) * (high - low) + current_index
            ) >= result_indices.shape[0]:
                new_result_indices = np.empty(
                    max_batch + result_indices.shape[0],
                    np.int64
                )
                new_result_indices[:result_indices.shape[0]] = result_indices
                result_indices = new_result_indices
            peptide_candidates = db_peptides[low: high]
            peptide_candidates_set = set(peptide_candidates)
            neighbors = indices[start: end]
            for i, neighbor in enumerate(neighbors):
                neighbor_low = low_peptide_indices[neighbor]
                neighbor_high = high_peptide_indices[neighbor]
                if neighbor_low == neighbor_high:
                    result_indptr[start + i] = current_index
                    continue
                neighbor_peptide_candidates = db_peptides[neighbor_low: neighbor_high]
                for neighbor_peptide_candidate in neighbor_peptide_candidates:
                    if neighbor_peptide_candidate in peptide_candidates_set:
                        result_indices[current_index] = neighbor_peptide_candidate
                        current_index += 1
                result_indptr[start + i] = current_index
        result_indptr[1:] = result_indptr[:-1]
        result_indptr[0] = 0
        return result_indptr, result_indices[:current_index]

    def get_candidate_peptide_indices_for_nodes(
        self,
        database,
        parameters
    ):
        # TODO: Docstring
        self.logger.info(f"Writing node candidates to {self.file_name}")
        max_ppm = parameters["annotation_ppm"]
        self_mzs = self.ion_network.get_ion_coordinates("FRAGMENT_MZ")
        mz_order = np.argsort(self_mzs)
        database_mzs = database.get_fragment_coordinates("mz")
        low_limits = np.searchsorted(
            np.log(database_mzs) * 10**6,
            np.log(self_mzs[mz_order]) * 10**6 - max_ppm,
            "left"
        )
        high_limits = np.searchsorted(
            np.log(database_mzs) * 10**6,
            np.log(self_mzs[mz_order]) * 10**6 + max_ppm,
            "right"
        )
        inv_order = np.argsort(mz_order)
        return low_limits[inv_order], high_limits[inv_order]

    @property
    def original_file_name(self):
        # TODO: Docstring
        with h5py.File(self.file_name, "r") as hdf_file:
            database_file_name = hdf_file.attrs["database_file_name"]
            database = db.Database(database_file_name)
        return database
