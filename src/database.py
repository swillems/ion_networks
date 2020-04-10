#!python

# builtin
import logging
# external
import pandas as pd
import numpy as np
import ms2pip.ms2pipC
import ms2pip.retention_time
import pyteomics.parser
import pyteomics.fasta
import h5py
# local
import utils


class Database(object):
    # TODO: Docstring

    def __init__(
        self,
        database_file_name,
        fasta_file_names=None,
        parameters={},
        logger=logging.getLogger(),
    ):
        # TODO: Docstring
        self.logger = logger
        self.file_name = database_file_name
        if fasta_file_names is not None:
            self.create_from_fastas(fasta_file_names, parameters)

    def create_from_fastas(
        self,
        fasta_file_names,
        parameters,
    ):
        # TODO: Docstring
        self.logger.info(f"Creating database {self.file_name}")
        proteins, peptides = self.read_proteins_and_peptides_from_fasta(
            fasta_file_names,
            **parameters
        )
        with h5py.File(self.file_name, "w") as hdf_file:
            self.write_proteins(proteins, hdf_file, **parameters)
            peprec = self.write_peptides(hdf_file, peptides, **parameters)
            self.write_fragments(hdf_file, peprec, **parameters)
            self.write_parameters(hdf_file, fasta_file_names, parameters)

    def write_parameters(self, hdf_file, fasta_file_names, parameters):
        # TODO: Docstring
        self.logger.info(f"Writing parameters to {self.file_name}")
        # TODO: Include prot, pep, frag counts
        hdf_file.attrs["fasta_file_names"] = str(fasta_file_names)
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

    def write_dataframe_to_hdf_group_name(
        self,
        hdf_file,
        group_name,
        data,
        **kwargs,
    ):
        # TODO: Docstring
        if group_name in hdf_file:
            del hdf_file[group_name]
        group = hdf_file.create_group(group_name)
        for column in data.columns:
            try:
                group.create_dataset(
                    column,
                    data=data[column],
                    compression="lzf",
                )
            except TypeError:
                group.create_dataset(
                    column,
                    data=data[column],
                    compression="lzf",
                    dtype=h5py.string_dtype()
                )

    def read_proteins_and_peptides_from_fasta(
        self,
        fasta_file_names,
        protease,
        missed_cleavages,
        min_peptide_length,
        max_peptide_length,
        standard_amino_acids,
        create_targets,
        create_decoys,
        decoy_prefix,
        decoy_suffix,
        **kwargs,
    ):
        # TODO: Docstring
        proteins = {}
        peptides = {}
        if not (create_targets or create_decoys):
            raise ValueError("No targets or decoys to create")
        for fasta_file_name in fasta_file_names:
            self.logger.info(f"Reading {fasta_file_name}")
            if create_targets:
                reversed_protein_decoy = False
                proteins, peptides = self.__read_proteins_and_peptides_from_fasta(
                    fasta_file_name,
                    protease,
                    missed_cleavages,
                    min_peptide_length,
                    max_peptide_length,
                    standard_amino_acids,
                    reversed_protein_decoy,
                    proteins,
                    peptides,
                    decoy_prefix,
                    decoy_suffix,
                    **kwargs,
                )
            if create_decoys:
                reversed_protein_decoy = True
                proteins, peptides = self.__read_proteins_and_peptides_from_fasta(
                    fasta_file_name,
                    protease,
                    missed_cleavages,
                    min_peptide_length,
                    max_peptide_length,
                    standard_amino_acids,
                    reversed_protein_decoy,
                    proteins,
                    peptides,
                    decoy_prefix,
                    decoy_suffix,
                    **kwargs,
                )
        return proteins, peptides

    def __read_proteins_and_peptides_from_fasta(
        self,
        fasta_file_name,
        protease,
        missed_cleavages,
        min_peptide_length,
        max_peptide_length,
        standard_amino_acids,
        reversed_protein_decoy,
        proteins,
        peptides,
        decoy_prefix,
        decoy_suffix,
        **kwargs,
    ):
        # TODO: Docstring
        protein_index = 0
        for description, sequence in pyteomics.fasta.FASTA(fasta_file_name):
            protein_info = pyteomics.fasta.parse(description)
            protein = protein_info["entry"]
            del protein_info["entry"]
            if reversed_protein_decoy:
                protein = f"{decoy_prefix}{protein}{decoy_suffix}"
                sequence = sequence[::-1]
            protein_info["sequence"] = sequence
            protein_info["index"] = protein_index
            if protein in proteins:
                continue
            proteins[protein] = protein_info
            for peptide in pyteomics.parser.cleave(
                sequence,
                pyteomics.parser.expasy_rules[protease],
                missed_cleavages
            ):
                if not (min_peptide_length <= len(peptide) <= max_peptide_length):
                    continue
                if standard_amino_acids:
                    if len(set(peptide) - set("ARNDBCEQZGHILKMFPSTWYV")) != 0:
                        continue
                if peptide not in peptides:
                    peptides[peptide] = []
                peptides[peptide].append(protein_index)
            protein_index += 1
        return proteins, peptides

    def write_proteins(self, proteins, hdf_file, **parameters):
        # TODO: Docstring
        self.logger.info(f"Writing proteins to {self.file_name}")
        columns = sorted(
            {
                v for values in proteins.values() for v in values
            }
        )
        protrec = pd.DataFrame(
            [
                tuple(
                    [protein_name] + [
                        protein_info[column] if column in protein_info else "" for column in columns
                    ]
                ) for protein_name, protein_info in sorted(
                    proteins.items(),
                    key=lambda p: p[1]["index"]
                )
            ],
            columns=["protein"] + columns
        )
        protrec.set_index("index", inplace=True)
        self.write_dataframe_to_hdf_group_name(hdf_file, "proteins", protrec)

    def write_peptides(
        self,
        hdf_file,
        peptides,
        variable_ptms,
        fixed_ptms,
        ptm_dict,
        **kwargs,
    ):
        # TODO: Docstring
        self.logger.info(f"Writing peptidoforms to {self.file_name}")
        peptide_list = [
            (
                peptide,
                pyteomics.mass.fast_mass(peptide),
                ";".join([str(p) for p in protein_list])
            ) for (peptide, protein_list) in sorted(peptides.items())
        ]
        columns = [
            "sequence",
            "mass",
            "proteins",
        ]
        if (len(variable_ptms) + len(fixed_ptms)) > 0:
            columns += ["modifications"]
            modified_peptide_list = []
            for peptide, mass, proteins in peptide_list:
                for ptm_combination in self.generate_ptm_combinations(
                    f".{peptide}.",
                    [[]] * (len(peptide) + 2),
                    variable_ptms,
                    fixed_ptms,
                    static_ptms=False
                ):
                    parsed_ptm_combination = "|".join(
                        [
                            f"{i}|{ptm_dict[ptm][0]}" for i, ptm in enumerate(ptm_combination) if ptm != ""
                        ]
                    )
                    ptm_combo_mass = np.sum(
                        [
                            ptm_dict[ptm][1] for ptm in ptm_combination if ptm != ""
                        ]
                    )
                    if parsed_ptm_combination == "":
                        parsed_ptm_combination = "-"
                    modified_peptide_list.append(
                        (
                            peptide,
                            mass + ptm_combo_mass,
                            proteins,
                            parsed_ptm_combination
                        )
                    )
            peptide_list = modified_peptide_list
        peprec = pd.DataFrame(
            peptide_list,
            columns=columns
        )
        peprec["index"] = np.arange(peprec.shape[0])
        peprec.set_index("index", inplace=True)
        self.write_dataframe_to_hdf_group_name(hdf_file, "peptides", peprec)
        return peprec

    def write_fragments(
        self,
        hdf_file,
        peprec,
        model,
        cpu_count,
        charges,
        ptm_dict,
        variable_ptms,
        fixed_ptms,
        batch_size,
        **kwargs,
    ):
        # TODO: Docstring
        self.logger.info(f"Predicting fragments")
        ms2pip_params = {
            "ms2pip": {
                "model": model,
                "frag_error": 0,
                "ptm": [
                    ",".join([str(s) for s in ptm_values]) for ptm_values in ptm_dict.values()
                ],
                "sptm": [],
                "gptm": [],
            }
        }
        for start in range(0, peprec.shape[0], batch_size):
            end = start + batch_size
            charged_peprec = peprec[["sequence"]][start:end].rename(
                columns={'sequence': 'peptide'}
            )
            charged_peprec["spec_id"] = peprec.index[start:end]
            if (len(variable_ptms) + len(fixed_ptms)) == 0:
                charged_peprec["modifications"] = "-"
            else:
                charged_peprec["modifications"] = peprec["modifications"][start:end]
            for charge in charges:
                charged_peprec["charge"] = charge
                charged_fragrec = ms2pip.ms2pipC.MS2PIP(
                    charged_peprec,
                    num_cpu=cpu_count,
                    params=ms2pip_params,
                    return_results=True,
                ).run()
                del charged_fragrec["charge"]
                charged_fragrec.set_index(["spec_id", "ion", "ionnumber", "mz"], inplace=True)
                try:
                    partial_fragrec = partial_fragrec.join(charged_fragrec)
                    partial_fragrec.rename(
                        columns={'prediction': f'prediction_charge_{charge}'},
                        inplace=True
                    )
                except NameError:
                    partial_fragrec = charged_fragrec.rename(
                        columns={'prediction': f'prediction_charge_{charge}'}
                    )
            try:
                fragrec = pd.concat([fragrec, partial_fragrec])
            except NameError:
                fragrec = partial_fragrec.copy()
            del partial_fragrec
        fragrec.reset_index(inplace=True)
        fragrec.sort_values(by="mz", inplace=True)
        fragrec.rename(
            columns={'spec_id': f'peptide_index'},
            inplace=True
        )
        fragrec["b_ion"] = fragrec["ion"] == "B"
        fragrec["y_ion"] = fragrec["ion"] == "Y"
        del fragrec["ion"]
        fragrec["index"] = np.arange(fragrec.shape[0])
        fragrec.set_index("index", inplace=True)
        self.logger.info(f"Writing fragments to {self.file_name}")
        self.write_dataframe_to_hdf_group_name(hdf_file, "fragments", fragrec)

    @staticmethod
    def generate_ptm_combinations_recursively(ptms, selected=[]):
        # TODO: Docstring
        if len(selected) == len(ptms):
            yield selected
        else:
            for ptm in ptms[len(selected)]:
                for ptm_combination in Database.generate_ptm_combinations_recursively(
                    ptms,
                    selected + [ptm]
                ):
                    yield ptm_combination

    @staticmethod
    def generate_ptm_combinations(
        sequence,
        ptms,
        variable_ptms,
        fixed_ptms,
        static_ptms=False,
        **kwargs,
    ):
        # TODO: Docstring
        local_ptms = [[] for i in sequence]
        if sequence[0] == "n":
            local_ptms[0] += ptms[0]
        if sequence[-1] == "c":
            local_ptms[-1] = ptms[-1]
        for i, ptm in enumerate(ptms[1:-1]):
            local_ptms[i + 1] += ptm
        for i, aa in enumerate(f"n{sequence[1:-1]}c"):
            if (not static_ptms) or (len(local_ptms[i]) == 0):
                if aa in variable_ptms:
                    local_ptms[i] += variable_ptms[aa]
                if aa in fixed_ptms:
                    local_ptms[i] += fixed_ptms[aa]
                else:
                    local_ptms[i].append("")
        for ptm_combination in Database.generate_ptm_combinations_recursively(local_ptms):
            yield ptm_combination

    def get_fragment_coordinates(self, dimensions=None, indices=...):
        # TODO: Docstring
        return utils.read_hdf_coordinates(self, "fragments", dimensions, indices)

    def get_peptide_coordinates(self, dimensions=None, indices=...):
        # TODO: Docstring
        return utils.read_hdf_coordinates(self, "peptides", dimensions, indices)
