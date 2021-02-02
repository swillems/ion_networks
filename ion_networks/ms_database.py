#!python

# external
import pandas as pd
import numpy as np
import scipy.signal
try:
    import ms2pip.ms2pipC
    import ms2pip.retention_time
except ModuleNotFoundError:
    print(
        "MS2PIP not found, see "
        "https://github.com/compomics/ms2pip_c or "
        "https://github.com/compomics/ms2pip_c/issues/19 "
        "for manual installation."
    )
import pyteomics.parser
import pyteomics.fasta
# local
from ion_networks import ms_utils


class HDF_Database_File(ms_utils.HDF_File):
    # TODO: Docstring

    def create_from_fastas(
        self,
        fasta_file_names,
        parameters,
    ):
        # TODO: Docstring
        ms_utils.LOGGER.info(f"Creating database {self.file_name}")
        proteins, peptides = self.read_proteins_and_peptides_from_fasta(
            fasta_file_names,
            **parameters
        )
        self.write_proteins(proteins, **parameters)
        peprec, mod_mass_list = self.write_peptides(peptides, **parameters)
        self.write_fragments(peprec, mod_mass_list, **parameters)
        self.write_parameters(fasta_file_names, parameters)

    def write_parameters(self, fasta_file_names, parameters):
        # TODO: Docstring
        ms_utils.LOGGER.info(f"Writing parameters to {self.file_name}")
        self.write_attr("fasta_file_names", fasta_file_names)
        for key, value in parameters.items():
            self.write_attr(key, value)

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
            ms_utils.LOGGER.info(f"Reading {fasta_file_name}")
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
                    protein_index=len(proteins),
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
        protein_index=0,
        **kwargs,
    ):
        # TODO: Docstring
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

    def write_proteins(self, proteins, **parameters):
        # TODO: Docstring
        ms_utils.LOGGER.info(f"Writing proteins to {self.file_name}")
        columns = sorted(
            {
                v for values in proteins.values() for v in values
            }
        )
        protrec = pd.DataFrame(
            [
                tuple(
                    [
                        protein_name,
                        protein_name.startswith(
                            parameters["decoy_prefix"]
                        ) & protein_name.endswith(
                            parameters["decoy_suffix"]
                        )
                    ] + [
                        protein_info[column] if column in protein_info else "" for column in columns
                    ]
                ) for protein_name, protein_info in sorted(
                    proteins.items(),
                    key=lambda p: p[1]["index"]
                )
            ],
            columns=["protein", "decoy"] + columns
        )
        protrec.set_index("index", inplace=True)
        self.write_dataset("proteins", protrec)

    def write_peptides(
        self,
        peptides,
        variable_ptms,
        fixed_ptms,
        ptm_dict,
        **kwargs,
    ):
        # TODO: Docstring
        ms_utils.LOGGER.info(f"Writing peptidoforms to {self.file_name}")
        decoys = self.read_dataset("decoy", parent_group_name="proteins")
        peptide_list = [
            (
                peptide,
                pyteomics.mass.fast_mass(peptide),
                ";".join([str(p) for p in protein_list]),
                np.all(decoys[protein_list])
            ) for (peptide, protein_list) in sorted(peptides.items())
        ]
        columns = [
            "sequence",
            "mass",
            "proteins",
            "decoy",
        ]
        if (len(variable_ptms) + len(fixed_ptms)) > 0:
            columns += ["modifications"]
            modified_peptide_list = []
            mod_mass_list = []
            for peptide, mass, proteins, decoy in peptide_list:
                for ptm_combination in self.generate_ptm_combinations(
                    f".{peptide}.",
                    [[]] * (len(peptide) + 2),
                    variable_ptms,
                    fixed_ptms,
                    static_ptms=False
                ):
                    parsed_ptm_combination = "|".join(
                        [
                            f"{i}|{ptm_dict[ptm][0]}" for i, ptm in enumerate(
                                ptm_combination
                            ) if ptm != ""
                        ]
                    )
                    mod_masses = [
                        ptm_dict[ptm][1] if ptm != "" else 0 for ptm in ptm_combination
                    ]
                    ptm_combo_mass = np.sum(mod_masses)
                    if parsed_ptm_combination == "":
                        parsed_ptm_combination = "-"
                    modified_peptide_list.append(
                        (
                            peptide,
                            mass + ptm_combo_mass,
                            proteins,
                            decoy,
                            parsed_ptm_combination
                        )
                    )
                    mod_mass_list.append(mod_masses)
            peptide_list = modified_peptide_list
        peprec = pd.DataFrame(
            peptide_list,
            columns=columns
        )
        peprec["index"] = np.arange(peprec.shape[0])
        peprec.set_index("index", inplace=True)
        self.write_dataset("peptides", peprec)
        return peprec, mod_mass_list

    def write_fragments(
        self,
        peprec,
        mod_mass_list,
        model,
        charges,
        ptm_dict,
        variable_ptms,
        fixed_ptms,
        batch_size,
        predict_intensities,
        **kwargs,
    ):
        # TODO: Docstring
        if predict_intensities:
            fragrec = self.predict_with_ms2pip(
                peprec,
                model,
                charges,
                ptm_dict,
                variable_ptms,
                fixed_ptms,
                batch_size,
                **kwargs,
            )
        else:
            fragrec = self.fragment_peptides(
                peprec,
                ptm_dict,
                variable_ptms,
                fixed_ptms,
                mod_mass_list,
                **kwargs,
            )
        ms_utils.LOGGER.info(f"Writing fragments to {self.file_name}")
        self.write_dataset("fragments", fragrec)

    @staticmethod
    def fragment_peptides(
        peprec,
        ptm_dict,
        variable_ptms,
        fixed_ptms,
        mod_mass_list,
        **kwargs,
    ):
        ms_utils.LOGGER.info("Generating fragments")
        mzs = []
        peptide_indices = []
        ion_numbers = []
        for peptide_index, (peptide, mods) in enumerate(
            zip(
                peprec["sequence"],
                mod_mass_list,
            )
        ):
            peptide_indices += [peptide_index] * (len(peptide) - 1) * 2
            mod_masses = np.array(mods)
            nt_mod_masses = np.cumsum(mod_masses)
            ct_mod_masses = np.cumsum(mod_masses[::-1])[::-1]
            for ion_number in range(1, len(peptide)):
                mz = pyteomics.mass.fast_mass(
                    peptide[:ion_number],
                    ion_type="b",
                    charge=1
                ) + nt_mod_masses[ion_number]
                mzs.append(mz)
                mz = pyteomics.mass.fast_mass(
                    peptide[ion_number:],
                    ion_type="y",
                    charge=1
                ) + ct_mod_masses[ion_number + 1]
                mzs.append(mz)
                ion_numbers += [ion_number, len(peptide) - ion_number]
        b_ions = np.array([True, False] * (len(mzs) // 2))
        prediction_charge_2 = np.ones(len(mzs))
        fragrec = pd.DataFrame(
            {
                "b_ion": b_ions,
                "ionnumber": ion_numbers,
                "mz": mzs,
                "peptide_index": peptide_indices,
                "prediction_charge_2": prediction_charge_2,
                "y_ion": ~b_ions,
            }
        )
        return fragrec.sort_values(by=['mz'])

    @staticmethod
    def predict_with_ms2pip(
        peprec,
        model,
        charges,
        ptm_dict,
        variable_ptms,
        fixed_ptms,
        batch_size,
        **kwargs,
    ):
        ms_utils.LOGGER.info("Predicting fragments with MS2PIP")
        ms2pip_params = {
            "ms2pip": {
                "model": model,
                "frag_error": 0,
                "ptm": [
                    ",".join(
                        [str(s) for s in ptm_values]
                    ) for ptm_values in ptm_dict.values()
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
                    num_cpu=ms_utils.MAX_THREADS,
                    params=ms2pip_params,
                    return_results=True,
                ).run()
                del charged_fragrec["charge"]
                charged_fragrec.set_index(
                    ["spec_id", "ion", "ionnumber", "mz"],
                    inplace=True
                )
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
            columns={'spec_id': 'peptide_index'},
            inplace=True
        )
        fragrec["b_ion"] = fragrec["ion"] == "B"
        fragrec["y_ion"] = fragrec["ion"] == "Y"
        del fragrec["ion"]
        fragrec["index"] = np.arange(fragrec.shape[0])
        fragrec.set_index("index", inplace=True)
        return fragrec

    @staticmethod
    def generate_ptm_combinations_recursively(ptms, selected=[]):
        # TODO: Docstring
        if len(selected) == len(ptms):
            yield selected
        else:
            for ptm in ptms[len(selected)]:
                for ptm_combination in HDF_Database_File.generate_ptm_combinations_recursively(
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
        for ptm_combination in HDF_Database_File.generate_ptm_combinations_recursively(local_ptms):
            yield ptm_combination

    def get_fragment_coordinates(self, dimensions=None, indices=...):
        # TODO: Docstring
        if isinstance(dimensions, str):
            return self.read_dataset(
                dimensions,
                parent_group_name="fragments",
                indices=indices
            )
        elif dimensions is None:
            dimensions = self.read_group(
                parent_group_name="fragments"
            )
        return [
            self.read_dataset(
                dimension,
                parent_group_name="fragments",
                indices=indices
            ) for dimension in dimensions
        ]

    def get_peptide_coordinates(self, dimensions=None, indices=...):
        # TODO: Docstring
        if isinstance(dimensions, str):
            return self.read_dataset(
                dimensions,
                parent_group_name="peptides",
                indices=indices
            )
        elif dimensions is None:
            dimensions = self.read_group(
                parent_group_name="peptides"
            )
        return [
            self.read_dataset(
                dimension,
                parent_group_name="peptides",
                indices=indices
            ) for dimension in dimensions
        ]

    def get_db_targets(
        self,
        max_ppm=100,
        min_distance=0.5,
        ms_level=2,
    ):
        if ms_level == 1:
            db_mzs_ = self.read_dataset(
                "mass",
                parent_group_name="peptides",
            )
        elif ms_level == 2:
            db_mzs_ = self.read_dataset(
                "mz",
                parent_group_name="fragments",
            )
        else:
            raise ValueError(f"{ms_level} is not a valid ms level")
        tmp_result = np.bincount(
            (
                np.log10(
                    db_mzs_[
                        np.isfinite(db_mzs_) & (db_mzs_ > 0)
                    ].flatten()
                ) * 10**6
            ).astype(np.int64)
        )
        db_mz_distribution = np.zeros_like(tmp_result)
        for i in range(1, max_ppm):
            db_mz_distribution[i:] += tmp_result[:-i]
            db_mz_distribution[:-i] += tmp_result[i:]
        peaks = scipy.signal.find_peaks(db_mz_distribution, distance=max_ppm)[0]
        db_targets = 10 ** (peaks / 10**6)
    #     db_vals = db_mz_distribution[peaks]
    #     plt.vlines(db_targets, 0, db_vals)
        db_array = np.zeros(int(db_targets[-1]) + 1, dtype=np.float64)
        last_int_mz = -1
        last_mz = -1
        for mz in db_targets:
            mz_int = int(mz)
            if (mz_int != last_int_mz) & (mz > (last_mz + min_distance)):
                db_array[mz_int] = mz
            else:
                db_array[mz_int] = 0
            last_int_mz = mz_int
            last_mz = mz
        return db_array

    def align_mz_values(
        self,
        mzs,
        rts,
        max_ppm_distance=np.inf,
        rt_step_size=0.1,
    ):
        selected = mzs.astype(np.int64)
        ds = np.zeros((3, len(selected)))
        db_array = self.get_db_targets()
        if len(db_array) < len(selected) + 1:
            tmp = np.zeros(len(selected) + 1)
            tmp[:len(db_array)] = db_array
            db_array = tmp
        ds[0] = mzs - db_array[selected - 1]
        ds[1] = mzs - db_array[selected]
        ds[2] = mzs - db_array[selected + 1]
        min_ds = np.take_along_axis(
            ds,
            np.expand_dims(np.argmin(np.abs(ds), axis=0), axis=0),
            axis=0
        ).squeeze(axis=0)
        ppm_ds = min_ds / mzs * 10**6
        selected = np.abs(ppm_ds) < max_ppm_distance
        selected &= np.isfinite(rts)
        rt_order = np.argsort(rts)
        rt_order = rt_order[selected[rt_order]]
        ordered_rt = rts[rt_order]
        ordered_ppm = ppm_ds[rt_order]
        rt_idx_break = np.searchsorted(
            ordered_rt,
            np.arange(ordered_rt[0], ordered_rt[-1], rt_step_size),
            "left"
        )
        median_ppms = np.empty(len(rt_idx_break) - 1)
        for i in range(len(median_ppms)):
            median_ppms[i] = np.median(
                ordered_ppm[rt_idx_break[i]: rt_idx_break[i + 1]]
            )
        estimated_errors = scipy.interpolate.griddata(
            rt_step_size / 2 + np.arange(
                ordered_rt[0],
                ordered_rt[-1] - 2 * rt_step_size,
                rt_step_size
            ),
            median_ppms,
            rts,
            fill_value=0,
            method="linear",
            rescale=True
        )
        estimated_errors[~np.isfinite(estimated_errors)] = 0
        return mzs * (1 + estimated_errors / 10**6)
