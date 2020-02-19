#!python

import numpy as np
import numba


@numba.njit
def determine_isotopic_rt_difference(
    mzs,
    rts,
    isotopic_distance,
    ppm,
    best_hit_only=False
):
    mz_order = np.argsort(mzs)
    mzs_in_mz_order = mzs[mz_order]
    rts_in_mz_order = rts[mz_order]
    if isotopic_distance > 0:
        lower_limits = np.searchsorted(
            mzs_in_mz_order,
            (mzs_in_mz_order + isotopic_distance) / (1 + ppm * 10**-6),
            "left"
        )
    else:
        lower_limits = np.arange(len(mzs)) + 1
    upper_limits = np.searchsorted(
        mzs_in_mz_order,
        (mzs_in_mz_order + isotopic_distance) * (1 + ppm * 10**-6),
        "right"
    )
    if not best_hit_only:
        first_rts = np.repeat(rts_in_mz_order, upper_limits - lower_limits)
        second_rts = np.concatenate(
            [
                rts_in_mz_order[l: h] for l, h in zip(
                    lower_limits,
                    upper_limits
                )
            ]
        )
    else:
        isotopes_found = upper_limits > lower_limits
        first_rts = rts_in_mz_order[isotopes_found]
        lower_limits = lower_limits[isotopes_found]
        upper_limits = upper_limits[isotopes_found]
        second_rts = np.array(
            [
                np.min(rts_in_mz_order[l: h]) for l, h in zip(
                    lower_limits,
                    upper_limits
                )
            ]
        )
    isotope_rt_diffs = np.abs(second_rts - first_rts)
    isotope_rt_diffs, isotope_count = np.unique(
        isotope_rt_diffs,
        return_counts=True
    )
    isotope_count = np.cumsum(isotope_count)
    return isotope_rt_diffs, isotope_count


#
#
# def getPairwiseAlignment(
#     self,
#     sample1,
#     sample2,
#     max_rt_diff=.5,
#     max_dt_diff=3,
#     ppm=30,
#     batch_size=10**7,
#     filter_unique_alignments=True,
#     parallel_processes=8,
#     sample1_preloaded_data=None,
#     sample2_preloaded_data=None,
# ):
#     sample2 = self.getSample(sample2)
#     if sample1_preloaded_data is None:
#         sample1 = self.getSample(sample1)
#         mz1 = self.getArray("MZ/values", sample1)
#         dt1 = self.getArray("DT/raw", sample1)
#         rt1 = self.getArray("RT/raw", sample1)
#         mz_identifiers1 = self.getArray("MZ/identifiers", sample1)
#     else:
#         mz1 = sample1_preloaded_data["MZ"]
#         dt1 = sample1_preloaded_data["MZ"]
#         rt1 = sample1_preloaded_data["MZ"]
#         mz_identifiers1 = sample1_preloaded_data["MZ"]
#     if sample2_preloaded_data is None:
#         sample2 = self.getSample(sample2)
#         mz2 = self.getArray("MZ/values", sample2)
#         dt2 = self.getArray("DT/raw", sample2)
#         rt2 = self.getArray("RT/raw", sample2)
#         mz_identifiers2 = self.getArray("MZ/identifiers", sample2)
#     else:
#         mz2 = sample2_preloaded_data["MZ"]
#         dt2 = sample2_preloaded_data["MZ"]
#         rt2 = sample2_preloaded_data["MZ"]
#         mz_identifiers2 = sample2_preloaded_data["MZ"]
#     max_mz_diff = (10**6 + ppm) / 10**6
#     low_limits = np.searchsorted(
#         mz2,
#         mz1 / max_mz_diff,
#         "left"
#     )
#     high_limits = np.searchsorted(
#         mz2,
#         mz1 * max_mz_diff,
#         "right"
#     )
#     batches = np.cumsum(high_limits - low_limits)
#     batch_size /= parallel_processes
#     batches = np.searchsorted(
#         batches,
#         np.arange(0, batches[-1] + 1 + batch_size, batch_size)
#     )
#     batch_iterations = zip(batches[:-1], batches[1:])
#     pairs = np.concatenate(
#         [
#             batch_pairs for batch_pairs in parallelizedGenerator(
#                 parallel_processes,
#                 batchedPairwiseAlign,
#                 batch_iterations,
#                 mz_identifiers1,
#                 mz_identifiers2,
#                 high_limits,
#                 low_limits,
#                 dt1,
#                 dt2,
#                 rt1,
#                 rt2,
#                 max_dt_diff,
#                 max_rt_diff,
#             )
#         ]
#     )
#     if filter_unique_alignments:
#         occurences = np.bincount(pairs[:, 0])
#         good = np.isin(pairs[:, 0], np.flatnonzero(occurences == 1))
#         occurences = np.bincount(pairs[:, 1])
#         good &= np.isin(pairs[:, 1], np.flatnonzero(occurences == 1))
#         pairs = pairs[good]
#     return pairs


#
# def batchedPairwiseAlign(
#     batch_iteration,
#     mz_identifiers1,
#     mz_identifiers2,
#     high_limits,
#     low_limits,
#     dt1,
#     dt2,
#     rt1,
#     rt2,
#     max_dt_diff,
#     max_rt_diff,
# ):
#     batch_low, batch_high = batch_iteration
#     ind1 = np.concatenate(
#         [
#             np.repeat(i, h) for i, h in zip(
#                 mz_identifiers1[batch_low: batch_high],
#                 high_limits[batch_low: batch_high] - low_limits[batch_low: batch_high]
#             )
#         ]
#     )
#     ind2 = np.concatenate(
#         [
#             mz_identifiers2[l: h] for l, h in zip(
#                 low_limits[batch_low: batch_high],
#                 high_limits[batch_low: batch_high]
#             )
#         ]
#     )
#     good = np.abs(dt1[ind1] - dt2[ind2]) < max_dt_diff
#     good &= np.abs(rt1[ind1] - rt2[ind2]) < max_rt_diff
#     batch_pairs = np.stack([ind1[good], ind2[good]]).T
#     return batch_pairs
# #
#
# def align(mz1, rt1, dt1, mz2, rt2, dt2, ppm, max_rt_diff, max_dt_diff):
#     @numba.njit
#     def numba_wrapper():
#         max_mz_diff = (1 + ppm * 10**-6)
#         mz1_order = np.argsort(mz1)
#         mz2_order = np.argsort(mz2)
#         low_limits = np.searchsorted(
#             mz2[mz2_order],
#             mz1[mz1_order] / max_mz_diff,
#             "left"
#         )
#         high_limits = np.searchsorted(
#             mz2[mz2_order],
#             mz1[mz1_order] * max_mz_diff,
#             "right"
#         )
#         rt2o = rt2[mz2_order]
#         dt2o = dt2[mz2_order]
#         results = []
#         for rt, dt, low, high in zip(
#             rt1[mz1_order],
#             dt1[mz1_order],
#             low_limits,
#             high_limits
#         ):
#             rtd = rt2o[low: high] - rt
#             dtd = dt2o[low: high] - dt
#             good = (np.abs(rtd) < max_rt_diff) & (np.abs(dtd) < max_dt_diff)
#             matches = mz2_order[low + np.flatnonzero(good)]
#             results.append(matches)
#         return results, mz1_order
#     results, mz1_order = numba_wrapper()
#     first_indices = np.repeat(mz1_order, [len(l) for l in results])
#     second_indices = np.concatenate(results)
#     return first_indices, second_indices


def align_edges(self, other, alignment):
    self_edges = self.get_edges()
    other_edges = other.get_edges()
    alignment_matrix = alignment.get_node_alignment(
        self,
        other,
        return_as_scipy_csr=True
    )
    indirect = (ali.T * edges1 * ali).tocsr()
    positive = indirect.multiply(edges2)
    available = edges2 * ali.T * ali
    negative = available - positive
    return self_edges, other_edges, alignment_matrix


    def getSecondaryCoelutionEvidencePairs(
        self,
        sample,
        evidence_sample,
        coelution_matrix=None,
        reorder_by_rt=True
    ):
        if coelution_matrix is None:
            coelution_matrix = createMatrixFromPairedArray(
                self.getArray(
                    "coeluting_ions",
                    evidence_sample,
                ),
                (
                    self.getSampleSize(evidence_sample),
                    self.getSampleSize(evidence_sample)
                )
            )
        alignment_matrix = createMatrixFromPairedArray(
            self.getArray(
                "pairwise/alignments",
                sample,
                evidence_sample,
            ),
            (
                self.getSampleSize(sample),
                self.getSampleSize(evidence_sample),
            )
        )
        result = alignment_matrix * coelution_matrix * alignment_matrix.T
        pairs = np.stack(result.nonzero()).T
        if reorder_by_rt:
            rt = self.getArray("RT/raw", sample)
            pairs = np.concatenate(
                [pairs, pairs[:, ::-1]]
            )
            pairs = pairs[np.diff(rt[pairs], axis=1).flatten() > 0]
        return pairs

    def getAvailableCoelutionEvidencePairs(
        self,
        sample,
        evidence_sample,
        coelution_matrix=None,
    ):
        if coelution_matrix is None:
            coelution_matrix = createMatrixFromPairedArray(
                self.getArray(
                    "coeluting_ions",
                    sample,
                ),
                (
                    self.getSampleSize(sample),
                    self.getSampleSize(sample)
                )
            )
        alignment_matrix = createMatrixFromPairedArray(
            self.getArray(
                "pairwise/alignments",
                sample,
                evidence_sample,
            ),
            (
                self.getSampleSize(sample),
                self.getSampleSize(evidence_sample),
            )
        )
        result = (coelution_matrix * alignment_matrix) * alignment_matrix.T
        pairs = np.stack(result.nonzero()).T
        return pairs

    def getEvidencedCoelutionMatrixWithOverlapHitArrays(self, sample):
        sample = self.getSample(sample)
        available_evidence_matrix = scipy.sparse.csr_matrix(
            (
                self.getSampleSize(sample),
                self.getSampleSize(sample)
            ),
            dtype=np.int
        )
        for evidence_sample in self.samples:
            if evidence_sample != sample:
                available_evidence_matrix += createMatrixFromPairedArray(
                    self.getArray(
                        "pairwise/evidence/available",
                        sample,
                        evidence_sample,
                    ),
                    (
                        self.getSampleSize(sample),
                        self.getSampleSize(sample)
                    )
                )
        secondary_evidence_matrix = scipy.sparse.csr_matrix(
            (
                self.getSampleSize(sample),
                self.getSampleSize(sample)
            ),
            dtype=np.int
        )
        for evidence_sample in self.samples:
            if evidence_sample != sample:
                secondary_evidence_matrix += createMatrixFromPairedArray(
                    self.getArray(
                        "pairwise/evidence/secondary",
                        sample,
                        evidence_sample,
                    ),
                    (
                        self.getSampleSize(sample),
                        self.getSampleSize(sample)
                    )
                )
        hits = secondary_evidence_matrix.copy().multiply(
            available_evidence_matrix.astype(np.bool)
        )
        overlap = available_evidence_matrix.copy().multiply(
            secondary_evidence_matrix.astype(np.bool)
        )
        if True:
            # target_edges = ((hits * 2) > overlap)  # & (overlap > 1)
            target_edges = (hits.data * 2 - overlap.data) > 0
            positive_evidence = hits.data[target_edges]
            available_evidence = overlap.data[target_edges]
            hits.data = target_edges
            hits.eliminate_zeros()
            hits.sort_indices()


    @functools.lru_cache
    def get_alignment(
        self,
        first_ion_network,
        second_ion_network,
        nodes=True,
        edges=True,
        return_as_scipy_csr=False
    ):
        # @functools.lru_cache
        def create_matrix(array1, array2, size1, size2, swap):
            if swap:
                matrix = scipy.sparse.csr_matrix(
                    (
                        np.ones(len(array1), dtype=np.bool),
                        (array2, array1)
                    ),
                    shape=(size2, size1)
                )
            else:
                matrix = scipy.sparse.csr_matrix(
                    (
                        np.ones(len(array1), dtype=np.bool),
                        (array1, array2)
                    ),
                    shape=(size1, size2)
                )
            return matrix
        results = []
        if first_ion_network > second_ion_network:
            first_ion_network, second_ion_network = second_ion_network, first_ion_network
            swap = True
        else:
            swap = False
        with h5py.File(self.file_name, "r") as alignment_file:
            alignment_group = alignment_file[
                first_ion_network.key
            ][
                second_ion_network.key
            ]
            if nodes:
                nodes_first_indices = alignment_group["nodes/first_indices"][...]
                nodes_second_indices = alignment_group["nodes/second_indices"][...]
                if return_as_scipy_csr:
                    node_result = create_matrix(
                        nodes_first_indices,
                        nodes_second_indices,
                        first_ion_network.node_count,
                        second_ion_network.node_count,
                        swap
                    )
                else:
                    if swap:
                        node_result = np.stack(
                            [nodes_second_indices, nodes_first_indices]
                        ).T
                    else:
                        node_result = np.stack(
                            [nodes_first_indices, nodes_second_indices]
                        ).T
                results.append(node_result)
            if edges:
                edges_first_indices = alignment_group["edges/first_indices"][...]
                edges_second_indices = alignment_group["edges/second_indices"][...]
                if return_as_scipy_csr:
                    edge_result = create_matrix(
                        edges_first_indices,
                        edges_second_indices,
                        first_ion_network.edge_count,
                        second_ion_network.edge_count,
                        swap
                    )
                else:
                    if swap:
                        edge_result = np.stack(
                            [edges_second_indices, edges_first_indices]
                        ).T
                    else:
                        edge_result = np.stack(
                            [edges_first_indices, edges_second_indices]
                        ).T
                results.append(edge_result)
        return results


    # def get_evidence(
    #     self,
    #     network_keys=None,
    #     kind=["positive_edges", "negative_edges", "nodes"],
    #     return_total=False
    # ):
    #     # TODO: Docstring
    #     if network_keys is None:
    #         network_keys = self.network_keys
    #     if isinstance(network_keys, Network):
    #         # FIXME: Does not work?
    #         network_keys = [network_keys.key]
    #     elif isinstance(network_keys, str):
    #         network_keys = [network_keys]
    #     else:
    #         network_keys = [
    #             network_key if isinstance(
    #                 network_key, str
    #             ) else network_key.key for network_key in network_keys
    #         ]
    #     arrays = {
    #         "positive_edges": [],
    #         "negative_edges": [],
    #         "nodes": [],
    #     }
    #     single_kind = False
    #     if isinstance(kind, str):
    #         kind = [kind]
    #         single_kind = True
    #     with h5py.File(self.file_name, "r") as evidence_file:
    #         for network_key in network_keys:
    #             if "positive_edges" in kind:
    #                 arrays["positive_edges"].append(
    #                     evidence_file[network_key]["positive_edge_mask"][...]
    #                 )
    #             if "negative_edges" in kind:
    #                 arrays["negative_edges"].append(
    #                     evidence_file[network_key]["negative_edge_mask"][...]
    #                 )
    #             if "nodes" in kind:
    #                 arrays["nodes"].append(
    #                     evidence_file[network_key]["node_mask"][...]
    #                 )
    #     if return_total:
    #         # TODO: earlier summation saves memory!
    #         arrays = {
    #             key: np.sum(value, axis=0) for key, value in arrays.items()
    #         }
    #     if (len(network_keys) == 1) or return_total:
    #         if return_total:
    #             if single_kind:
    #                 return arrays[kind[0]]
    #             return tuple(arrays[k] for k in kind)
    #         if single_kind:
    #             return arrays[kind[0]][0]
    #         return tuple(arrays[k][0] for k in kind)
    #     else:
    #         if single_kind:
    #             return arrays[kind[0]]
    #         return tuple(arrays[k] for k in kind)


@numba.njit(fastmath=True)
def longest_increasing_subsequence(sequence):
    M = np.repeat(0, len(sequence) + 1)
    P = np.repeat(0, len(sequence))
    max_subsequence_length = 0
    for current_index, current_element in enumerate(sequence):
        low_bound = 1
        high_bound = max_subsequence_length
        while low_bound <= high_bound:
            mid = (low_bound + high_bound) // 2
            if sequence[M[mid]] <= current_element:
                low_bound = mid + 1
            else:
                high_bound = mid - 1
        subsequence_length = low_bound
        P[current_index] = M[subsequence_length - 1]
        M[subsequence_length] = current_index
        if subsequence_length > max_subsequence_length:
            max_subsequence_length = subsequence_length
    longest_increasing_subsequence = np.repeat(0, max_subsequence_length)
    index = M[max_subsequence_length]
    for current_index in range(max_subsequence_length - 1, -1, -1):
        longest_increasing_subsequence[current_index] = index
        index = P[index]
    # good = np.repeat(True, max_subsequence_length)
    # good[1:] = sequence[longest_increasing_subsequence[:-1]] != sequence[longest_increasing_subsequence[1:]]
    # return longest_increasing_subsequence[good]
    return longest_increasing_subsequence


def quick_align(self, other, ppm):
    self_mzs = self.get_ion_coordinates("FRAGMENT_MZ")
    other_mzs = other.get_ion_coordinates("FRAGMENT_MZ")
    self_mz_order = np.argsort(self_mzs)
    other_mz_order = np.argsort(other_mzs)
    max_mz_diff = 1 + ppm * 10**-6
    low_limits = np.searchsorted(
        self_mzs[self_mz_order],
        other_mzs[other_mz_order] / max_mz_diff,
        "left"
    )
    high_limits = np.searchsorted(
        self_mzs[self_mz_order],
        other_mzs[other_mz_order] * max_mz_diff,
        "right"
    )
    other_rt_order = np.argsort(other_mz_order)
    self_indices = np.concatenate(
        [
            self_mz_order[l:h] for l, h in zip(
                low_limits[other_rt_order],
                high_limits[other_rt_order]
            )
        ]
    )
    other_indices = np.repeat(
        np.arange(len(other_rt_order)),
        high_limits[other_rt_order] - low_limits[other_rt_order]
    )
    selection = longest_increasing_subsequence(self_indices)
    self_indices_mask = np.empty(len(selection) + 2, dtype=int)
    self_indices_mask[0] = 0
    self_indices_mask[1: -1] = self_indices[selection]
    self_indices_mask[-1] = len(self_mzs) - 1
    other_indices_mask = np.empty(len(selection) + 2, dtype=int)
    other_indices_mask[0] = 0
    other_indices_mask[1: -1] = other_indices[selection]
    other_indices_mask[-1] = len(other_mzs) - 1
    return self_indices_mask, other_indices_mask


def calibrate_precursor_rt(self, other, ppm=10):
    self_indices, other_indices = quick_align(self, other, ppm=10)
    self_rts = self.get_ion_coordinates("PRECURSOR_RT")
    other_rts = other.get_ion_coordinates("PRECURSOR_RT", indices=other_indices)
    new_self_rts = []
    for self_start_index, self_end_index, other_rt_start, other_rt_end in zip(
        self_indices[:-1],
        self_indices[1:],
        other_rts[:-1],
        other_rts[1:]
    ):
        self_rt_start = self_rts[self_start_index]
        self_rt_end = self_rts[self_end_index]
        if self_rt_start == self_rt_end:
            new_rts = np.repeat(other_rt_start, self_end_index - self_start_index)
        else:
            slope = (other_rt_end - other_rt_start) / (self_rt_end - self_rt_start)
            new_rts = other_rt_start + slope * (
                self_rts[self_start_index: self_end_index] - self_rt_start
            )
        new_self_rts.append(new_rts)
    new_self_rts.append([other_rts[-1]])
    new_self_rts = np.concatenate(new_self_rts)
    return new_self_rts
