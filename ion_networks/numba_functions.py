#!python

# external
import numpy as np
import numba


@numba.njit(nogil=True, cache=True)
def longest_increasing_subsequence(sequence):
    # TODO:Docstring
    M = np.zeros(len(sequence) + 1, np.int64)
    P = np.zeros(len(sequence), np.int64)
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
    return longest_increasing_subsequence


@numba.njit(nogil=True, cache=True)
def increase_buffer(buffer, max_batch=10**7):
    new_buffer = np.empty(buffer.shape[0] + max_batch, np.int64)
    new_buffer[:len(buffer)] = buffer
    return new_buffer


@numba.njit(nogil=True, cache=True)
def quick_align(
    self_mzs,
    other_mzs,
    self_mz_order,
    other_mz_order,
    other_rt_order,
    ppm
):
    # TODO: Docstring
    max_mz_diff = 1 + ppm * 10**-6
    low_limits = np.searchsorted(
        self_mzs[self_mz_order],
        other_mzs[other_mz_order] / max_mz_diff,
        "left"
    )[other_rt_order]
    high_limits = np.searchsorted(
        self_mzs[self_mz_order],
        other_mzs[other_mz_order] * max_mz_diff,
        "right"
    )[other_rt_order]
    diffs = high_limits - low_limits
    ends = np.cumsum(diffs)
    self_indices = np.empty(ends[-1], np.int64)
    for l, h, e, d in zip(low_limits, high_limits, ends, diffs):
        self_indices[e - d: e] = self_mz_order[l: h]
    selection = longest_increasing_subsequence(self_indices)
    self_indices_mask = np.empty(len(selection) + 2, np.int64)
    self_indices_mask[0] = 0
    self_indices_mask[1: -1] = self_indices[selection]
    self_indices_mask[-1] = len(self_mzs) - 1
    other_indices_mask = np.empty(len(selection) + 2, np.int64)
    other_indices_mask[0] = 0
    other_indices = np.repeat(
        np.arange(len(other_rt_order)),
        high_limits - low_limits
    )
    other_indices_mask[1: -1] = other_indices[selection]
    other_indices_mask[-1] = len(other_mzs) - 1
    return self_indices_mask, other_indices_mask


@numba.njit(nogil=True, cache=True)
def align_coordinates(
    queries,
    lower_limits,
    upper_limits,
    self_coordinates,
    other_coordinates,
    max_errors,
    # kind="euclidean"
):
    indptr = np.zeros(len(queries), np.int64)
    indices = np.empty(10**7, np.int64)
    total = 0
    for index, query in enumerate(queries):
        low_limit = lower_limits[query]
        high_limit = upper_limits[query]
        candidate_count = high_limit - low_limit
        if candidate_count == 0:
            continue
        elif (candidate_count + total) >= len(indices):
            indices = increase_buffer(indices)
        dists = other_coordinates[low_limit: high_limit] - self_coordinates[query]
        # TODO: what if error==0?
        # if kind == "euclidean":
        dists /= max_errors
        dists = dists**2
        projected_dists = np.sum(dists, axis=1)
        projected_dists = np.sqrt(projected_dists)
        candidates = low_limit + np.flatnonzero(projected_dists <= 1)
        # elif kind == "manhattan":
        #     projected_dists = np.all(dists < max_errors, axis=1)
        #     candidates = low_limit + np.flatnonzero(projected_dists)
        candidate_count = len(candidates)
        indices[total: total + candidate_count] = candidates
        indptr[index] = candidate_count
        total += candidate_count
    return (indptr, indices[:total])


@numba.njit(nogil=True, cache=True)
def make_symmetric(
    indptr,
    indices,
):
    # TODO: multithread?
    offsets = np.cumsum(np.bincount(indices))
    indptr_ = indptr.copy()
    indptr_[1:1 + offsets.shape[0]] += offsets
    indptr_[1 + offsets.shape[0]:] += offsets[-1]
    indices_ = np.empty(indptr_[-1], np.int64)
    pointers_ = np.empty_like(indices_)
    offsets = indptr_[:-1] + np.diff(indptr)
    for index in range(indptr.shape[0] - 1):
        start = indptr[index]
        end = indptr[index + 1]
        current_indices = indices[start: end]
        pointers = np.arange(start, end)
        start_ = indptr_[index]
        end_ = start_ + current_indices.shape[0]
        indices_[start_: end_] = current_indices
        pointers_[start_: end_] = pointers
        current_offsets = offsets[current_indices]
        indices_[current_offsets] = index
        pointers_[current_offsets] = pointers
        offsets[current_indices] += 1
    return indptr_, indices_, pointers_


@numba.njit(nogil=True, cache=True)
def align_edges(
    queries,
    self_indptr,
    self_indices,
    self_pointers,
    other_indptr,
    other_indices,
    alignment,
    alignment_mask,
):
    self_pointers_ = np.empty(10**7, np.int64)
    other_pointers_ = np.empty(10**7, np.int64)
    pointer_offset = 0
    for index in queries:
        possible_start = self_indptr[index]
        possible_end = self_indptr[index + 1]
        if possible_start == possible_end:
            continue
        current_index = alignment[index]
        current_start = other_indptr[current_index]
        current_end = other_indptr[current_index + 1]
        if current_start == current_end:
            continue
        possible_indices = self_indices[possible_start: possible_end]
        possible_mask = alignment_mask[possible_indices]
        if not np.any(possible_mask):
            continue
        possible_indices = alignment[possible_indices[possible_mask]]
        possible_pointers = self_pointers[possible_start: possible_end][
            possible_mask
        ]
        current_indices = other_indices[current_start: current_end]
        candidates1 = np.searchsorted(
            current_indices,
            possible_indices,
            "left"
        )
        candidates2 = np.searchsorted(
            current_indices,
            possible_indices,
            "right"
        )
        overlap = np.flatnonzero(candidates2 != candidates1)
        overlap_count = len(overlap)
        if len(overlap) == 0:
            continue
        elif (overlap_count + pointer_offset) >= len(self_pointers_):
            self_pointers_ = increase_buffer(self_pointers_)
            other_pointers_ = increase_buffer(other_pointers_)
        self_pointers_[
            pointer_offset: pointer_offset + overlap_count
        ] = possible_pointers[overlap]
        other_pointers_[
            pointer_offset: pointer_offset + overlap_count
        ] = current_start + candidates1[overlap]
        pointer_offset += overlap_count
    return self_pointers_[:pointer_offset], other_pointers_[:pointer_offset]


@numba.njit(cache=True)
def find_peak_indices(
    input_array,
    output_array,
    max_distance,
):
    peaks = np.zeros(int(input_array[-1]), np.int64)
    current_max_mz = 0
    current_max_int = 0
    current_max_index = 0
    for index, (intensity, mz) in enumerate(zip(output_array, input_array)):
        if mz > current_max_mz + max_distance:
            peaks[int(current_max_mz)] = current_max_index
            current_max_mz = mz
            current_max_int = intensity
            current_max_index = index
        elif intensity > current_max_int:
            current_max_mz = mz
            current_max_int = intensity
            current_max_index = index
    return peaks


@numba.njit(nogil=True, cache=True)
def get_unique_apex_and_count(
    ordered_self_indices,
    ordered_other_indices,
    return_all_counts=True
):
    counts = np.zeros_like(ordered_self_indices)
    self_max = np.max(ordered_self_indices)
    other_max = np.max(ordered_other_indices)
    unique_pair = np.zeros(counts.shape[0], np.bool_)
    self_frequencies = np.zeros(self_max + 1, np.int64)
    other_frequencies = np.zeros(other_max + 1, np.int64)
    self_indptr = np.empty(self_max + 2, np.int64)
    self_indptr[0] = 0
    self_indptr[1:] = np.cumsum(np.bincount(ordered_self_indices))
    self_order = np.argsort(ordered_self_indices)
    other_indptr = np.empty(other_max + 2, np.int64)
    other_indptr[0] = 0
    other_indptr[1:] = np.cumsum(np.bincount(ordered_other_indices))
    other_order = np.argsort(ordered_other_indices)
    unique_count = 0
    max_count = 0
    apex = 0
    for i in range(counts.shape[0]):
        self_index = ordered_self_indices[i]
        other_index = ordered_other_indices[i]
        if (
            self_frequencies[self_index] == 0
        ) & (
            other_frequencies[other_index] == 0
        ):
            unique_count += 1
            unique_pair[i] = True
            if unique_count > max_count:
                apex = i
                max_count = unique_count
        else:
            self_locs = self_order[
                self_indptr[self_index]: self_indptr[self_index + 1]
            ]
            if np.any(unique_pair[self_locs]):
                unique_count -= 1
            other_locs = other_order[
                other_indptr[other_index]: other_indptr[other_index + 1]
            ]
            if np.any(unique_pair[other_locs]):
                unique_count -= 1
            unique_pair[self_locs] = False
            unique_pair[other_locs] = False
        self_frequencies[self_index] += 1
        other_frequencies[other_index] += 1
        counts[i] = unique_count
    if not return_all_counts:
        counts = counts[apex: apex + 1]
    return apex, counts


@numba.njit
def cluster_network(
    indptr,
    indices,
    edge_pointers,
    selected_edges,
):
    node_count = indptr.shape[0] - 1
    clusters = np.zeros(node_count, np.int64)
    cluster_number = 0
    for index in range(node_count):
        if clusters[index] != 0:
            continue
        current_cluster = set()
        new_indices = set()
        new_indices.add(index)
        while len(new_indices) > 0:
            new_index = new_indices.pop()
            current_cluster.add(new_index)
            neighbors = indices[indptr[new_index]: indptr[new_index + 1]]
            pointers = edge_pointers[indptr[new_index]: indptr[new_index + 1]]
            selected = selected_edges[pointers]
            new_indices |= set(neighbors[selected]) - current_cluster
        cluster_number += 1
        for i in current_cluster:
            clusters[i] = cluster_number
    return clusters


@numba.njit()
def __get_candidate_peptide_indices_for_edges(
    indptr,
    indices,
    low_peptide_indices,
    high_peptide_indices,
    database_peptides,
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
        peptide_candidates = database_peptides[low: high]
        peptide_candidates_set = set(peptide_candidates)
        neighbors = indices[start: end]
        for i, neighbor in enumerate(neighbors):
            neighbor_low = low_peptide_indices[neighbor]
            neighbor_high = high_peptide_indices[neighbor]
            if neighbor_low == neighbor_high:
                result_indptr[start + i] = current_index
                continue
            neighbor_peptide_candidates = database_peptides[
                neighbor_low: neighbor_high
            ]
            for neighbor_peptide_candidate in neighbor_peptide_candidates:
                if neighbor_peptide_candidate in peptide_candidates_set:
                    result_indices[
                        current_index
                    ] = neighbor_peptide_candidate
                    current_index += 1
            result_indptr[start + i] = current_index
    result_indptr[1:] = result_indptr[:-1]
    result_indptr[0] = 0
    return result_indptr, result_indices[:current_index]


@numba.njit(cache=True, nogil=True)
def annotate_mgf(
    queries,
    spectra_indptr,
    low_limits,
    high_limits,
    peptide_pointers,
):
    peptide_count = np.max(peptide_pointers) + 1
    count = 0
    for s in queries:
        count += spectra_indptr[s + 1] - spectra_indptr[s]
    score_results = np.empty(count, np.float64)
    fragment_results = np.empty(count, np.int64)
    index_results = np.empty(count, np.int64)
    count_results = np.empty(count, np.int64)
    candidate_counts = np.empty(count, np.int64)
    spectrum_sizes = np.empty(count, np.int64)
    current_i = 0
    for spectrum_index in queries:
        spectrum_start = spectra_indptr[spectrum_index]
        spectrum_end = spectra_indptr[spectrum_index + 1]
        candidates = np.zeros(peptide_count, np.int64)
        for ion_index in range(spectrum_start, spectrum_end):
            peptide_low = low_limits[ion_index]
            peptide_high = high_limits[ion_index]
            if peptide_low == peptide_high:
                continue
            peptides = peptide_pointers[peptide_low: peptide_high]
            candidates[peptides] += 1
        for ion_index in range(spectrum_start, spectrum_end):
            peptide_low = low_limits[ion_index]
            peptide_high = high_limits[ion_index]
            if peptide_low == peptide_high:
                continue
            (
                score,
                max_count,
                max_fragment,
                candidate_count
            ) = score_regression_estimator(
                candidates[peptide_pointers[peptide_low: peptide_high]],
                peptide_low,
                peptide_count
            )
            if score > 0:
                score_results[current_i] = score
                fragment_results[current_i] = max_fragment
                index_results[current_i] = ion_index
                count_results[current_i] = max_count
                candidate_counts[current_i] = candidate_count
                spectrum_sizes[current_i] = spectrum_end - spectrum_start
                current_i += 1
    return (
        score_results[:current_i],
        fragment_results[:current_i],
        index_results[:current_i],
        count_results[:current_i],
        candidate_counts[:current_i],
        spectrum_sizes[:current_i],
    )


@numba.njit(cache=True, nogil=True)
def annotate_network(
    queries,
    indptr,
    indices,
    edge_pointers,
    selected_edges,
    low_limits,
    high_limits,
    peptide_pointers,
):
    peptide_count = np.max(peptide_pointers) + 1
    count = len(queries)
    score_results = np.empty(count, np.float64)
    fragment_results = np.empty(count, np.int64)
    index_results = np.empty(count, np.int64)
    count_results = np.empty(count, np.int64)
    candidate_counts = np.empty(count, np.int64)
    neighbor_counts = np.empty(count, np.int64)
    current_i = 0
    for ion_index in queries:
        peptide_low = low_limits[ion_index]
        peptide_high = high_limits[ion_index]
        if peptide_low == peptide_high:
            continue
        ion_start = indptr[ion_index]
        ion_end = indptr[ion_index + 1]
        good_neighbors = selected_edges[edge_pointers[ion_start: ion_end]]
        neighbor_count = np.sum(good_neighbors)
        if neighbor_count == 0:
            continue
        neighbors = indices[ion_start: ion_end][good_neighbors]
        candidates = np.zeros(peptide_count, np.int64)
        for neighbor_ion_index in neighbors:
            neighbor_peptide_low = low_limits[neighbor_ion_index]
            neighbor_peptide_high = high_limits[neighbor_ion_index]
            if neighbor_peptide_low == neighbor_peptide_high:
                continue
            peptides = peptide_pointers[
                neighbor_peptide_low: neighbor_peptide_high
            ]
            candidates[peptides] += 1
        (
            score,
            max_count,
            max_fragment,
            candidate_count
        ) = score_regression_estimator(
            candidates[peptide_pointers[peptide_low: peptide_high]] + 1,
            peptide_low,
            peptide_count
        )
        if score > 0:
            score_results[current_i] = score
            fragment_results[current_i] = max_fragment
            index_results[current_i] = ion_index
            count_results[current_i] = max_count
            candidate_counts[current_i] = candidate_count
            neighbor_counts[current_i] = neighbor_count
            current_i += 1
    return (
        score_results[:current_i],
        fragment_results[:current_i],
        index_results[:current_i],
        count_results[:current_i],
        candidate_counts[:current_i],
        neighbor_counts[:current_i],
    )


@numba.njit(cache=True, nogil=True)
def score_regression_estimator(candidates, offset, peptide_count):
    frequencies = np.bincount(candidates)
    frequencies = np.cumsum(frequencies[::-1])[::-1]
    max_count = len(frequencies) - 1
    max_fragment = offset + np.flatnonzero(candidates == max_count)[0]
    if frequencies[-1] != 1:
        score = 0
    elif frequencies[1] == 1:
        # score = 1 - 2**(-np.log2(peptide_count) * (max_count - 1))
        score = 1 - peptide_count**(1 - max_count)
    else:
        x0 = 2 + np.flatnonzero(frequencies[2:] == 1)[0]
        y0 = np.log2(frequencies[1])
        slope = y0 / (x0 - 1)
        score = 1 - 2**(-slope * (max_count - x0))
    return score, max_count, max_fragment, len(candidates)
