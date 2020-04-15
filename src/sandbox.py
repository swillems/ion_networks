#!python

import numpy as np
import numba
import multiprocessing as mp


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


def matrix_multiplication(a, b):
    aa_indices, aa_values = a.nonzero()
    bb_values, bb_indices = b.nonzero()
    return matrix_multiplication_numba_wrapper(
        aa_indices,
        bb_indices,
        aa_values,
        bb_values,
        np.argsort(aa_values),
        np.argsort(bb_values),
    )


@numba.njit(
    numba.i4[:, :](
        numba.i4[:],
        numba.i4[:],
        numba.i4[:],
        numba.i4[:],
        numba.i8[:],
        numba.i8[:]
    )
)
def matrix_multiplication_numba_wrapper(
    aa_indices,
    bb_indices,
    aa_values,
    bb_values,
    aa_order,
    bb_order
):
    aa_indices = aa_indices[aa_order]
    aa_values = aa_values[aa_order]
    bb_indices = bb_indices[bb_order]
    bb_values = bb_values[bb_order]
    low_limits = np.searchsorted(
        aa_values,
        bb_values,
        "left"
    )
    high_limits = np.searchsorted(
        aa_values,
        bb_values,
        "right"
    )
    diffs = high_limits - low_limits
    ends = np.cumsum(diffs)
    result = np.empty((2, ends[-1]), np.int32)
    result[1, :] = np.repeat(
        bb_indices,
        high_limits - low_limits
    )
    for l, h, e, d in zip(low_limits, high_limits, ends, diffs):
        result[0, e - d: e] = aa_indices[l: h]
    return result


def parallelizedGenerator(
    process_count,
    function,
    iterable,
    *args,
    **kwargs
):
    in_queue = mp.Queue()
    out_queue = mp.Queue()

    def _parallel_function():
        while True:
            element = in_queue.get()
            if element is None:
                out_queue.put(None)
                return
            result = function(
                element,
                *args,
                **kwargs
            )
            out_queue.put(result)

    for element in iterable:
        in_queue.put(element)
    processes = []
    for process_index in range(process_count):
        process = mp.Process(
            target=_parallel_function,
        )
        in_queue.put(None)
        process.start()
        processes.append(process)
    running_processes = process_count
    while running_processes != 0:
        result = out_queue.get()
        if result is None:
            running_processes -= 1
        else:
            yield result
    in_queue.close()
    out_queue.close()
    in_queue.join_thread()
    out_queue.join_thread()
    while processes:
        process = processes.pop()
        process.join()


import functools
import time
import multiprocessing as mp


def generate_in_parallel(cpu_count=0):
    max_cpu_count = mp.cpu_count()
    if cpu_count > max_cpu_count:
        cpu_count = max_cpu_count
    elif cpu_count <= 0:
        if cpu_count + max_cpu_count <= 0:
            cpu_count = 1
        else:
            cpu_count += max_cpu_count
    def decorator_repeat(func):
        @functools.wraps(func)
        def wrapper_generate_in_parallel(*args, **kwargs):
            in_queue = mp.Queue()
            out_queue = mp.Queue()
            def _parallel_function():
                while True:
                    element = in_queue.get()
                    if element is None:
                        out_queue.put(None)
                        return
                    result = func(
                        element,
                        *args[1:],
                        **kwargs
                    )
                    out_queue.put(result)
            for element in args[0]:
                in_queue.put(element)
            processes = []
            for process_index in range(cpu_count):
                process = mp.Process(
                    target=_parallel_function,
                )
                in_queue.put(None)
                process.start()
                processes.append(process)
            running_processes = cpu_count
            while running_processes != 0:
                result = out_queue.get()
                if result is None:
                    running_processes -= 1
                else:
                    yield result
            in_queue.close()
            out_queue.close()
            in_queue.join_thread()
            out_queue.join_thread()
            while processes:
                process = processes.pop()
                process.join()
        return wrapper_generate_in_parallel
    return decorator_repeat
