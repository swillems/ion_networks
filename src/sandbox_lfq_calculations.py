#!python

import numpy as np


# Calculating edge accuracy
with log.newSection("Calculating consistent coeluton accuracy"):
    classes = np.round(logfcs)
    a_indices, b_indices = neighbors.nonzero()
    overlap = neighbors.data
    a_classes = classes[a_indices]
    b_classes = classes[b_indices]
    c = np.in1d(a_classes, [-2, 0, 1])
    c &= np.in1d(b_classes, [-2, 0, 1])
    a_classes = a_classes[c]
    b_classes = b_classes[c]
    a_indices = a_indices[c]
    b_indices = b_indices[c]
    overlap = neighbors.data[c]
    eco1 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    yea1 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    hum1 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    eco2 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    yea2 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    hum2 = np.zeros(parameters["SAMPLE_COUNT"] + 1)
    for i in range(2, parameters["SAMPLE_COUNT"] + 1):
        full_overlap = overlap == i
        a_full = a_classes[full_overlap]
        b_full = b_classes[full_overlap]
        ecoli, human, yeast = np.unique(np.concatenate([a_full, b_full]), return_counts=True)[1]
        total = ecoli + human + yeast
        ecoli_ratio = ecoli / total
        human_ratio = human / total
        yeast_ratio = yeast / total
        ecoli_expected = ecoli_ratio * ecoli_ratio
        human_expected = human_ratio * human_ratio
        yeast_expected = yeast_ratio * yeast_ratio
        equal = a_full == b_full
        ecoli_hits = equal[(a_full == -2) | (a_full == -2)]
        human_hits = equal[(a_full == 0) | (a_full == 0)]
        yeast_hits = equal[(a_full == 1) | (a_full == 1)]
        ecoli_hit_ratio = 2 * np.sum(ecoli_hits) / total
        yeast_hit_ratio = 2 * np.sum(yeast_hits) / total
        human_hit_ratio = 2 * np.sum(human_hits) / total
        eco1[i] = ecoli_hit_ratio / ecoli_ratio
        yea1[i] = yeast_hit_ratio / yeast_ratio
        hum1[i] = human_hit_ratio / human_ratio
        eco2[i] = ecoli_expected / ecoli_ratio
        yea2[i] = yeast_expected / yeast_ratio
        hum2[i] = human_expected / human_ratio

# Plotting edge accuracy
with log.newSection("Plotting consitent coelution accuracy"):
    e_ex, = plt.plot(100 * eco1[2:], marker="o", linestyle="-", c="blue")
    h_ex, = plt.plot(100 * hum1[2:], marker="o", linestyle="-", c="red")
    y_ex, = plt.plot(100 * yea1[2:], marker="o", linestyle="-", c="green")
    e_th, = plt.plot(100 * eco2[2:], marker=".", linestyle=":", c="blue")
    h_th, = plt.plot(100 * hum2[2:], marker=".", linestyle=":", c="red")
    y_th, = plt.plot(100 * yea2[2:], marker=".", linestyle=":", c="green")
    tmp = plt.legend(
        (
            e_ex,
            h_ex,
            y_ex,
            e_th,
            h_th,
            y_th
        ), (
            "LogFC = -2 Experimental (ECOLI)",
            "LogFC = 0 Experimental (HUMAN)",
            "LogFC = 1 Experimental (YEAST)",
            "LogFC = -2 Expected (ECOLI)",
            "LogFC = 0 Expected (HUMAN)",
            "LogFC = 1 Expected (YEAST)"
        ),
        loc="upper left"
    )
    tmp = plt.xlabel("Consistently Co-eluting Runs")
    tmp = plt.xticks(
        np.arange(parameters["SAMPLE_COUNT"] - 1),
        np.arange(2, parameters["SAMPLE_COUNT"] + 1)
    )
    tmp = plt.ylabel("% Consistent Co-eluting Aggregates With Equal LogFC")
    # tmp = plt.show()
    tmp = plt.savefig(parameters["PLOTS_PATH"] + extension + "_lfq_edge_accuracy.pdf", bbox_inches='tight')
    tmp = plt.close()
