#!/usr/bin/python

import numpy as np
#import sklearn.datasets as data
import re
import sys
import hdbscan

if __name__ == '__main__':
    input_file = sys.argv[1]
    cluster_size = eval(sys.argv[2])
    output_file = sys.argv[3]

    f = open(input_file, 'r')
    out1 = open(output_file, 'w')
    t1 = []
    t2 = []
    for line in f:
        line = line.strip('\n')
        sent = line.split('\t')
        t1.append([eval(sent[0]),eval(sent[2])])

    test_data1 = np.array(t1)
    #test_data1 = np.vstack(t1)
    #test_data = np.vstack([moons, blobs])
    f.close()
    #moons, _ = data.make_moons(n_samples=50, noise=0.05)
    #blobs, _ = data.make_blobs(n_samples=50, centers=[(-0.75,2.25), (1.0, 2.0)], cluster_std=0.25)
    #test_data = np.vstack([moons, blobs])

    #print test_data

    clusterer = hdbscan.HDBSCAN(min_cluster_size=cluster_size, prediction_data=True, gen_min_span_tree=True)
    cluster_labels = clusterer.fit_predict(test_data1)

    #print len(cluster_labels),

    for i in range(len(cluster_labels)):
        out1.write("%d\t%d\t%d\n" % (test_data1[i,0],test_data1[i,1],cluster_labels[i]+1))

    out1.close()
