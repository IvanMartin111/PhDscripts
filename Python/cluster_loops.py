from sklearn.datasets import make_blobs
import pandas as pd
import hdbscan

blobs, labels = make_blobs(n_samples=2000, n_features=10)
pd.DataFrame(blobs).head()

clusterer = hdbscan.HDBSCAN()
clusterer.fit(blobs)

HDBSCAN(algorithm='best', alpha=1.0, approx_min_span_tree=True,
    gen_min_span_tree=False, leaf_size=40, memory=Memory(cachedir=None),
    metric='euclidean', min_cluster_size=5, min_samples=None, p=None)