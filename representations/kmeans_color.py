"""
K-means Color Clustering
========================
Clusters image pixels in LAB color space into K groups. Because H&E (and
other histological) stains produce reproducible, constrained color palettes,
the clusters naturally correspond to tissue compartments:

  K=3  →  background/lumen | nuclear/epithelial | stroma/cytoplasm
  K=5  →  further splits lumen from background and stroma from cytoplasm

No labels or AI models required. Purely unsupervised.
"""

import numpy as np
from skimage.color import rgb2lab
from sklearn.cluster import KMeans, MiniBatchKMeans


# Fixed palette for false-colour label maps (up to 8 clusters)
_PALETTE = np.array([
    [255, 255, 255],   # 0 – white  (typically background)
    [100, 149, 237],   # 1 – cornflower blue
    [220,  20,  60],   # 2 – crimson
    [ 50, 205,  50],   # 3 – lime green
    [255, 165,   0],   # 4 – orange
    [148,   0, 211],   # 5 – dark violet
    [  0, 206, 209],   # 6 – dark turquoise
    [255, 215,   0],   # 7 – gold
], dtype=np.uint8)


class KMeansColor:
    """Segment a histology image by colour using K-means in LAB space."""

    def __init__(self, k_values: tuple = (3, 5), random_state: int = 42,
                 use_minibatch: bool = True):
        """
        Parameters
        ----------
        k_values      : number of clusters to try (one result per K)
        random_state  : for reproducibility
        use_minibatch : use MiniBatchKMeans (faster on large images)
        """
        self.k_values = k_values
        self.random_state = random_state
        self.use_minibatch = use_minibatch

    def compute(self, rgb: np.ndarray) -> dict:
        """
        Parameters
        ----------
        rgb : np.ndarray  shape (H, W, 3)  dtype uint8

        Returns
        -------
        dict  keyed by k value, each containing:
            'label_map'      – (H, W) int array of cluster indices
            'false_colour'   – (H, W, 3) uint8 false-colour image
            'cluster_means'  – (K, 3) mean LAB colour per cluster
            'cluster_sizes'  – (K,) pixel count per cluster
        """
        H, W = rgb.shape[:2]
        lab = rgb2lab(rgb.astype(np.float32) / 255.0)  # (H, W, 3)
        pixels = lab.reshape(-1, 3)                     # (H*W, 3)

        results = {}
        for k in self.k_values:
            Cls = MiniBatchKMeans if self.use_minibatch else KMeans
            km = Cls(
                n_clusters=k,
                random_state=self.random_state,
                n_init=10,
            )
            labels = km.fit_predict(pixels)             # (H*W,)
            label_map = labels.reshape(H, W)

            # Sort clusters by lightness so label 0 is always brightest
            # (background) — makes cross-image comparison easier
            cluster_means = km.cluster_centers_          # (K, 3) in LAB
            order = np.argsort(-cluster_means[:, 0])     # descending L*
            remap = np.empty(k, dtype=int)
            remap[order] = np.arange(k)
            label_map = remap[label_map]
            cluster_means = cluster_means[order]

            # False-colour map
            palette = _PALETTE[:k]
            false_colour = palette[label_map]            # (H, W, 3) uint8

            sizes = np.bincount(label_map.ravel(), minlength=k)

            results[k] = {
                "label_map": label_map,
                "false_colour": false_colour,
                "cluster_means": cluster_means,
                "cluster_sizes": sizes,
            }

        return results
