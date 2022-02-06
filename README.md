# Cluster Center Consistency Guided Sampling Learning for Multiple Kernel Clustering
---

authors: Jiali You, Yanglei Hou, Zhenwen Ren, Xiaojian You, Quansen Sun, Xiaobo Shen, Yuancheng Yao
year: 2021

---

## abstract
Multiple kernel clustering (MKC) is popular in processing non-linear data over the past years, where the main challenge is that the kernel matrix with the size $n\times n$ leads to high memory and CPU consumption ($n$ denotes the number of samples). Although the widely used preprocessed anchor sampling can mitigate such a challenging considerably, how to dynamically select anchor points with a learning manner from kernel matrix is a difficult problem. To address these issues, this paper proposes a novel method dubbed as cluster center consistency guided sampling learning (3CSL) for multiple kernel clustering (3CSL-MKC). Specifically, by taking the cluster center consistency between the original partial kernel data and anchor points into consideration, 3CSL-MKC learns the shared anchor sampling matrix gradually. With the help of high-quality anchors, the essential clustering information of each kernel partition can be transformed largely into a concentrated low-dimensional representation matrix. Meanwhile, based on dictionary learning, 3CSL-MKC fuses these candidate representation matrixes to produce the resulting consensus representation, such that the clustering assignments can be directly obtained relying on simple $k$-means. A large number of experiments are conducted on different multiple kernel datasets to verify the effectiveness and efficiency of the proposed method.
