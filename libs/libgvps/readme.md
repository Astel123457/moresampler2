libgvps
---

libgvps stands for *Generic Viterbi Path Searcher*.

This library is an implementation of Viterbi algorithm in a *probabilistic* formulation that helps you **find the most likely sequence given a set of states or candidates over a certain time period and the state transition probabilities**.

libgvps allows the states/candidates space to be dense or sparse or even circular; transition probabilities can be either static or dynamic (along time axis).

| Function | Geometric description | Examples (if any)
| --- | --- | --- |
| `gvps_sparse_sampled` | path searching on a plane
| `gvps_sparse_sampled_hidden` | path searching on a double-sided plane with only one side exposed
| `gvps_sparse_sampled_static` | path searching on a plane assuming time locality of trans. prob | F0 tracking
| `gvps_sparse_sampled_hidden_static` | path searching on a double-sided plane with only one side exposed, assuming time locality of trans. prob | F0 tracking, (lowpass-like) filtering on a discrete/quantized signal
| `gvps_sparse_circular` | path searching on a cylinder
| `gvps_sparse_circular_hidden` | path searching on a double-sided cylinder with only one side exposed
| `gvps_sparse_circular_static` | path searching on a cylinder assuming time locality of trans. prob
| `gvps_sparse_circular_hidden_static` | path searching on a double-sided cylinder with only one side exposed, assuming time locality of trans. prob
| `gvps_full_sampled` | path searching on a plane
| `gvps_full_sampled_hidden` | path searching on a double-sided plane with only one side exposed
| `gvps_full_sampled_static` | path searching on a plane assuming time locality of trans. prob | trajectory tracking in a gray-scale image, (lowpass-like) filtering on a discrete/quantized signal
| `gvps_full_sampled_hidden_static` | path searching on a double-sided plane with only one side exposed, assuming time locality of trans. prob | trajectory tracking in a gray-scale image, allowing discontinuities
| `gvps_full_circular` | path searching on a cylinder
| `gvps_full_circular_hidden` | path searching on a double-sided cylinder with only one side exposed
| `gvps_full_circular_static` | path searching on a cylinder assuming time locality of trans. prob
| `gvps_full_circular_hidden_static` | path searching on a double-sided cylinder with only one side exposed, assuming time locality of trans. prob
| `gvps_full` | path searching on a (topological) plane with arbitrary geometric assumptions
| `gvps_full_static` | path searching on a (topological) plane with arbitrary geometric assumptions, assuming time locality of trans. prob | MAP inference for Hidden Markov Models
| `gvps_variable` | path searching in an arbitrary space | unit selection speech synthesis