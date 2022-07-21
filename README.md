# Primitive-based Shape Abstraction via Nonparametric Bayesian Inference

This repo provides the source code for the ECCV 2022 paper:

> [**Primitive-based Shape Abstraction via Nonparametric Bayesian Inference**](https://arxiv.org/pdf/2203.14714.pdf "ArXiv version of the paper.")  
> Yuwei Wu, Weixiao Liu, [Sipu Ruan](https://ruansp.github.io/), [Gregory S. Chirikjian](https://cde.nus.edu.sg/me/staff/chirikjian-gregory-s/)

We propose an algorithm to infer a superquadric-based abstraction from a point cloud with good accuracy and concise primitives.  

<img src="/figures/SQ.PNG" alt="superquadrics" width="600" align="center"/>

## Abstration

3D shape abstraction has drawn great interest over the years. Apart from low-level representations such as meshes and voxels, researchers also seek to semantically abstract complex objects with basic geometric primitives. 
Recent deep learning methods rely heavily on datasets, with limited generality to unseen categories.
Furthermore, abstracting an object accurately yet with a small number of primitives still remains a challenge.
In this paper, we propose a novel non-parametric Bayesian statistical method to infer an abstraction, consisting of an unknown number of geometric primitives, from a point cloud.
We model the generation of points as observations sampled from an infinite mixture of Gaussian Superquadric Taper Models (GSTM).
Our approach formulates the abstraction as a clustering problem, in which: 1) each point is assigned to a cluster via the Chinese Restaurant Process (CRP); 2) a primitive representation is optimized for each cluster, and 3) a merging post-process is incorporated to provide a concise representation.
We conduct extensive experiments on two datasets.
The results indicate that our method outperforms the state-of-the-art in terms of accuracy and is generalizable to various types of objects.

## Implementations

This repo provides the MATLAB implementation.

### Installation

The code is written and tested on MATALB R2022a and
should be compatible with releases newer than R2018a.

### File Structure

- `superquadricSegment.m` is the main code for implementing our shape abstraction algorithm.
- Add `/algorithm` and `/utility_functions` to MATLAB path before running the code

### Run Demo

- Demo `.ply` files are located inside `/example`
- To run demo example, run `test_script.m`