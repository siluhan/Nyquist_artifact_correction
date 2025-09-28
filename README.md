# **Coil-signature-based Reconstruction for Multi-shot EPI Nyquist Aritfact Correction**

## 📘 **Description**
This repository contains the implementation of the method described in our manuscript: "Correction of Aliasing Artifact in Accelerated Echo-planar Imaging". Submitted to *Magnetic Resonance Imaging*.

The code is written in MATLAB R2024a and is intended to reproduce the experimental results presented in the paper.

## ⚠️ **Important Notes**
- Unoptimized Code:

  This code is a direct and unoptimized implementation of our method.
  As a result, runtime may be relatively long, especially for large datasets. Nevertheless, the results are consistent with those reported in the manuscript.

- Exclusion of MUSE Reconstruction Code:
  
  Please note that MUSE reconstruction code is not included in this repository.
  The reconstruction was performed externally using a separate pipeline.
  The provided code assumes that the MUSE-reconstructed data is already available as input.
