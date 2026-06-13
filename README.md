# MEGN MATLAB Demo

This repository provides a minimal MATLAB implementation of the finite-memory EGN (MEGN) model for fiber-optic channels.

The main entry point is:

```matlab
MEGN_Demo
```

The demo compares the standard EGN estimate with the MEGN correction for different CCDM shaping blocklengths and mapping strategies. It reproduces the calculations behind Fig. 9(b) of the accompanying arXiv paper and produces two figures:

- nonlinear coefficient eta versus CCDM shaping blocklength;
- optimum effective SNR versus CCDM shaping blocklength.

## Repository Contents

- `MEGN_Demo.m` - minimal working example and reader-facing script.
- `MEGNModel.m` - reusable MEGN model class with kernel loading, kernel computation, and SCI correction.
- `calculateEGNSCI.m` - isolated baseline EGN calculation.
- `calculateSymbolCorrelations.m` - analytical CCDM symbol-correlation calculation.
- `etaToSNR.m` - eta-to-optimum-SNR conversion.
- `kernel_cache/` - precomputed kernel file for the demo configuration.

## Requirements

- MATLAB R2022a or newer is recommended.
- No external MATLAB toolboxes are required for the current demo.

## Quick Start

1. Open MATLAB.
2. Set the current folder to this repository.
3. Run:

```matlab
MEGN_Demo
```

The first run should load the 10-span demo kernel from `kernel_cache/`. If kernels are missing for a changed parameter set, `MEGNModel` computes them and stores them for future runs.

## Notes on Kernel Caching

Kernel computation is the main runtime bottleneck. The repository includes only the kernel file required by the minimum working example, where `sys.Ns = 10` and `model.tau_max = 50`. Cached kernels are organized by physical link and model-memory parameters.

If you change parameters such as symbol rate, span count, or memory length, the model may create a new cache entry. Longer-distance links generally require a larger memory length (`model.tau_max`) to capture more delayed symbol correlations, which also increases kernel-computation time and cache size.

## Troubleshooting

If MATLAB cannot find core functions such as `addpath`, `num2str`, or `fullfile`, reset the MATLAB path:

```matlab
restoredefaultpath
rehash toolboxcache
savepath
matlabrc
```

## Citation

If you use this code, please cite:

Kaiquan Wu, Gabriele Liga, Marco Secondini, Stella Civelli, Hussam Batshon, Greg Raybon, Xi Chen, and Alex Alvarado, "The Memory-Enhanced Gaussian Noise (MEGN) Model for Fiber-Optic Channels," arXiv:2604.11601, 2026. https://arxiv.org/abs/2604.11601

```bibtex
@misc{wu2026megn,
  title = {The Memory-Enhanced Gaussian Noise (MEGN) Model for Fiber-Optic Channels},
  author = {Wu, Kaiquan and Liga, Gabriele and Secondini, Marco and Civelli, Stella and Batshon, Hussam and Raybon, Greg and Chen, Xi and Alvarado, Alex},
  year = {2026},
  eprint = {2604.11601},
  archivePrefix = {arXiv},
  primaryClass = {eess.SP},
  doi = {10.48550/arXiv.2604.11601},
  url = {https://arxiv.org/abs/2604.11601}
}
```
