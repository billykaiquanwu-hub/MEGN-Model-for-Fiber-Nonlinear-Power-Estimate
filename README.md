# MEGN MATLAB Demo

This repository provides a minimal MATLAB implementation of the finite-memory EGN (MEGN) model for self-channel interference correction.

The main entry point is:

```matlab
MEGN_Demo
```

The demo compares the standard EGN estimate with the MEGN correction for different CCDM shaping blocklengths and mapping dimensions. It produces two figures:

- nonlinear coefficient eta versus CCDM shaping blocklength;
- optimum effective SNR versus CCDM shaping blocklength.

## Repository Contents

- `MEGN_Demo.m` - minimal working example and reader-facing script.
- `MEGNModel.m` - reusable MEGN model class with kernel loading, kernel computation, and SCI correction.
- `calculateEGNSCI.m` - isolated baseline EGN calculation.
- `calculateSymbolCorrelations.m` - analytical CCDM symbol-correlation calculation.
- `etaToSNR.m` - eta-to-optimum-SNR conversion.
- `kernel_cache/` - precomputed kernel files used to avoid repeating the slow kernel calculation.

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

The first run should load kernels from `kernel_cache/`. If kernels are missing for a changed parameter set, `MEGNModel` computes them and stores them for future runs.

## Notes on Kernel Caching

Kernel computation is the main runtime bottleneck. Cached kernels are organized by physical link and model-memory parameters. Keeping the cache in the repository makes the minimum demo fast and reproducible for readers.

If you change parameters such as symbol rate, span count, or memory length, the model may create a new cache entry.

## Troubleshooting

If MATLAB cannot find core functions such as `addpath`, `num2str`, or `fullfile`, reset the MATLAB path:

```matlab
restoredefaultpath
rehash toolboxcache
savepath
matlabrc
```

## Citation

Citation details will be added after the manuscript information is finalized.
