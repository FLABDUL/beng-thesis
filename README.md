# Computational Geometry CAD Filter

An interactive restoration of my BEng thesis project, **Development of a Machine-Learning Computer-Aided Design Filter Using Computational Geometry**.

The project investigates how the medial axis and local feature size of a sampled CAD surface can produce geometric signals for downstream filtering. It extends the TU Delft [`masbcpp`](https://github.com/tudelft3d/masbcpp) shrinking-ball implementation with explicit medial-ball radius output, tighter convergence controls and a mesh-to-NumPy workflow.

> The original repository contains the computational-geometry stage of the thesis workflow. It does not contain a trained machine-learning model. The browser experience is an educational 2D reconstruction of the algorithm; the native C++ tools remain the reference implementation.

## Try the interactive demo

The portfolio demo lives in [`demo/`](demo). It turns the core shrinking-ball idea into an explorable visual: choose a CAD-like profile, tune the initial radius and feature threshold, then inspect the retained medial centres and radius distribution.

**[Open the private live demo](https://computational-geometry-cad-filter.flabdul.chatgpt.site)**

```bash
cd demo
npm ci
npm run dev
```

The site uses Node.js 22 or newer. Run `npm test` to build it and execute its server-rendering checks.

## How the native pipeline works

1. Sample a CAD mesh into oriented surface points.
2. Load `coords.npy` and `normals.npy` as `N x 3` arrays.
3. Roll shrinking balls along both normal directions to approximate interior and exterior medial-axis points.
4. Record each ball centre, its opposing surface-point index and its radius.
5. Use medial geometry to estimate local feature size and simplify the point cloud while preserving small features.

The radius arrays added in this thesis fork make the geometric scale signal directly available to notebooks and downstream filtering experiments.

## Native build

### Ubuntu / Debian

```bash
sudo apt-get install build-essential cmake ninja-build libeigen3-dev libpcl-dev zlib1g-dev
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

This produces three command-line tools:

| Tool | Purpose |
| --- | --- |
| `compute_normals` | Estimate surface normals with PCA when only coordinates are available. |
| `compute_ma` | Approximate interior and exterior medial balls. |
| `simplify` | Compute local feature size and feature-aware point-cloud simplification. |

Run any tool with `--help` for its complete options.

### Docker

The container uses Ubuntu 24.04 and distribution-provided PCL rather than compiling a historical PCL release from source.

```bash
docker compose build
docker compose run --rm geometry-filter compute_ma /data/input /data/output
```

Place inputs under `data/input`; results will be written to `data/output` on the host.

## Input and output

`compute_ma` accepts NumPy `float32` or `float64` arrays:

```text
input/
├── coords.npy      # N x 3 point coordinates
└── normals.npy     # N x 3 consistently oriented unit normals
```

Example:

```bash
./build/compute_ma --radius 200 --convergence 1e-7 --iterations 200 input output
```

It writes:

```text
output/
├── ma_coords_in.npy   # interior medial-ball centres, N x 3 float32
├── ma_coords_out.npy  # exterior medial-ball centres, N x 3 float32
├── ma_qidx_in.npy     # opposing input point indices, N int32
├── ma_qidx_out.npy
├── ma_rad_in.npy      # interior ball radii, N float32
├── ma_rad_out.npy     # exterior ball radii, N float32
└── compute_ma.txt     # parameters used for the run
```

The historical research workflow is preserved in [`notebooks/debug_masbcpp.ipynb`](notebooks/debug_masbcpp.ipynb). It uses `trimesh` to sample an STL, exports the NumPy arrays and inspects the resulting radii.

## What was revitalised

- Reframed the repository around the thesis contribution and documented its provenance honestly.
- Added a responsive, accessible browser demo suitable for linking from a portfolio.
- Replaced the 2016-era CMake and Ubuntu image with target-based CMake, bundled `cnpy`, modern dependencies and CI smoke tests.
- Cleaned generated binaries, IDE state and notebook checkpoints from version control.
- Made NumPy loading accept both `float32` and `float64` geometry.
- Made radius output consistent with the final medial-ball centre and exposed convergence/iteration limits as CLI options.

For a ready-to-paste portfolio summary, see [`PORTFOLIO.md`](PORTFOLIO.md).

## Limitations

- Runtime increases substantially for dense meshes; pre-sampling or subdivision choices matter.
- Input and output use NumPy arrays rather than reading CAD formats directly.
- Results depend on point density, normal orientation, initial radius and denoising thresholds.
- The interactive demo is deliberately 2D and illustrative, not a WebAssembly port of PCL.

## Credits and licence

The native implementation is derived from Ravi Peters and the [TU Delft 3D geoinformation group](https://github.com/tudelft3d/masbcpp), with an intermediate fork by [Drew Sherlock](https://github.com/drewsherlock/masbcpp). My thesis fork adds the research workflow and medial-ball radius signal described above.

Released under the [MIT License](LICENSE). Original copyright notices are retained in source files.
