# Portfolio copy

## Short card

**Computational Geometry CAD Filter**

C++ · PCL · Eigen · NumPy · React

I extended a shrinking-ball medial-axis pipeline to expose geometric scale features from sampled CAD surfaces. The revitalised project includes a modern native build and an interactive browser reconstruction that makes the underlying geometry explorable.

[Interactive demo](https://computational-geometry-cad-filter.flabdul.chatgpt.site) · [Source](https://github.com/FLABDUL/beng-thesis)

## Case-study summary

For my BEng thesis, *Development of a Machine-Learning Computer-Aided Design Filter Using Computational Geometry*, I investigated whether shape-aware geometric measurements could support automated CAD filtering. Starting from TU Delft's `masbcpp` research implementation, I adapted the C++ pipeline to process mesh samples exported as NumPy arrays and added medial-ball radii as an explicit signal for downstream analysis.

The key idea is local feature size: a ball shrinks along a surface normal until it touches the sampled shape in at least two places. Repeating that process across the surface approximates its medial axis. Small medial balls identify narrow or detailed regions; larger balls describe broader mass. That gives a model-aware alternative to simplifying every part of a shape uniformly.

### My contribution

- Adapted the geometry I/O to the `trimesh`/NumPy research workflow.
- Added interior and exterior medial-ball radius outputs.
- Tightened convergence and increased the iteration budget for detailed CAD samples.
- Built notebook analysis around the resulting radius distributions.
- Restored the project with reproducible builds, safer I/O, CI checks and an interactive portfolio demo.

### Honest scope

The surviving repository represents the computational-geometry feature-extraction stage, not a complete trained ML product. The browser visualisation is a faithful educational reconstruction of the shrinking-ball concept, while the C++/PCL pipeline is the research implementation.
