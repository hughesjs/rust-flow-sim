# Rust Flow Sim

This repo contains a work-in-progress smoothed particle hydrodynamics flow simulation written in rust.

It serves the dual purpose of letting me teach myself rust, and SPH simultaneously.

Currently, it's using a simplified model that assumes incompressible flow.

In all likelihood, this will change when I need to consider hypervelocity collisions.

It also needs converting to OpenCL because currently it just runs on your CPU which will cook it.