use std::time::Duration;
use cgmath::Vector4;
use crate::smoothed_particle_hydrodynamics::smoothing_kernels::smoothing_kernel::SmoothingKernel;

pub struct SimulationDefinition<TKernel: SmoothingKernel> {
    pub(crate) time_step: Duration,
    pub(crate) sim_length: Duration,
    pub(crate) num_particles: u64,
    pub(crate) smoothing_radius: f64,
    pub(crate) gravity: Vector4<f64>,
    pub(crate) particle_mass_kg: f64,
    pub(crate) fluid_constant: f64,
    pub(crate) viscous_constant: f64,
    pub(crate) rho_zero: f64,
    pub(crate) smoothing_kernel: TKernel,
}

impl<TKernel: SmoothingKernel> SimulationDefinition<TKernel> {
    pub fn new(time_step: Duration,
               sim_length: Duration,
               num_particles: u64, smoothing_radius: f64, gravity: Vector4<f64>,
               particle_mass_kg: f64,
               fluid_constant: f64,
               viscous_constant: f64,
               rho_zero: f64
    ) -> Self {
        SimulationDefinition {
            time_step,
            sim_length,
            num_particles,
            smoothing_radius,
            gravity,
            particle_mass_kg,
            fluid_constant,
            viscous_constant,
            rho_zero,
            smoothing_kernel: TKernel::new(smoothing_radius),
        }
    }
}