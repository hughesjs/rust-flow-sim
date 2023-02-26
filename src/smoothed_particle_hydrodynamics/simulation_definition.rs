use std::time::Duration;
use cgmath::Vector4;

pub struct SimulationDefinition {
    pub(crate) time_step: Duration,
    pub(crate) sim_length: Duration,
    pub(crate) num_particles: u64,
    pub(crate) smoothing_radius: f64,
    pub(crate) gravity: Vector4<f64>,
    pub(crate) particle_mass_kg: f64,
    pub(crate) fluid_constant: f64,
    pub(crate) viscous_constant: f64,
    pub(crate) rho_zero: f64
}