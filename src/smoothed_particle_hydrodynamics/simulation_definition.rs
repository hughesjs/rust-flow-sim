use std::time::Duration;
use cgmath::Vector4;
use crate::SimulationFloat;
use crate::smoothed_particle_hydrodynamics::smoothing_kernels::smoothing_kernel::SmoothingKernel;

pub struct SimulationDefinition<TKernel: SmoothingKernel> {
    pub(crate) time_step: Duration,
    pub(crate) time_step_in_secs_as_float: SimulationFloat,
    pub(crate) sim_length: Duration,
    pub(crate) num_particles: u64,
    pub(crate) smoothing_radius: SimulationFloat,
    pub(crate) gravity: Vector4<SimulationFloat>,
    pub(crate) particle_mass_kg: SimulationFloat,
    pub(crate) fluid_constant: SimulationFloat,
    pub(crate) viscous_constant: SimulationFloat,
    pub(crate) rho_zero: SimulationFloat,
    pub(crate) smoothing_kernel: TKernel,
}

impl<TKernel: SmoothingKernel> SimulationDefinition<TKernel> {
    pub fn new(time_step: Duration,
               sim_length: Duration,
               num_particles: u64, smoothing_radius: SimulationFloat, gravity: Vector4<SimulationFloat>,
               particle_mass_kg: SimulationFloat,
               fluid_constant: SimulationFloat,
               viscous_constant: SimulationFloat,
               rho_zero: SimulationFloat
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
            time_step_in_secs_as_float: time_step.as_secs_f64() as SimulationFloat,
            smoothing_kernel: TKernel::new(smoothing_radius),
        }
    }
}