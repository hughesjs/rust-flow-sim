use std::time::Duration;
use cgmath::Vector4;

mod smoothed_particle_hydrodynamics;
use crate::smoothed_particle_hydrodynamics::simulation_definition::*;
use crate::smoothed_particle_hydrodynamics::simulation_data::*;

const TIME_STEP_SECONDS: Duration = Duration::from_micros(1);
const SIM_LENGTH_SECONDS: Duration = Duration::from_secs(10);
const NUM_PARTICLES: u64 = 12000; // should be 12000 = 1.2 Litres
const SMOOTHING_RADIUS: f64 = 0.001; // 1mm
const GRAVITY: Vector4<f64> = Vector4::new(0.0, -9.8, 0.0, 0.0);
const PARTICLE_MASS_KG: f64 = 0.000001; //1uL of water
const FLUID_CONST: f64 = 2200000000.0;
const VISCOUS_CONST: f64 = 0.0016;
const RHO_ZERO: f64 = 1000.0;

//TODO - Voxelise

fn main() {
    let simulation_definition = SimulationDefinition{
        time_step: TIME_STEP_SECONDS,
        sim_length: SIM_LENGTH_SECONDS,
        num_particles: NUM_PARTICLES,
        smoothing_radius: SMOOTHING_RADIUS,
        gravity: GRAVITY,
        particle_mass_kg: PARTICLE_MASS_KG,
        fluid_constant: FLUID_CONST,
        viscous_constant: VISCOUS_CONST,
        rho_zero: RHO_ZERO,
    };
    let mut sim_data: SimulationData = SimulationData::new(simulation_definition);

    while !sim_data.is_finished {
        let start = std::time::Instant::now();
        SimulationData::step(&mut sim_data);
        SimulationData::trace_first_pixel(&sim_data);
        let end = std::time::Instant::now();
        println!("Iteration Time: {:?}", end - start);
    }
}

