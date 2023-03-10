use std::time::Duration;
use cgmath::Vector3;

mod smoothed_particle_hydrodynamics;

use crate::smoothed_particle_hydrodynamics::simulation_definition::*;
use crate::smoothed_particle_hydrodynamics::simulation_data::*;
use crate::smoothed_particle_hydrodynamics::smoothing_kernels::poly_six_kernel::PolySixKernel;

pub type SimulationFloat = f64;

const TIME_STEP_SECONDS: Duration = Duration::from_micros(1);
const SIM_LENGTH_SECONDS: Duration = Duration::from_secs(10);
const NUM_PARTICLES: u64 = 12000;
// should be 12000 = 1.2 Litres
const SMOOTHING_RADIUS: SimulationFloat = 0.001;
// 1mm
const GRAVITY: Vector3<SimulationFloat> = Vector3::new(0.0, -9.8, 0.0);
const PARTICLE_MASS_KG: SimulationFloat = 0.000001;
//1uL of water
const FLUID_CONST: SimulationFloat = 2200000000.0;
const VISCOUS_CONST: SimulationFloat = 0.0016;
const RHO_ZERO: SimulationFloat = 1000.0;

//TODO - Voxelise

fn main() {
    let simulation_definition: SimulationDefinition<PolySixKernel> = SimulationDefinition::new(
        TIME_STEP_SECONDS,
        SIM_LENGTH_SECONDS,
        NUM_PARTICLES,
        SMOOTHING_RADIUS,
        GRAVITY,
        PARTICLE_MASS_KG,
        FLUID_CONST,
        VISCOUS_CONST,
        RHO_ZERO,
    );

    let mut sim_data: SimulationData<PolySixKernel> = SimulationData::new(simulation_definition);

    while !sim_data.is_finished {
        let start = std::time::Instant::now();
        SimulationData::step(&mut sim_data);
        SimulationData::trace_first_pixel(&sim_data);
        let end = std::time::Instant::now();
        println!("Iteration Time: {:?}", end - start);
    }
}

