//use std::collections::HashMap;
use cgmath::{InnerSpace, Vector4, Zero};
use rand::Rng;
use std::f64::consts::PI;
use ndarray::{Array1, Zip};


const NUM_PARTICLES: u32 = 1000; // should be 12000 = 1.2 Litres
const BUFFER_LEN: usize = NUM_PARTICLES as usize;
const SMOOTHING_RADIUS: f64 = 0.001; // 1mm
const TIME_STEP_SECONDS: f64 = 0.01;
const GRAVITY: Vector4<f64> = Vector4::new(0.0, -9.8, 0.0, 0.0);
const SIM_LENGTH_SECONDS: f64 = 10.0;
const PARTICLE_MASS_KG: f64 = 0.000001; //1uL of water
const FLUID_CONST: f64 = 2200000000.0;
const VISCOUS_CONST: f64 = 0.0016;
const RHO_ZERO: f64 = 1000.0;

//TODO - Voxelise

fn main() {
    let mut sim_data: SimulationData = SimulationData::new(SIM_LENGTH_SECONDS, TIME_STEP_SECONDS);
    while !sim_data.is_finished {
        let start = std::time::Instant::now();
        SimulationData::step(&mut sim_data);
        SimulationData::trace_first_pixel(&sim_data);
        let end = std::time::Instant::now();
        println!("Iteration Time: {:?}", end - start);
    }
}

struct SimulationSpace {
    positions: Array1<Vector4<f64>>,
    velocities:  Array1<Vector4<f64>>,
    accelerations:  Array1<Vector4<f64>>,
}

struct SimulationData {
    simulation_space: SimulationSpace,
    current_time_seconds: f64,
    time_step_seconds: f64,
    sim_length_seconds: f64,
    is_finished: bool
    //voxel_pixel_map: HashMap<u16, Vec<u32>>
}

impl SimulationData {
    fn new(sim_length_seconds: f64, time_step_seconds: f64) -> SimulationData {
        let sim_space = SimulationSpace {
            positions: Self::get_initial_positions(),
            velocities: Self::get_initial_velocities(),
            accelerations: Self::get_initial_accelerations(),
        };

        let sim_data = SimulationData {
            simulation_space: sim_space,
            current_time_seconds: 0.0,
            time_step_seconds,
            sim_length_seconds,
            is_finished: false,

            //voxel_pixel_map: HashMap::new()
        };

        return sim_data;
    }

    fn step(sim_data: &mut SimulationData) {
        let densities = Self::calculate_densities(&sim_data.simulation_space);
        let pressures = Self::calculate_pressures(&densities);
        let pressure_grad_terms = Self::calculate_pressure_grad_terms(&sim_data.simulation_space, &pressures, &densities);
        let viscosity_terms = Self::calculate_viscosity_terms(&sim_data.simulation_space, &densities);
        let accelerations = Self::calculate_accelerations(&pressure_grad_terms, &viscosity_terms);
        let new_velocities = Self::calculate_velocities(&sim_data.simulation_space, &accelerations);
        let new_positions = Self::calculate_positions(&sim_data.simulation_space, &new_velocities);

        //TODO - Check for collisions with walls and floor

        sim_data.simulation_space.velocities = new_velocities;
        sim_data.simulation_space.positions = new_positions;
        sim_data.simulation_space.accelerations = accelerations;

        sim_data.current_time_seconds += sim_data.time_step_seconds;
        sim_data.is_finished = sim_data.current_time_seconds > sim_data.sim_length_seconds;
    }

    fn calculate_positions(simulation_space: &SimulationSpace, velocities: &Array1<Vector4<f64>>) -> Array1<Vector4<f64>> {
        Zip::from(&simulation_space.positions)
            .and(velocities)
            .map_collect(|&current_position, &current_velocity| current_position + current_velocity * TIME_STEP_SECONDS)
    }


    fn calculate_velocities(simulation_space: &SimulationSpace, accelerations: &Array1<Vector4<f64>>) -> Array1<Vector4<f64>> {
        Zip::from(&simulation_space.velocities)
            .and(accelerations)
            .map_collect(|&current_velocity, &current_acceleration| current_velocity + current_acceleration * TIME_STEP_SECONDS)
    }

    fn calculate_accelerations(pressure_terms: &Array1<Vector4<f64>>, viscosity_terms: &Array1<Vector4<f64>>) -> Array1<Vector4<f64>> {
        Zip::from(pressure_terms)
            .and(viscosity_terms)
            .map_collect(|&pressure_term, &viscosity_term| GRAVITY + pressure_term + viscosity_term)
    }

    fn calculate_viscosity_terms(simulation_space: &SimulationSpace, densities: &Array1<f64>) -> Array1<Vector4<f64>> {
        Zip::from(&simulation_space.positions)
            .and(densities)
            .and(&simulation_space.velocities)
            .map_collect(|&current_position, &current_density, &current_velocity | {
                Zip::from(&simulation_space.positions)
                    .and(densities)
                    .and(&simulation_space.velocities)
                    .fold(Vector4::zero(), | acc, &other_position, &other_density, &other_velocity | {
                        if Self::is_in_interaction_radius_and_not_self(current_position, other_position) {
                            acc + ((VISCOUS_CONST / current_density) * PARTICLE_MASS_KG * ((other_velocity - current_velocity) / other_density) * Self::laplacian_smooth(current_position, other_position))
                        } else {
                            acc
                        }
                    })
            })
    }

    fn calculate_pressure_grad_terms(simulation_space: &SimulationSpace, densities: &Array1<f64>, pressures: &Array1<f64>) -> Array1<Vector4<f64>> {
        Zip::from(&simulation_space.positions)
            .and(densities)
            .and(pressures)
            .map_collect(|&current_position, &current_density, &current_pressure | {
                Zip::from(&simulation_space.positions)
                    .and(densities)
                    .and(pressures)
                    .fold(Vector4::zero(), | acc, &other_position, &other_density, &other_pressure | {
                        if Self::is_in_interaction_radius_and_not_self(current_position, other_position) {
                            acc + PARTICLE_MASS_KG * ((current_pressure / current_density.powi(2)) + (other_pressure / other_density.powi(2)))* Self::grad_smooth(current_position, other_position)
                        } else {
                            acc
                        }
                    })
            })
    }

    fn calculate_pressures(densities: &Array1<f64>) -> Array1<f64> {
        densities.map(|&density| FLUID_CONST * (density - RHO_ZERO))
    }

    fn calculate_densities(sim_space: &SimulationSpace) -> Array1<f64> {
        sim_space.positions.map(|&current_position| {
            sim_space.positions.fold(0.0, |acc, &other_position| {
                if Self::is_in_interaction_radius_and_not_self(current_position, other_position) {
                    acc + PARTICLE_MASS_KG * Self::smooth(current_position, other_position)
                } else {
                    acc
                }
            })
        })
    }
    fn is_in_interaction_radius_and_not_self(current: Vector4<f64>, other: Vector4<f64>) -> bool {
        current != other && (current - other).magnitude() < SMOOTHING_RADIUS
    }

    // There are some terms reused between smooth, grad_smooth and laplacian_smooth, this could be optimised
    // Not to mention, some of these terms are constants...
    fn smooth(current_position: Vector4<f64>, other_position: Vector4<f64>) -> f64 {
        (315.0 / (64.0 * PI * (SMOOTHING_RADIUS.powi(9)))) * (SMOOTHING_RADIUS.powi(2) - (current_position - other_position).magnitude2()).powi(3)
    }

    fn grad_smooth(current_position: Vector4<f64>, other_position: Vector4<f64>) -> Vector4<f64> {
        (-45.0 / (PI * (SMOOTHING_RADIUS.powi(6))))
            * (SMOOTHING_RADIUS - (current_position - other_position).magnitude2()).powi(2)
            * ((current_position - other_position) / (current_position - other_position).magnitude2())
    }

    fn laplacian_smooth(current_position: Vector4<f64>, other_position: Vector4<f64>) -> f64 {
        (45.0 / (PI * (SMOOTHING_RADIUS.powi(6))))
            * (SMOOTHING_RADIUS - (current_position - other_position).magnitude())
    }

    fn get_initial_positions() -> Array1<Vector4<f64>>  {
        let mut rng = rand::thread_rng();
        let mut positions: Array1<Vector4<f64>> = Array1::from(vec![Vector4::zero(); BUFFER_LEN]);
        for i in 0..BUFFER_LEN {
            positions[i] = Vector4::new(
                rng.gen_range(0.0..0.001),
                rng.gen_range(0.0..0.001),
                rng.gen_range(4.999..5.0),
                0.0,
            );
        }
        return positions;
    }

    fn get_initial_velocities() ->  Array1<Vector4<f64>> {
        Array1::from(vec![Vector4::zero(); BUFFER_LEN])
    }

    fn get_initial_accelerations() ->  Array1<Vector4<f64>> {
        Array1::from(vec![Vector4::zero(); BUFFER_LEN])
    }

    fn trace_first_pixel(sim_data: &SimulationData) {
        println!("{:?}", sim_data.simulation_space.positions[0])
    }
}






// struct Particle {
//     position: Vector3<f64>,
//     velocity: Vector3<f64>,
//     acceleration: Vector3<f64>,
//     mass: f64,
//     pressure_at_location: f64,
//     density_at_location: f64
// }
