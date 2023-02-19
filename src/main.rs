//use std::collections::HashMap;
use std::f32::consts::PI;
use cgmath::{InnerSpace, Vector3, Vector4, Zero};
use rand::Rng;

const NUM_PARTICLES: u32 = 1000; // 1.2 Litres
const BUFFER_LEN: usize = NUM_PARTICLES as usize;
const BOX_DIMENSIONS_M: Vector3<f32> = Vector3::new(0.30, 0.10, 0.10);
const SMOOTHING_RADIUS: f32 = NUM_PARTICLES as f32 / (0.3 * 0.1 * 0.1 * 60.0);
const TIME_STEP_SECONDS: f32 = 0.01;
const GRAVITY: Vector4<f32> = Vector4::new(0.0, -9.8, 0.0, 0.0);
const SIM_LENGTH_SECONDS: f32 = 10.0;
const PARTICLE_MASS_KG: f32 = 0.001; //1mL of water
const FLUID_CONST: f32 = 2.2; // Possibly needs a x10^9
const VISCOUS_CONST: f32 = 0.0001;
const RHO_ZERO: f32 = 1000.0;

//TODO - Voxelise

fn main() {
	let mut sim_data: SimulationData = setup();
	let mut current_time: f32 = 0.0;
	while current_time < SIM_LENGTH_SECONDS {
		sim_step(&mut sim_data);
		current_time += TIME_STEP_SECONDS;
		println!("Time: {}s", current_time);
	}
}

fn sim_step(sim_data: &mut SimulationData) {

	let densities = sim_data.simulation_space.positions.iter().enumerate()
		.map(|(i, _)| calculate_density_at_point(&sim_data.simulation_space, i))
		.collect::<Vec<f32>>();

	let pressures = sim_data.simulation_space.positions.iter().enumerate()
		.map(|(i, _)| calculate_pressure_at_point(densities[i]))
		.collect::<Vec<f32>>();

	let pressure_grad_terms = sim_data.simulation_space.positions.iter().enumerate()
		.map(|(i, _)| calculate_pressure_grad_term_at_point(&sim_data.simulation_space, &pressures, &densities, i))
		.collect::<Vec<Vector4<f32>>>();

	let viscosity_terms = sim_data.simulation_space.positions.iter().enumerate()
		.map(|(i, _)| calculate_viscosity_term_at_point(&sim_data.simulation_space, &densities, i))
		.collect::<Vec<Vector4<f32>>>();

	let accelerations = sim_data.simulation_space.positions.iter().enumerate()
		.map(|(i, _)| calculate_acceleration_at_point(&pressure_grad_terms, &viscosity_terms, i))
		.collect::<Vec<Vector4<f32>>>();

	let new_velocities = sim_data.simulation_space.velocities.iter().enumerate()
		.map(|(i, v)| v + accelerations[i] * TIME_STEP_SECONDS)
		.collect::<Vec<Vector4<f32>>>();

	let new_positions = sim_data.simulation_space.positions.iter().enumerate()
		.map(|(i, p)| p + new_velocities[i] * TIME_STEP_SECONDS)
		.collect::<Vec<Vector4<f32>>>();

	//TODO - Check for collisions with walls and floor

	// This will make rustaceans cry
	sim_data.simulation_space.velocities = new_velocities.try_into().unwrap();
	sim_data.simulation_space.positions = new_positions.try_into().unwrap();
}

fn calculate_acceleration_at_point(pressure_terms: &Vec<Vector4<f32>>, viscosity_terms: &Vec<Vector4<f32>>, i: usize) -> Vector4<f32> {
	let acceleration: Vector4<f32> = GRAVITY + pressure_terms[i] + viscosity_terms[i];
	return acceleration;
}

fn calculate_viscosity_term_at_point(sim_space: &SimulationSpace, densities: &Vec<f32>, i: usize) -> Vector4<f32> {
	let viscosity_term: Vector4<f32> = sim_space.positions.iter()
		.filter(|&x| is_in_interaction_radius_and_not_self(x, sim_space.positions[i]))
		.enumerate()// exclude self
		.fold(Vector4::zero(), |acc: Vector4<f32>, (j, _)| {
			acc + (VISCOUS_CONST / densities[j]) * PARTICLE_MASS_KG * (sim_space.velocities[j] - sim_space.velocities[i]) / densities[j] * laplacian_smooth(sim_space.positions[i], sim_space.positions[j])
		});
	return viscosity_term;
}

fn calculate_pressure_at_point(density: f32) -> f32 {
	let pressure_at_point: f32 = FLUID_CONST * (density - RHO_ZERO);
	return pressure_at_point;
}

fn calculate_pressure_grad_term_at_point(sim_space: &SimulationSpace, pressures: &Vec<f32>, densities: &Vec<f32>, i: usize) -> Vector4<f32> {
	let pressure_grad_term: Vector4<f32> = sim_space.positions.iter()
		.filter(|&x| is_in_interaction_radius_and_not_self(x, sim_space.positions[i]))
		.enumerate()// exclude self
		.fold(Vector4::zero(), |mut acc: Vector4<f32>, (j, _)| {
			acc + PARTICLE_MASS_KG * (pressures[i] / (densities[i].powi(2)) + pressures[j] / (densities[j].powi(2))) * grad_smooth(sim_space.positions[i], sim_space.positions[j])
		});
	return pressure_grad_term;
}

fn calculate_density_at_point(sim_space: &SimulationSpace, i: usize) -> f32 {
	let density = sim_space.positions.iter()
		.filter(|&x| is_in_interaction_radius_and_not_self(x, sim_space.positions[i]))
		.enumerate()// exclude self
		.fold(0.0, |mut acc: f32, (j, _)| acc + PARTICLE_MASS_KG * smooth(sim_space.positions[i], sim_space.positions[j]));

	return density;
}

fn is_in_interaction_radius_and_not_self(current: &Vector4<f32>, other: Vector4<f32>) -> bool {
	return *current != other && (current - other).magnitude() < SMOOTHING_RADIUS;
}

// There are some terms reused between smooth, grad_smooth and laplacian_smooth, this could be optimised
// Not to mention, some of these terms are constants...
fn smooth(current_position: Vector4<f32>, other_position: Vector4<f32>) -> f32
{
	return (315.0/(64.0*PI*(SMOOTHING_RADIUS.powi(9)))) * (SMOOTHING_RADIUS.powi(2) - (current_position - other_position).magnitude2()).powi(3);
}

fn grad_smooth(current_position: Vector4<f32>, other_position: Vector4<f32>) -> Vector4<f32>
{

	return (-45.0/(PI*(SMOOTHING_RADIUS.powi(6)))) * (SMOOTHING_RADIUS - (current_position - other_position).magnitude2()).powi(2) * ((current_position - other_position) / (current_position - other_position).magnitude2());
}

fn laplacian_smooth(current_position: Vector4<f32>, other_position: Vector4<f32>) -> f32
{
	return (45.0/(PI*(SMOOTHING_RADIUS.powi(6)))) * (SMOOTHING_RADIUS - (current_position - other_position).magnitude());
}

fn setup() -> SimulationData {
	let sim_space = SimulationSpace{
		positions: get_initial_positions(),
		velocities: get_initial_velocities(),
		accelerations: get_initial_accelerations()
	};

	let sim_data = SimulationData {
		simulation_space: sim_space,
		//voxel_pixel_map: HashMap::new()
	};

	return sim_data;
}

fn get_initial_positions() -> [Vector4<f32>; BUFFER_LEN] {
	let mut rng = rand::thread_rng();
	let mut positions: [Vector4<f32>; BUFFER_LEN] = [Vector4::zero(); BUFFER_LEN];
	for i in 0..BUFFER_LEN {
		positions[i] = Vector4::new(
			rng.gen_range(0.0..0.1),
			rng.gen_range(0.0..0.01),
			rng.gen_range(4.9..5.0),
			0.0
		);
	}
	return positions;
}

fn get_initial_velocities() -> [Vector4<f32>; BUFFER_LEN] {
	return [ Vector4::zero(); BUFFER_LEN];
}

fn get_initial_accelerations() -> [Vector4<f32>; BUFFER_LEN] {
	return [ Vector4::zero(); BUFFER_LEN];
}

struct SimulationData {
	simulation_space: SimulationSpace,
	//voxel_pixel_map: HashMap<u16, Vec<u32>>
}

struct SimulationSpace {
	positions: [Vector4<f32>; BUFFER_LEN as usize],
	velocities: [Vector4<f32>; BUFFER_LEN as usize],
	accelerations: [Vector4<f32>; BUFFER_LEN as usize]
}