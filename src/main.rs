use std::collections::HashMap;
use std::iter::Map;
use cgmath::Vector4;

const NUM_PARTICLES: u32 = 12000;
const BUFFER_LEN: usize = NUM_PARTICLES as usize;
const SMOOTHING_RADIUS: u8 = 5;
const VOXEL_SIZE: u16 = 10;
const NUM_VOXELS: u16 = 100;

fn main() {
	let sim_data: SimulationData = setup();

}

fn setup() -> SimulationData {
	let zero_vector = Vector4::new(0.0, 0.0, 0.0, 0.0);
	let empty_vector_buffer: [Vector4<f32>; BUFFER_LEN] = [zero_vector; BUFFER_LEN];
	let sim_space = SimulationSpace{
		position: empty_vector_buffer,
		velocity: empty_vector_buffer,
		acceleration: empty_vector_buffer
	};

	let sim_data = SimulationData {
		simulation_space: sim_space,
		voxel_pixel_map: HashMap::new()
	};

	return sim_data;
}

struct SimulationData {
	simulation_space: SimulationSpace,
	voxel_pixel_map: HashMap<u16, Vec<u32>>
}

struct SimulationSpace {
	position: [Vector4<f32>; BUFFER_LEN as usize],
	velocity: [Vector4<f32>; BUFFER_LEN as usize],
	acceleration: [Vector4<f32>; BUFFER_LEN as usize]
}