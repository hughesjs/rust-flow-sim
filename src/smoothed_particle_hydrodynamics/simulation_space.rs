use cgmath::{Vector3, Zero};
use ndarray::Array1;
use rand::Rng;
use crate::SimulationFloat;

pub(crate) struct SimulationSpace {
    pub(crate) positions: Array1<Vector3<SimulationFloat>>,
    pub(crate) velocities:  Array1<Vector3<SimulationFloat>>,
    pub(crate) accelerations:  Array1<Vector3<SimulationFloat>>,
}

impl SimulationSpace {
    pub fn new(buffer_len: usize) -> Self {
        SimulationSpace {
            positions: Self::get_initial_positions(buffer_len),
            velocities: Self::get_initial_velocities(buffer_len),
            accelerations: Self::get_initial_accelerations(buffer_len),
        }
    }

    fn get_initial_positions(buffer_len: usize) -> Array1<Vector3<SimulationFloat>>  {
        let mut rng = rand::thread_rng();
        let mut positions: Array1<Vector3<SimulationFloat>> = Array1::from(vec![Vector3::zero(); buffer_len]);
        for i in 0..buffer_len {
            positions[i] = Vector3::new(
                rng.gen_range(0.0..0.1),
                rng.gen_range(0.0..0.1),
                rng.gen_range(4.9..5.0)
            );
        }
        positions
    }

    fn get_initial_velocities(buffer_len: usize) ->  Array1<Vector3<SimulationFloat>> {
        Array1::from(vec![Vector3::zero(); buffer_len])
    }

    fn get_initial_accelerations(buffer_len: usize) ->  Array1<Vector3<SimulationFloat>> {
        Array1::from(vec![Vector3::zero(); buffer_len])
    }
}