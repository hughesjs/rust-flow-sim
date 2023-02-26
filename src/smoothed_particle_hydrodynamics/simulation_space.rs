use cgmath::{Vector4, Zero};
use ndarray::Array1;
use rand::Rng;

pub(crate) struct SimulationSpace {
    pub(crate) positions: Array1<Vector4<f64>>,
    pub(crate) velocities:  Array1<Vector4<f64>>,
    pub(crate) accelerations:  Array1<Vector4<f64>>,
}

impl SimulationSpace {
    pub fn new(buffer_len: usize) -> Self {
        SimulationSpace {
            positions: Self::get_initial_positions(buffer_len),
            velocities: Self::get_initial_velocities(buffer_len),
            accelerations: Self::get_initial_accelerations(buffer_len),
        }
    }

    fn get_initial_positions(buffer_len: usize) -> Array1<Vector4<f64>>  {
        let mut rng = rand::thread_rng();
        let mut positions: Array1<Vector4<f64>> = Array1::from(vec![Vector4::zero(); buffer_len]);
        for i in 0..buffer_len {
            positions[i] = Vector4::new(
                rng.gen_range(0.0..0.1),
                rng.gen_range(0.0..0.1),
                rng.gen_range(4.9..5.0),
                0.0,
            );
        }
        positions
    }

    fn get_initial_velocities(buffer_len: usize) ->  Array1<Vector4<f64>> {
        Array1::from(vec![Vector4::zero(); buffer_len])
    }

    fn get_initial_accelerations(buffer_len: usize) ->  Array1<Vector4<f64>> {
        Array1::from(vec![Vector4::zero(); buffer_len])
    }
}