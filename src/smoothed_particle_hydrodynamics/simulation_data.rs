use cgmath::{InnerSpace, Vector4, Zero};
use std::f64::consts::PI;
use ndarray::{Array1, Zip};
use std::time::Duration;

use crate::smoothed_particle_hydrodynamics::simulation_definition::SimulationDefinition;
use crate::smoothed_particle_hydrodynamics::simulation_space::SimulationSpace;


pub struct SimulationData {
    simulation_space: SimulationSpace,
    simulation_definition: SimulationDefinition,
    current_time: Duration,
    pub is_finished: bool,
    //voxel_pixel_map: HashMap<u16, Vec<u32>>
}

impl SimulationData {
    pub fn new(simulation_definition: SimulationDefinition) -> SimulationData {
        SimulationData {
            simulation_space: SimulationSpace::new(simulation_definition.num_particles as usize),
            simulation_definition,
            current_time: Duration::default(),
            is_finished: false,

            //voxel_pixel_map: HashMap::new()
        }
    }

    pub fn step(&mut self) {
        let densities = self.calculate_densities();
        let pressures = self.calculate_pressures(&densities);
        let pressure_terms = self.calculate_pressure_grad_terms(&densities, &pressures);
        let viscosity_terms = self.calculate_viscosity_terms(&densities);
        let accelerations = self.calculate_accelerations(&pressure_terms, &viscosity_terms);
        let new_velocities = self.calculate_velocities( &accelerations);
        let new_positions = self.calculate_positions(&new_velocities);

        //TODO - Check for collisions with walls and floor

        self.simulation_space.velocities = new_velocities;
        self.simulation_space.positions = new_positions;
        self.simulation_space.accelerations = accelerations;

        self.current_time += self.simulation_definition.time_step;
        self.is_finished = self.current_time > self.simulation_definition.sim_length;
    }

    pub fn trace_first_pixel(&self) {
        println!("{:?}", self.simulation_space.positions[0])
    }

    fn calculate_positions(&self, velocities: &Array1<Vector4<f64>>) -> Array1<Vector4<f64>> {
        Zip::from(&self.simulation_space.positions)
            .and(velocities)
            .par_map_collect(|&current_position, &current_velocity| current_position + current_velocity * self.simulation_definition.time_step.as_secs_f64())
    }

    fn calculate_velocities(&self, accelerations: &Array1<Vector4<f64>>) -> Array1<Vector4<f64>> {
        Zip::from(&self.simulation_space.positions)
            .and(accelerations)
            .par_map_collect(|&current_velocity, &current_acceleration| current_velocity + current_acceleration * self.simulation_definition.time_step.as_secs_f64())
    }

    fn calculate_accelerations(&self, pressure_terms: &Array1<Vector4<f64>>, viscosity_terms: &Array1<Vector4<f64>>) -> Array1<Vector4<f64>> {
        Zip::from(pressure_terms)
            .and(viscosity_terms)
            .par_map_collect(|&pressure_term, &viscosity_term| self.simulation_definition.gravity + pressure_term + viscosity_term)
    }

    fn calculate_viscosity_terms(&self, densities: &Array1<f64>) -> Array1<Vector4<f64>> {
        Zip::from(&self.simulation_space.positions)
            .and(densities)
            .and(&self.simulation_space.velocities)
            .par_map_collect(|&current_position, &current_density, &current_velocity | {
                Zip::from(&self.simulation_space.positions)
                    .and(densities)
                    .and(&self.simulation_space.velocities)
                    .fold(Vector4::zero(), | acc, &other_position, &other_density, &other_velocity | {
                        if self.is_in_interaction_radius_and_not_self(current_position, other_position) {
                            acc + ((self.simulation_definition.viscous_constant / current_density) * self.simulation_definition.particle_mass_kg * ((other_velocity - current_velocity) / other_density) * self.laplacian_smooth(current_position, other_position))
                        } else {
                            acc
                        }
                    })
            })
    }

    fn calculate_pressure_grad_terms(&self, densities: &Array1<f64>, pressures: &Array1<f64>) -> Array1<Vector4<f64>> {
        Zip::from(&self.simulation_space.positions)
            .and(densities)
            .and(pressures)
            .par_map_collect(|&current_position, &current_density, &current_pressure | {
                Zip::from(&self.simulation_space.positions)
                    .and(densities)
                    .and(pressures)
                    .fold(Vector4::zero(), | acc, &other_position, &other_density, &other_pressure | {
                        if self.is_in_interaction_radius_and_not_self(current_position, other_position) {
                            acc + self.simulation_definition.particle_mass_kg * ((current_pressure / current_density.powi(2)) + (other_pressure / other_density.powi(2)))* self.grad_smooth(current_position, other_position)
                        } else {
                            acc
                        }
                    })
            })
    }

    fn calculate_pressures(&self, densities: &Array1<f64>) -> Array1<f64> {
        densities.map(|&density| self.simulation_definition.fluid_constant * (density - self.simulation_definition.rho_zero))
    }

    fn calculate_densities(&self) -> Array1<f64> {
        // I don't like zipping this single thing but...
        Zip::from(&self.simulation_space.positions).par_map_collect(|&current_position| {
            self.simulation_space.positions.fold(0.0, |acc, &other_position| {
                if self.is_in_interaction_radius_and_not_self(current_position, other_position) {
                    acc + self.simulation_definition.particle_mass_kg * self.smooth(current_position, other_position)
                } else {
                    acc
                }
            })
        })
    }

    fn is_in_interaction_radius_and_not_self(&self, current: Vector4<f64>, other: Vector4<f64>) -> bool {
        current != other && (current - other).magnitude() < self.simulation_definition.smoothing_radius
    }

    // There are some terms reused between smooth, grad_smooth and laplacian_smooth, this could be optimised
    // Not to mention, some of these terms are constants...
    fn smooth(&self, current_position: Vector4<f64>, other_position: Vector4<f64>) -> f64 {
        (315.0 / (64.0 * PI * (self.simulation_definition.smoothing_radius.powi(9)))) * (self.simulation_definition.smoothing_radius.powi(2) - (current_position - other_position).magnitude2()).powi(3)
    }

    fn grad_smooth(&self, current_position: Vector4<f64>, other_position: Vector4<f64>) -> Vector4<f64> {
        (-45.0 / (PI * (self.simulation_definition.smoothing_radius.powi(6))))
            * (self.simulation_definition.smoothing_radius - (current_position - other_position).magnitude2()).powi(2)
            * ((current_position - other_position) / (current_position - other_position).magnitude2())
    }

    fn laplacian_smooth(&self, current_position: Vector4<f64>, other_position: Vector4<f64>) -> f64 {
        (45.0 / (PI * (self.simulation_definition.smoothing_radius.powi(6))))
            * (self.simulation_definition.smoothing_radius - (current_position - other_position).magnitude())
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
