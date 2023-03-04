use std::f64::consts::PI;
use cgmath::{InnerSpace, Vector4, Zero};
use crate::{SimulationFloat};
use crate::smoothed_particle_hydrodynamics::smoothing_kernels::smoothing_kernel::SmoothingKernel;

pub struct PolySixKernel {
    kernel_coefficient:  SimulationFloat,
    grad_coefficient: SimulationFloat,
    laplacian_coefficient: SimulationFloat,
    radius: SimulationFloat,
    radius_squared: SimulationFloat
}

impl SmoothingKernel for PolySixKernel  {
    fn new(smoothing_radius: SimulationFloat) -> PolySixKernel {
        PolySixKernel {
            kernel_coefficient: (315.0 / (64.0 * PI * (smoothing_radius.powi(9)))),
            grad_coefficient: (-45.0 / (PI * (smoothing_radius.powi(6)))),
            laplacian_coefficient: (45.0 / (PI * (smoothing_radius.powi(6)))),
            radius: 0.0,
            radius_squared: 0.0,
        }
    }

    fn kernel(&self, current_location: Vector4<SimulationFloat>, other_location: Vector4<SimulationFloat>) -> SimulationFloat {
        let displacement = current_location - other_location;
        let distance = displacement.magnitude();
        let distance_squared = distance.powi(2);

        if distance > self.radius {
            0.0
        } else {
            self.kernel_coefficient * (self.radius_squared - distance_squared).powi(3)
        }

    }

    fn kernel_grad(&self, current_location: Vector4<SimulationFloat>, other_location: Vector4<SimulationFloat>) -> Vector4<SimulationFloat> {
        let displacement = current_location - other_location;
        let unit_vector = displacement.normalize();
        let distance = displacement.magnitude();

        if distance > self.radius {
            Vector4::zero()
        } else {
            self.grad_coefficient * (self.radius - distance).powi(2) * unit_vector
        }
    }

    fn kernel_laplacian(&self, current_location: Vector4<SimulationFloat>, other_location: Vector4<SimulationFloat>) -> SimulationFloat {
        let displacement = current_location - other_location;
        let distance = displacement.magnitude();

        if distance > self.radius {
            0.0
        } else {
            self.laplacian_coefficient * (self.radius - distance)
        }
    }
}