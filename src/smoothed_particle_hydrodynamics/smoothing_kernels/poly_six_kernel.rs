use std::f64::consts::PI;
use cgmath::{InnerSpace, Vector4};
use crate::float_precision;
use crate::smoothed_particle_hydrodynamics::smoothing_kernels::smoothing_kernel::SmoothingKernel;

pub struct PolySixKernel {
    kernel_coefficient:  float_precision!(),
    grad_coefficient: float_precision!(),
    laplacian_coefficient: float_precision!(),
    radius: float_precision!(),
    radius_squared: float_precision!()
}

impl SmoothingKernel for PolySixKernel  {
    fn new(smoothing_radius: float_precision!()) -> PolySixKernel {
        PolySixKernel {
            kernel_coefficient: (315.0 / (64.0 * PI * (smoothing_radius.powi(9)))),
            grad_coefficient: (-45.0 / (PI * (smoothing_radius.powi(6)))),
            laplacian_coefficient: (45.0 / (PI * (smoothing_radius.powi(6)))),
            radius: 0.0,
            radius_squared: 0.0,
        }
    }

    fn kernel(&self, current_location: Vector4<float_precision!()>, other_location: Vector4<float_precision!()>) -> float_precision!() {
        let displacement = current_location - other_location;
        let distance_squared = displacement.magnitude2();

         self.kernel_coefficient * (self.radius_squared - distance_squared).powi(3)
    }

    fn kernel_grad(&self, current_location: Vector4<float_precision!()>, other_location: Vector4<float_precision!()>) -> Vector4<float_precision!()> {
        let displacement = current_location - other_location;
        let unit_vector = displacement.normalize();
        let distance = displacement.magnitude();

        self.grad_coefficient * (self.radius - distance).powi(2) * unit_vector
    }

    fn kernel_laplacian(&self, current_location: Vector4<float_precision!()>, other_location: Vector4<float_precision!()>) -> float_precision!() {
        let displacement = current_location - other_location;
        let distance = displacement.magnitude();

        self.laplacian_coefficient * (self.radius - distance)
    }
}