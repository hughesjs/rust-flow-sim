use cgmath::{Vector4};
use crate::SimulationFloat;

pub trait SmoothingKernel: Sync + Send{
    fn new(smoothing_radius: SimulationFloat) -> Self;
    fn kernel(&self, current_location: Vector4<SimulationFloat>, other_location: Vector4<SimulationFloat>) -> SimulationFloat;
    fn kernel_grad(&self, current_location: Vector4<SimulationFloat>, other_location: Vector4<SimulationFloat>) -> Vector4<SimulationFloat>;
    fn kernel_laplacian(&self, current_location: Vector4<SimulationFloat>, other_location: Vector4<SimulationFloat>) -> SimulationFloat;
}