use cgmath::{Vector3};
use crate::SimulationFloat;

pub trait SmoothingKernel: Sync + Send{
    fn new(smoothing_radius: SimulationFloat) -> Self;
    fn kernel(&self, current_location: Vector3<SimulationFloat>, other_location: Vector3<SimulationFloat>) -> SimulationFloat;
    fn kernel_grad(&self, current_location: Vector3<SimulationFloat>, other_location: Vector3<SimulationFloat>) -> Vector3<SimulationFloat>;
    fn kernel_laplacian(&self, current_location: Vector3<SimulationFloat>, other_location: Vector3<SimulationFloat>) -> SimulationFloat;
}