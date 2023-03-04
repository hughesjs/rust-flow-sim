use cgmath::{Vector4};

#[macro_export] macro_rules! float_precision {
    () => {
        f64
    };
}

pub trait SmoothingKernel {
    fn new(smoothing_radius: float_precision!()) -> Self;
    fn kernel(&self, current_location: Vector4<float_precision!()>, other_location: Vector4<float_precision!()>) -> float_precision!();
    fn kernel_grad(&self, current_location: Vector4<float_precision!()>, other_location: Vector4<float_precision!()>) -> Vector4<float_precision!()>;
    fn kernel_laplacian(&self, current_location: Vector4<float_precision!()>, other_location: Vector4<float_precision!()>) -> float_precision!();
}