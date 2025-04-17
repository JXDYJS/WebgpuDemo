fn clamp01(x: f32) -> f32 {
    return clamp(x, 0.0, 1.0);
}

fn remap(value: f32, fromMin: f32, fromMax: f32, toMin: f32, toMax: f32) -> f32 {
    return (value - fromMin) * (toMax - toMin) / (fromMax - fromMin) + toMin;
}