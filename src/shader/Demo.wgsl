
struct vo {
    @builtin(position) position: vec4<f32>;
    @location(0) uv: vec2<f32>;
    @location(1) color: vec4<f32>;
};

@vertex
fn vertex_main(@location(0) position: vec3<f32>,@location(1) color:vec4<f32) -> vo{
    vo.position = vec4<f32>(position * 2.0 - 1.0,0.0,1.0);
    vo.uv = vec2<f32>(position.x,position.y);
    vo.color = color;
    return vo;
}

@fragment
fn fragment_main(vo: vo) -> @location(0) vec4<f32> {
    return vo.color;
}