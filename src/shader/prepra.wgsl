fn get_pixel(uv:vec2<f32>,sun_color: vec3<f32>,moon_color: vec3<f32>,ambient:vec3<f32>)->vec3<f32>{
    let dir = unproject_sky(uv);
    let pos = get_ball_pos();
    let color = get_color_without_sun(pos,dir,sun_color,moon_color,false);
    return color;
}


//VERTEX

struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) uv: vec2<f32>,
    @location(1) color: vec4<f32>,
    @location(2) sun_color:vec3<f32>,
    @location(3) moon_color:vec3<f32>,
    @location(4) ambient:vec3<f32>
};

@vertex
fn vertex_main(@location(0) position: vec2<f32>,@location(1) color:vec4<f32>) -> VertexOutput{
    var output: VertexOutput;
    output.position = vec4<f32>(position * 2.0 - 1.0,0.0,1.0);
    output.uv = vec2<f32>(position.x,position.y);
    output.color = color;
    output.sun_color = get_sun_exposure() * get_sun_tint();
    output.moon_color = get_moon_exposure() * get_moon_tint();
    output.ambient = vec3(0.15,0.2,0.1);
    return output;
}

//VERTEX


@fragment
fn fragment_main(@location(0) uv: vec2<f32>,@location(1) color: vec4<f32>,@location(2) sun_color:vec3<f32>,@location(3) moon_color:vec3<f32>,@location(4) ambient:vec3<f32>) -> @location(0) vec4<f32> {
    return vec4<f32>(pow(get_pixel(uv,sun_color,moon_color,ambient),vec3<f32>(1.0)),1.0);

}


//FRAGMENT