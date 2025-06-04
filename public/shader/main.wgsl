fn get_pixel(uv:vec2<f32>,sun_color: vec3<f32>,moon_color: vec3<f32>,ambient:vec3<f32>)->vec3<f32>{
    let dir = calculate_dir(uv);
    let res = is_hit_ball(POS(),dir);
    var r:vec3<f32>;
    var p:vec3<f32>;
    var nor:vec3<f32>;
    var size :i32;
    var trash:bool;
    if(res.hit_type == 2u){
        nor = normalize(res.position - get_ball_pos());
        r = reflect(dir,nor);
        p = res.position;
        trash = false;
    }
    else{
        r = dir;
        p = POS();
        trash = true;
    }
    let scene_color = get_color(p,r,sun_color,moon_color,false);

    var color:vec3<f32>;
    if(res.hit_type == 2u){
        let base_color = ball_material.albedo;
        let F0 = ball_material.f0.x;
        let roughness = ball_material.roughness;
        let metallic = ball_material.metallic;

        let pbr_res = pbr_shading(base_color,roughness,F0,metallic,nor,uniforms.sun_dir,-dir,normalize(-dir + uniforms.sun_dir));
        let scene_pbr_res = pbr_shading(base_color,roughness,F0,metallic,nor,r,-dir,normalize(-dir + r));
        let n = dot(nor,uniforms.sun_dir) * 0.5 + 0.5;
        let H = normalize(r - dir);
        let VoH = max(dot(-dir, H), 0.0);
        let t = dot(-dir,nor);
        color = pbr_res.diffuse * sun_color * 0.15 + pbr_res.specular * sun_color +  scene_pbr_res.specular * scene_color * smoothstep(0.02,0.3,t);
        let ambient = sun_color * (n + 0.5) * 0.3 * base_color;
        color += ambient * 0.05;
        let env_res = pbr_env(base_color,roughness,F0,metallic,nor,r,-dir,normalize(-dir + r));
        color = ACESToneMapping(color,0.5);
    }
    else{
        color = scene_color;
    }
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
    return vec4<f32>(pow(get_pixel(uv,sun_color,moon_color,ambient),vec3<f32>(0.65)),1.0);

}


//FRAGMENT