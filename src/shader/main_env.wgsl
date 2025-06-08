fn get_pixel(uv:vec2<f32>,sun_color: vec3<f32>,moon_color: vec3<f32>,ambient:vec3<f32>)->vec3<f32>{
    let ball = get_ball_material();
    let dir = calculate_dir(uv);
    let res = is_hit_ball(POS(),dir);
    let camera_pos = POS();
    var r:vec3<f32>;
    var p:vec3<f32>;
    var nor:vec3<f32>;
    var size :i32;
    var trash:bool;
    if(res.hit_type == 2u){
        nor = normalize(res.position - get_ball_pos());
        r = reflect(dir,nor);
        p = res.position;
        trash = false; //只有球体需要计算环境光
    }
    else{
        r = dir;
        p = camera_pos;
        trash = true;
    }
    let scene_color = get_color(p,r,sun_color,moon_color,false,trash);
    let env_color = get_env_color(p,nor,dir,sun_color,moon_color,trash);

    var color = vec3(0.0);
    if(res.hit_type == 2u){
        let base_color = ball.albedo;
        let F0 = ball.f0.x;
        let roughness = ball.roughness;
        let metallic = ball.metallic;

        let light_data = get_light_data(get_ball_pos());//点光源光照
        let light_pos = light_data[0];
        let light_color = light_data[1];
        let dist = length(light_pos - p);
        let light_pbr = pbr_shading(base_color,roughness,F0,metallic,nor,normalize(light_pos - p),-dir,normalize(-dir + normalize(light_pos - p)));
        color += light_pbr.diffuse * light_color / (dist * dist * 0.07 + dist * 0.14 + 1.0) + light_pbr.specular * light_color / (dist * dist * 0.07 + dist * 0.14 + 1.0);

        let pbr_res = pbr_shading(base_color,roughness,F0,metallic,nor,uniforms.sun_dir,-dir,normalize(-dir + uniforms.sun_dir));
        let scene_pbr_res = pbr_shading(base_color,roughness,F0,metallic,nor,r,-dir,normalize(-dir + r));
        let n = dot(nor,uniforms.sun_dir) * 0.5 + 0.5;
        let H = normalize(r - dir);
        let VoH = max(dot(-dir, H), 0.0);
        let t = dot(-dir,nor);
        color += pbr_res.diffuse * sun_color * 0.15 + pbr_res.specular * sun_color;
        let ambient = sun_color * (n + 0.5) * 0.3 * base_color;
        color += ambient * 0.05 * max((1 - metallic),0.2);
        let env_res = pbr_env(base_color,roughness,ball.f0,metallic,nor,normalize(r),-dir,normalize(-dir + r));
        color += (mix(scene_color * env_res,env_color,roughness)) * (1 - roughness);

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
    @location(4) ambient:vec3<f32>,
    @location(5) ball_pos:vec3<f32>,
    @location(6) camera_pos:vec3<f32>,
    @location(7) light_pos:vec3<f32>,
    @location(8) light_color:vec3<f32>,
};

@vertex
fn vertex_main(@location(0) position: vec2<f32>,@location(1) color:vec4<f32>) -> VertexOutput{
    var output: VertexOutput;
    output.position = vec4<f32>(position * 2.0 - 1.0,0.0,1.0);
    output.uv = vec2<f32>(position.x,position.y);
    output.color = color;
    output.sun_color = get_sun_exposure() * get_sun_tint();
    output.moon_color = get_moon_exposure() * get_moon_tint();
    output.ambient = vec3<f32>(0.0);
    return output;
}

//VERTEX


@fragment
fn fragment_main(@location(0) uv: vec2<f32>,@location(1) color: vec4<f32>,
                 @location(2) sun_color:vec3<f32>,@location(3) moon_color:vec3<f32>,
                 @location(4) ambient:vec3<f32>)
                 -> @location(0) vec4f {
    let pixel = pow(
        get_pixel(uv, sun_color, moon_color, ambient),
        vec3f(0.65)
    );
    return vec4f(pixel, 1.0);
}


//FRAGMENT