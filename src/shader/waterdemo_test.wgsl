#ifndef WATERDEMO_TEST
#define WATERDEMO_TEST

#include "global.wgsl"
#include "atmosphere.wgsl"
#include "pbr.wgsl"

// ================================
// 光线追踪/海洋函数
// ================================
fn sea_octave(uv: vec2<f32>, choppy: f32) -> f32 {
    var uv_mut = uv;
    uv_mut += noise(uv);
    let wv = 1.0 - abs(sin(uv_mut));
    let swv = abs(cos(uv_mut));
    return pow(1.0 - pow(mix(wv, swv, wv).x * mix(wv, swv, wv).y, 0.65), choppy);
}

fn map(p: vec3<f32>) -> f32 {
    var freq: f32 = SEA_FREQ;
    var amp: f32 = SEA_HEIGHT;
    var choppy: f32 = SEA_CHOPPY;
    var uv: vec2<f32> = p.xz + vec2(0.5);
    uv.x *= 0.75;
    
    var h: f32 = 0.0;
    for(var i: i32 = 0; i < ITER_GEOMETRY; i++) { 
        let d = sea_octave((uv + SEA_TIME()) * freq, choppy) +
                sea_octave((uv - SEA_TIME()) * freq, choppy);
        h += d * amp;
        uv = uv * octave_m;
        freq *= 1.9;
        amp *= 0.22;
        choppy = mix(choppy, 1.0, 0.2);
    }
    return p.y - h;
}

fn map_detail(p: vec3<f32>) -> f32 {
    var freq: f32 = SEA_FREQ;
    var amp: f32 = SEA_HEIGHT;
    var choppy: f32 = SEA_CHOPPY;
    var uv: vec2<f32> = p.xz + vec2(0.5);
    uv.x *= 0.75;
    
    var h: f32 = 0.0;
    for(var i: i32 = 0; i < ITER_FRAGMENT; i++) { 
        let d = sea_octave((uv + SEA_TIME()) * freq, choppy) +
                sea_octave((uv - SEA_TIME()) * freq, choppy);
        h += d * amp;
        uv = uv * octave_m;
        freq *= 1.9;
        amp *= 0.22;
        choppy = mix(choppy, 1.0, 0.2);
    }
    return p.y - h;
}

fn getNormal(p:vec3<f32>,eps:f32) -> vec3<f32> {
    var n:vec3<f32>;
    n.y = map_detail(p);
    n.x = map_detail(vec3<f32>(p.x + eps, p.y, p.z)) - n.y;
    n.z = map_detail(vec3<f32>(p.x, p.y, p.z + eps)) - n.y;
    n.y = eps;
    return normalize(n);
}

fn heightMapTracing(ori: vec3<f32>, dir: vec3<f32>) -> TracingResult {
    var tm: f32 = 0.0;
    var tx: f32 = TX;
    var p: vec3<f32>;
    var hx = map(ori + dir * tx);
    if(hx > 0.0) {
        return TracingResult(tx, ori + dir * tx, 1u);
    }
    
    tx = min((POS().y + 2.0) / dir.y,tx);
    tm = max((POS().y - 3.0) / dir.y,0.0);
    var hm = map(ori + dir * tm);
    hx = map(ori + dir * tx);
    for(var i: i32 = 0; i < NUM_STEPS; i++) {
        let tmid: f32 = mix(tm, tx, hm / (hm - hx));
        p = ori + dir * tmid;
        let hmid = map(p);
        
        if(hmid < 0.0) {
            tx = tmid;
            hx = hmid;
        } else {
            tm = tmid;
            hm = hmid;
        }
        if(abs(hmid) < EPSILON / 2.0) { break; }
    }
    return TracingResult(
        mix(tm, tx, hm / (hm - hx)),  
        p,
        0u                         
    );
}

fn getSeaColor(p: vec3<f32>, n: vec3<f32>, l: vec3<f32>, eye: vec3<f32>, dist: vec3<f32>,
                sun_color:vec3<f32>,moon_color:vec3<f32>,
                ) -> vec3<f32> {  

    let fresnel = fresnel_schlick(dot(n, -eye), vec3(0.02));
    //反射和折射成分
    let rl = reflect(eye, n);
    let h = normalize(-eye + l);
    let pbr_res = pbr_shading_water(SEA_WATER_COLOR * sun_color * 0.12,0.002,0.02,0.0,n,l,-eye,h);
    let reflected: vec3<f32> = atmosphere_scattering(rl,sun_color,uniforms.sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0)) + draw_sun(rl,sun_color);
    let refracted: vec3<f32> = SEA_BASE  + (pbr_res.diffuse * sun_color) * 0.05; 
    
    // 基础颜色混合
    var color: vec3<f32> = mix(refracted, reflected, fresnel);
    
    // 距离衰减效果
    let atten: f32 = max(1.0 - dot(dist, dist) * 0.001, 0.0);
    color += (SEA_WATER_COLOR * (p.y - SEA_HEIGHT)) * 0.18 * atten;
    
    //高光添加
    color += ACESToneMapping(pbr_res.specular * sun_color,0.25);
    
    return ACESToneMapping(color,0.72);
}

// ================================
// 球体相关函数
// ================================
fn get_ball_pos() -> vec3<f32> {
    return ball_pos * fromEuler(ANGLE());
}

fn get_ball_material() -> Material{
    var ball:Material = ball_material;
    ball.metallic = uniforms.ball_data[1];
    ball.roughness = uniforms.ball_data[0];
    ball.f0 = mix(vec3(0.04),ball.albedo,ball.metallic);
    return ball;
}

fn is_hit_ball(ori: vec3<f32>, dir: vec3<f32>) -> TracingResult {
    let dis = length(ori - get_ball_pos());
    let d = dis * dot(normalize(get_ball_pos() - ori), dir);
    let l = sqrt(dis * dis - d * d);
    if(l < ball_radius && d > 0.0){
        let t = sqrt(sqr(ball_radius) - sqr(l));
        return TracingResult(d - t, ori + dir * (d - t), 2u);
    }
    return TracingResult(0.0,vec3<f32>(0.0,0.0,0.0),0u);
}

fn get_light_data(ball_pos: vec3<f32>) -> mat2x3<f32> {
    let dir = vec3(cos(uniforms.iTime * 0.1 * light_speed),sin(uniforms.iTime * 0.16 * light_speed),sin(uniforms.iTime * 0.1 * light_speed));
    let color = abs(normalize(dir)) * 3.0 * light_strength;
    return mat2x3<f32>(dir * light_dist + ball_pos, color);
}

// ================================
// 环境纹理采样
// ================================
fn fast_gaussian_5x5(uv: vec2<f32>) -> vec3<f32> {
    let texture_size = textureDimensions(env_texture);
    let inv_size = 1.0 / vec2<f32>(texture_size);
    var result = vec3<f32>(0.0);

    for (var i: i32 = -2; i <= 2; i++) {
        let offset = vec2<f32>(f32(i), 0.0) * inv_size;
        result += textureSample(env_texture, env_sampler, uv + offset).rgb * KERNEL_5[i + 2];
    }

    for(var i: i32 = -2; i <= 2; i++){
        let offset = vec2<f32>(0.0,f32(i)) * inv_size;
        result += textureSample(env_texture, env_sampler, uv + offset).rgb * KERNEL_5[i + 2];
    }
    result /= 2.0;
    return result;
}

fn sample_env_texture(dir: vec3<f32>, trash: bool) -> vec3<f32> {
    var uv: vec2<f32>;
    if (trash) {
        uv = vec2<f32>(0.0, 0.0);
    } else {
        uv = project_sky(dir);
    }
    uv.y = 1.0 - uv.y;
    
    return fast_gaussian_5x5(uv);
}

fn get_env_color(p:vec3<f32>,n:vec3<f32>,v:vec3<f32>,sun_color: vec3<f32>,moon_color:vec3<f32>,trash:bool) ->vec3<f32> {
    let ball = get_ball_material();
    var env_color = vec3<f32>(0.0);
    let tbn = get_tbn(n);
    let base_color = ball.albedo;
    let f0 =ball.f0;
    let roughness = ball.roughness;
    let metallic = ball.metallic;
    for(var i:f32 = 0; i < env_sample_size; i+=1){
        let hash_uv = vec2<f32>(hash(vec2(i / env_sample_size) + vec2(0.62,0.32)),hash(vec2(i / env_sample_size) + vec2(0.32,0.62)));
        let microfacet_normal = tbn * sample_ggx_vndf(-(v * tbn),vec2(1.0),hash_uv);
        let rl = reflect(v,microfacet_normal);
        let h = normalize(-v + rl);
        let pbr_env_res = pbr_env(base_color,roughness,f0,metallic,n,rl,-v,h);
        let t = (trash || dot(rl,n) < 0.001);
        let color = sample_env_texture(rl,t);
        env_color += pbr_env_res * color;
    }
    env_color /= env_sample_size;
    return env_color;
}

// ================================
// 主颜色计算函数
// ================================
fn calculate_dir(uv: vec2<f32>) -> vec3<f32>{
    var ndc_uv = vec2<f32>(uv.x * 2.0 - 1.0,uv.y * 2.0 - 1.0);
    ndc_uv.x *= uniforms.iResolution.x / uniforms.iResolution.y;
    var dir = normalize(vec3<f32>(ndc_uv.x,ndc_uv.y,-FOV));
    return normalize(normalize(dir) * fromEuler(ANGLE()));
}

fn get_color(p:vec3<f32>,dir:vec3<f32>,sun_color: vec3<f32>,moon_color:vec3<f32>,trash:bool,should_draw_sun:bool) -> vec3<f32>{
    var sky_color = atmosphere_scattering(dir,sun_color,uniforms.sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0));
    sky_color += draw_sun(dir,sun_color) * select(0.0,1.0,should_draw_sun);
    sky_color = ACESToneMapping(sky_color,0.5);

    var res :TracingResult;
    var dist :vec3<f32>;
    var n : vec3<f32>;

    if(trash){
        res = TracingResult(0.0,vec3<f32>(0.0,0.0,0.0),0u);
        dist = vec3<f32>(1.0);
        n = vec3<f32>(0.0);
    }
    else{
        res = heightMapTracing(p, dir);
        dist = res.position - p;
        n = getNormal(res.position,dot(dir, dist) * EPSILON_NRM());
    }
    var sea_color = getSeaColor(res.position, n, uniforms.sun_dir, dir, dist,sun_color,moon_color);
    let alp = 1 - smoothstep(-0.02,0.0,dir.y);

    let horizon_dir = normalize(vec3(dir.xz,min(dir.y,-0.1)).xzy);
    let horizon_color = atmosphere_scattering(horizon_dir,sun_color,uniforms.sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0));

    if(trash){
        return vec3(0.0);
    }

    let horizon_fac = linear_step(0.1,1.0,exp(-75 * sqr(uniforms.sun_dir.y + 0.0496)));
    let fog_color = mix(sky_color,horizon_color,sqr(horizon_fac));
    let fog = border_fog(dist,dir);
    
    sea_color = mix(fog_color,sea_color,fog);
    sky_color = mix(fog_color,sky_color,fog);
    return mix(sky_color,sea_color,alp);
}

fn get_color_without_sun(p:vec3<f32>,dir:vec3<f32>,sun_color: vec3<f32>,moon_color:vec3<f32>,trash:bool) -> vec3<f32>{
    var sky_color = atmosphere_scattering(dir,sun_color,uniforms.sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0));
    sky_color = ACESToneMapping(sky_color,0.5);

    var res :TracingResult;
    var dist :vec3<f32>;
    var n : vec3<f32>;

    if(trash){
        res = TracingResult(0.0,vec3<f32>(0.0,0.0,0.0),0u);
        dist = vec3<f32>(1.0);
        n = vec3<f32>(0.0);
    }
    else{
        res = heightMapTracing(p, dir);
        dist = res.position - p;
        n = getNormal(res.position,dot(dir, dist) * EPSILON_NRM());
    }
    var sea_color = getSeaColor(res.position, n, uniforms.sun_dir, dir, dist,sun_color,moon_color);
    let alp = 1 - smoothstep(-0.02,0.0,dir.y);

    let horizon_dir = normalize(vec3(dir.xz,min(dir.y,-0.1)).xzy);
    let horizon_color = atmosphere_scattering(horizon_dir,sun_color,uniforms.sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0));

    if(trash){
        return vec3(0.0);
    }

    let horizon_fac = linear_step(0.1,1.0,exp(-75 * sqr(uniforms.sun_dir.y + 0.0496)));
    let fog_color = mix(sky_color,horizon_color,sqr(horizon_fac));
    let fog = border_fog(dist,dir);
    
    sea_color = mix(fog_color,sea_color,fog);
    sky_color = mix(fog_color,sky_color,fog);
    return mix(sky_color,sea_color,alp);
}
#endif

