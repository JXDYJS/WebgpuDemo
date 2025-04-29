
//大气部分函数
const sun_dir = normalize(vec3<f32>(0.0, 0.1,1.0));
const moon_dir = normalize(vec3<f32>(0.0, -0.3,1.0));
const SUN_I: f32 = 1.0;
const MOON_I: f32 = 1.0;
const MOON_R: f32 = 0.75;
const MOON_G: f32 = 0.83;
const MOON_B: f32 = 1.0;
const PI: f32 = 3.141592;
const EPSILON: f32 = 1e-3;
const TAU: f32 = 6.283185307179586;
const DEGREE: f32 = PI / 180.0;
const TX:f32 = 300.0;

fn FADE()-> f32{
    if(sun_dir.y < 0.18){
        return 0.37 + 1.2 * max(0.0, -sun_dir.y);
    }
    return 0.17;
}

fn SUN_WEIGHT()-> f32{
    return pow(clamp(1.0 - FADE() * abs(sun_dir.y - 0.18) ,0.0,1.0),2.0);
}

fn SUN_RAISE()-> f32{
    return step(0.0001,sun_dir.x) * SUN_WEIGHT();
}

fn SUN_SET()-> f32{
    return step(0.0001,-sun_dir.x) * SUN_WEIGHT();
}

// 天体参数
const SUN_ANGULAR_RADIUS: f32 = 0.275 * DEGREE;
const MOON_ANGULAR_RADIUS: f32 = 0.278 * DEGREE;
const SUNLIGHT_COLOR: vec3<f32> = vec3(1.051, 0.985, 0.940);
const MOONLIGHT_COLOR: vec3<f32> = vec3(0.05, 0.1, 0.15);

// 大气参数
const PLANET_RADIUS: f32 = 6371e3;
const ATMOSPHERE_INNER_RADIUS: f32 = PLANET_RADIUS - 1e3;
const ATMOSPHERE_OUTER_RADIUS: f32 = PLANET_RADIUS + 110e3;
const ATMOSPHERE_THICKNESS: f32 = ATMOSPHERE_OUTER_RADIUS - ATMOSPHERE_INNER_RADIUS;

// LUT参数
const SCATTERING_RES: vec3<f32> = vec3(16.0, 64.0, 32.0);
const MIN_MU_S: f32 = -0.35;

fn linear_step(edge0: f32, edge1: f32, x: f32) -> f32 {
    return clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
}

fn pulse(a: f32, b: f32, x: f32) -> f32 {
    return step(a, x) - step(b, x);
}

fn from_srgb(srgb: vec3<f32>) -> vec3<f32> {
    let linear = select(
        srgb / 12.92,
        pow((srgb + 0.055) / 1.055, vec3<f32>(2.4)),
        srgb > vec3(0.04045)
    );
    return linear;
}

fn smoothstep(low: f32, high: f32, x: f32) -> f32 {
    let t = clamp((x - low) / (high - low), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
}

fn cube(a: f32) -> f32 {
    return a * a * a;
}

fn sqr(a:f32)->f32{
    return a * a;
}

fn fast_acos(x:f32) -> f32 {
	let C0 = 1.57018;
	let C1 = -0.201877;
	let C2 = 0.0464619;

	var res = (C2 * abs(x) + C1) * abs(x) + C0; // p(x)
	res *= sqrt(1.0 - abs(x));

	if(x>=0.0){
		return res;
    }
    else{
        return PI - res;
    }
}

fn cubic_length(v:vec2<f32>)->f32 {
	return pow(cube(abs(v.x)) + cube(abs(v.y)), 0.333333);
}

fn dampen(x:f32)->f32{
	let t = clamp(x,0.0,1.0);
	return t * (2.0 - t);
}

fn ACESToneMapping(color: vec3<f32>, adapted_lum: f32) -> vec3<f32> {
    const A: f32 = 2.51;
    const B: f32 = 0.03;
    const C: f32 = 2.43;
    const D: f32 = 0.59;
    const E: f32 = 0.14;

    let scaled_color = color * adapted_lum;
    return (scaled_color * (A * scaled_color + B)) / 
           (scaled_color * (C * scaled_color + D) + E);
}

fn get_uv_from_unit_range(val: f32, resolution: f32) -> f32 {
    return clamp(0.5 / resolution + val * (1.0 - 1.0 / resolution), 0.0, 1.0);
}

fn intersect_sphere(mu: f32, radius: f32, outer_radius: f32) -> vec2<f32> {
    let discriminant = radius * radius * (mu * mu - 1.0) + outer_radius * outer_radius;
    let sqrt_d = sqrt(max(discriminant, 0.0));
    return vec2(-radius * mu - sqrt_d, -radius * mu + sqrt_d);
}

fn atmosphere_mie_phase(nu: f32) -> f32 {
    const G: f32 = 0.76;
    let gg = G * G;
    return (3.0 * (1.0 - gg)) / (8.0 * PI * (2.0 + gg)) * (1.0 + nu * nu) / pow(1.0 + gg - 2.0 * G * nu, 1.5);
}

fn get_sun_exposure() -> f32 {
    let time_sunset = SUN_SET();
    let time_sunrise = SUN_RAISE();
    let base_scale = 7.0 * SUN_I;
    
    // 蓝调时刻计算
    let blue_hour = linear_step(
        0.05, 
        1.0, 
        exp(-190.0 * pow(sun_dir.y + 0.09604, 2.0))
    );
    
    let daytime_mul = 1.0 + 0.5 * (time_sunset + time_sunrise) + 40.0 * blue_hour;
    
    return base_scale * daytime_mul;
}

fn get_sun_tint() -> vec3<f32> {
    let time_sunset = SUN_SET();
    let time_sunrise = SUN_RAISE();
    // 蓝调时刻系数
    let blue_hour = linear_step(
        0.05, 
        1.0, 
        exp(-190.0 * pow(sun_dir.y + 0.09604, 2.0))
    );
    
    // 早晚色调
    var morning_evening_tint = vec3<f32>(1.05, 0.84, 0.93) * 1.2;
    morning_evening_tint = mix(
        vec3<f32>(1.0), 
        morning_evening_tint, 
        pow(pulse(0.17, 0.40, sun_dir.y), 2.0)
    );
    
    // 蓝调时刻色调
    var blue_hour_tint = vec3<f32>(0.95, 0.80, 1.0);
    blue_hour_tint = mix(vec3<f32>(1.0), blue_hour_tint, blue_hour);
    
    // 用户自定义色调（需传入以下常量）
    let tint_morning = from_srgb(vec3<f32>(1.0, 1.0, 1.0)); // 替换实际参数
    let tint_noon    = from_srgb(vec3<f32>(1.0, 1.0, 1.0));
    let tint_evening = from_srgb(vec3<f32>(1.0, 1.0, 1.0));
    
    var user_tint = mix(tint_noon, tint_morning, time_sunrise);
    user_tint = mix(user_tint, tint_evening, time_sunset);
    
    return morning_evening_tint * blue_hour_tint * user_tint;
}

// 月光计算 -------------------------------------------------
fn get_moon_exposure() -> f32 {
    return 0.66 * MOON_I * 1.0;
}

fn get_moon_tint() -> vec3<f32> {
    return from_srgb(vec3<f32>(MOON_R, MOON_G, MOON_B));
}

// 主散射函数 ---------------------------------------------------------------
fn atmosphere_scattering(
    ray_dir: vec3<f32>,
    sun_color: vec3<f32>,
    sun_dir: vec3<f32>,
    moon_color: vec3<f32>,
    moon_dir: vec3<f32>
) -> vec3<f32> {
    // 基础参数
    let mu = ray_dir.y;
    let nu_sun = dot(ray_dir, sun_dir);
    let nu_moon = dot(ray_dir, moon_dir);
    let mu_sun = sun_dir.y;
    let mu_moon = moon_dir.y;

    // 地平线优化
    var adjusted_mu = mu;
    let sun_step = smoothstep(-0.05, 0.1, mu_sun);
    let moon_step = smoothstep(0.05, 0.1, mu_moon);
    let horizon_mu = mix(-0.01, 0.03, clamp(sun_step + moon_step, 0.0, 1.0));
    adjusted_mu = max(mu, horizon_mu);

    // 太阳散射计算 ---------------------------------------------------------
    // nu映射
    let half_range_nu_sun = sqrt((1.0 - adjusted_mu * adjusted_mu) * (1.0 - mu_sun * mu_sun));
    let nu_min_sun = adjusted_mu * mu_sun - half_range_nu_sun;
    let nu_max_sun = adjusted_mu * mu_sun + half_range_nu_sun;
    let u_nu_sun = select(
        (nu_sun - nu_min_sun) / (nu_max_sun - nu_min_sun),
        nu_min_sun,
        nu_min_sun == nu_max_sun
    );
    let u_nu_sun_uv = get_uv_from_unit_range(u_nu_sun, SCATTERING_RES.x);

    // mu映射
    let r = PLANET_RADIUS;
    let H = sqrt(ATMOSPHERE_OUTER_RADIUS * ATMOSPHERE_OUTER_RADIUS - ATMOSPHERE_INNER_RADIUS * ATMOSPHERE_INNER_RADIUS);
    let rho = sqrt(max(r * r - ATMOSPHERE_INNER_RADIUS * ATMOSPHERE_INNER_RADIUS, 0.0));
    let rmu = r * adjusted_mu;
    let discriminant = rmu * rmu - r * r + ATMOSPHERE_INNER_RADIUS * ATMOSPHERE_INNER_RADIUS;

    var u_mu: f32;
    if (adjusted_mu < 0.0 && discriminant >= 0.0) {
        let d = -rmu - sqrt(discriminant);
        let d_min = r - ATMOSPHERE_INNER_RADIUS;
        let d_max = rho;
        u_mu = (d - d_min) / (d_max - d_min);
        u_mu = get_uv_from_unit_range(u_mu, SCATTERING_RES.y / 2.0);
        u_mu = 0.5 - 0.5 * u_mu;
    } else {
        let d = -rmu + sqrt(discriminant + H * H);
        let d_min = ATMOSPHERE_OUTER_RADIUS - r;
        let d_max = rho + H;
        u_mu = (d - d_min) / (d_max - d_min);
        u_mu = get_uv_from_unit_range(u_mu, SCATTERING_RES.y / 2.0);
        u_mu = 0.5 + 0.5 * u_mu;
    }

    // mu_s映射
    let intersect_sun = intersect_sphere(MIN_MU_S, ATMOSPHERE_INNER_RADIUS, ATMOSPHERE_OUTER_RADIUS);
    let D_sun = intersect_sun.y;
    let A_sun = (D_sun - ATMOSPHERE_THICKNESS) / (H - ATMOSPHERE_THICKNESS);
    let d_sun = intersect_sphere(mu_sun, ATMOSPHERE_INNER_RADIUS, ATMOSPHERE_OUTER_RADIUS).y;
    let a_sun = (d_sun - ATMOSPHERE_THICKNESS) / (H - ATMOSPHERE_THICKNESS);
    let u_mu_s_sun = max(1.0 - a_sun / A_sun, 0.0) / (1.0 + a_sun);
    let u_mu_s_sun_uv = get_uv_from_unit_range(u_mu_s_sun, SCATTERING_RES.z);

    // 月光散射计算 ---------------------------------------------------------
    // (结构相同，替换变量)
    let half_range_nu_moon = sqrt((1.0 - adjusted_mu * adjusted_mu) * (1.0 - mu_moon * mu_moon));
    let nu_min_moon = adjusted_mu * mu_moon - half_range_nu_moon;
    let nu_max_moon = adjusted_mu * mu_moon + half_range_nu_moon;
    let u_nu_moon = select(
        (nu_moon - nu_min_moon) / (nu_max_moon - nu_min_moon),
        nu_min_moon,
        nu_min_moon == nu_max_moon
    );
    let u_nu_moon_uv = get_uv_from_unit_range(u_nu_moon, SCATTERING_RES.x);

    // mu_s映射
    let intersect_moon = intersect_sphere(MIN_MU_S, ATMOSPHERE_INNER_RADIUS, ATMOSPHERE_OUTER_RADIUS);
    let D_moon = intersect_moon.y;
    let A_moon = (D_moon - ATMOSPHERE_THICKNESS) / (H - ATMOSPHERE_THICKNESS);
    let d_moon = intersect_sphere(mu_moon, ATMOSPHERE_INNER_RADIUS, ATMOSPHERE_OUTER_RADIUS).y;
    let a_moon = (d_moon - ATMOSPHERE_THICKNESS) / (H - ATMOSPHERE_THICKNESS);
    let u_mu_s_moon = max(1.0 - a_moon / A_moon, 0.0) / (1.0 + a_moon);
    let u_mu_s_moon_uv = get_uv_from_unit_range(u_mu_s_moon, SCATTERING_RES.z);

    //纹理采样 ------------------------------------------------------------
    let uv_sc_sun = vec3<f32>(u_nu_sun_uv * 0.5, u_mu, u_mu_s_sun_uv);
    let uv_sm_sun = vec3<f32>(u_nu_sun_uv * 0.5 + 0.5, u_mu, u_mu_s_sun_uv);
    let scattering_sc_sun = textureSample(atmosphere_scattering_lut, atmosphere_scattering_sampler, uv_sc_sun).rgb;
    let scattering_sm_sun = textureSample(atmosphere_scattering_lut, atmosphere_scattering_sampler, uv_sm_sun).rgb;

    let uv_sc_moon = vec3<f32>(u_nu_moon_uv * 0.5, u_mu, u_mu_s_moon_uv);
    let uv_sm_moon = vec3<f32>(u_nu_moon_uv * 0.5 + 0.5, u_mu, u_mu_s_moon_uv);
    let scattering_sc_moon = textureSample(atmosphere_scattering_lut, atmosphere_scattering_sampler, uv_sc_moon).rgb;
    let scattering_sm_moon = textureSample(atmosphere_scattering_lut, atmosphere_scattering_sampler, uv_sm_moon).rgb;

    // // 相位函数计算 ---------------------------------------------------------
    let mie_phase_sun = atmosphere_mie_phase(nu_sun);
    let mie_phase_moon = atmosphere_mie_phase(nu_moon);

    // 最终合成 ------------------------------------------------------------
    let sun_contribution = (scattering_sc_sun + scattering_sm_sun * mie_phase_sun) * sun_color;
    let moon_contribution = (scattering_sc_moon + scattering_sm_moon * mie_phase_moon) * moon_color;
    var atmosphere = sun_contribution + moon_contribution;

    return  atmosphere;
}

fn get_ambient_color(sun_color: vec3<f32>, moon_color: vec3<f32>) -> vec3<f32> {
    var sky_dir = normalize(vec3(0.0, 1.0, -0.8));
    var sky_color = atmosphere_scattering(sky_dir, sun_color, sun_dir, moon_color, moon_dir);
    sky_color = TAU * mix(sky_color,vec3(sky_color.b) * sqrt(2.0),1.0 / PI);
    return sky_color;
}

fn border_fog(scene_pos:vec3<f32>,world_dir:vec3<f32>)->f32 {
    var fog = cubic_length(scene_pos.xz) / TX;
    fog = exp2(-8.0 * pow(fog,8.0));
    fog = mix(fog, 1.0, 0.75 * dampen(linear_step(0.0, 0.2, world_dir.y)));
    return fog;
}

fn draw_sun(ray_dir:vec3<f32>,sun_color:vec3<f32>)->vec3<f32> {
	let nu = dot(ray_dir, sun_dir);
	let alpha = vec3(0.429, 0.522, 0.614);
	let center_to_edge = max(2 * (TAU / 360.0) - fast_acos(nu),0.0);
	let limb_darkening = pow(vec3(1.0 - sqr(1.0 - center_to_edge)), 0.5 * alpha);

	return 40.0 * sun_color * step(0.0, center_to_edge) * limb_darkening;
}
//ATOMOSPHERE

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

//FRAGMENT
const NUM_STEPS: i32 = 8;
const ITER_GEOMETRY: i32 = 4;
const ITER_FRAGMENT: i32 = 8;
const SEA_HEIGHT: f32 = 0.6;
const SEA_CHOPPY: f32 = 4.0;
const SEA_SPEED: f32 = 0.8;
const SEA_FREQ: f32 = 0.16;
const FOV:f32 = tan(radians(60.0));
const SEA_BASE: vec3<f32> = vec3<f32>(0.0, 0.09, 0.18);
const SEA_WATER_COLOR: vec3<f32> = vec3<f32>(0.8, 0.9, 0.6) * 0.6;
const octave_m: mat2x2<f32> = mat2x2<f32>(1.6, 1.2, -1.2, 1.6);

struct Uniforms {
    iTime: f32,
    @align(8) iResolution: vec2<f32>,
};

@group(0) @binding(0) var<uniform> uniforms: Uniforms;
@group(0) @binding(1) var atmosphere_scattering_lut: texture_3d<f32>;
@group(0) @binding(2) var atmosphere_scattering_sampler: sampler;

fn SEA_TIME() -> f32 { return 1.0 + uniforms.iTime * SEA_SPEED; }
fn EPSILON_NRM() -> f32 { return 0.05 / uniforms.iResolution.x; }
fn TIME() -> f32 {return uniforms.iTime * 0.1;}
fn ANGLE() -> vec3<f32> {return vec3<f32>(sin(TIME()*3.0)*0.1,sin(TIME())*0.2+0.3,TIME());}
fn POS() -> vec3<f32> {return vec3<f32>(0.0,5.0,TIME() * 5.0);}

fn fromEuler(ang: vec3<f32>) -> mat3x3<f32> {
    let a1 = vec2<f32>(sin(ang.x), cos(ang.x));
    let a2 = vec2<f32>(sin(ang.y), cos(ang.y));
    let a3 = vec2<f32>(sin(ang.z), cos(ang.z));
    return mat3x3<f32>(
        vec3<f32>(a1.y*a3.y + a1.x*a2.x*a3.x, a1.y*a2.x*a3.x + a3.y*a1.x, -a2.y*a3.x),
        vec3<f32>(-a2.y*a1.x, a1.y*a2.y, a2.x),
        vec3<f32>(a3.y*a1.x*a2.x + a1.y*a3.x, a1.x*a3.x - a1.y*a3.y*a2.x, a2.y*a3.y)
    );
}

fn hash(p: vec2<f32>) -> f32 {
    let h = dot(p, vec2<f32>(127.1, 311.7));    
    return fract(sin(h) * 43758.5453123);
}


fn noise(p: vec2<f32>) -> f32 {
    let i = floor(p);
    let f = fract(p);
    let u = f * f * (3.0 - 2.0 * f);
    
    return -1.0 + 2.0 * mix(
        mix(hash(i + vec2<f32>(0.0, 0.0)), hash(i + vec2<f32>(1.0, 0.0)), u.x),
        mix(hash(i + vec2<f32>(0.0, 1.0)), hash(i + vec2<f32>(1.0, 1.0)), u.x),
        u.y
    );
}

fn diffuse(n: vec3<f32>, l: vec3<f32>, p: f32) -> f32 {
    return pow(dot(n, l) * 0.4 + 0.6, p);
}

fn specular(n: vec3<f32>, l: vec3<f32>, e: vec3<f32>, s: f32) -> f32 {
    //let nrm = (s + 8.0) / (PI * 8.0); 
    return pow(max(dot(reflect(e, n), l), 0.0), s);
}

fn getSkyColor(e: vec3<f32>) -> vec3<f32> {
    let ey = (max(e.y, 0.0) * 0.8 + 0.2) * 0.8;
    return vec3<f32>(
        pow(1.0 - ey, 2.0),
        1.0 - ey,
        0.6 + (1.0 - ey) * 0.4
    ) * 1.1;
}

fn getSeaColor(p: vec3<f32>, n: vec3<f32>, l: vec3<f32>, eye: vec3<f32>, dist: vec3<f32>,
                sun_color:vec3<f32>,moon_color:vec3<f32>,
                ambient:vec3<f32>) -> vec3<f32> {  
    // 菲涅尔项计算
    var fresnel: f32 = clamp(1.0 - dot(n, -eye), 0.0, 1.0);
    fresnel = min(fresnel * fresnel, 0.5);
    
    // 反射和折射成分
    //let reflected: vec3<f32> = getSkyColor(reflect(eye, n));
    let reflected: vec3<f32> = atmosphere_scattering(reflect(eye, n),sun_color,sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0));
    let refracted: vec3<f32> = SEA_BASE  + (diffuse(n, l, 80.0) * sun_color * SEA_WATER_COLOR) * 0.12; 
    
    // 基础颜色混合
    var color: vec3<f32> = mix(refracted, reflected, fresnel);
    
    // 距离衰减效果
    let atten: f32 = max(1.0 - dot(dist, dist) * 0.001, 0.0);
    color += (SEA_WATER_COLOR * (p.y - SEA_HEIGHT)) * 0.18 * atten;
    
    // 高光添加
    color += ACESToneMapping(vec3(specular(n, l, eye, 60.0)) * sun_color,0.25);
    
    return color;
}

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
    var uv: vec2<f32> = p.xz;
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
    var uv: vec2<f32> = p.xz;
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

struct TracingResult {
    distance: f32,
    position: vec3<f32>,
    hit_sky: bool
};

fn heightMapTracing(ori: vec3<f32>, dir: vec3<f32>) -> TracingResult {
    var tm: f32 = 0.0;
    var tx: f32 = TX;
    var p: vec3<f32>;
    
    var hx = map(ori + dir * tx);
    if(hx > 0.0) {
        return TracingResult(tx, ori + dir * tx, true);
    }
    
    //var hm: f32 = map(ori);
    tx = min(POS().y / dir.y,tx);
    tm = max(POS().y / dir.y - 5.0,0.0);
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
        if(abs(hmid) < EPSILON) { break; }
    }
    return TracingResult(
        mix(tm, tx, hm / (hm - hx)),  
        p,
        false                             
    );
}

fn calculate_dir(uv: vec2<f32>) -> vec3<f32>{
    let ndc_uv = vec2<f32>(uv.x * 2.0 - 1.0,uv.y * 2.0 - 1.0);
    var dir = normalize(vec3<f32>(ndc_uv.x,ndc_uv.y,-FOV));
    //dir.z += length(dir.xy) * 0.14;
    dir = normalize(dir + vec3(0.0,0.2,0.0));
    //return normalize(dir);
     return normalize(dir) * fromEuler(ANGLE());
} 
//const pos:vec3<f32> = vec3<f32>(3.0,3.0,3.0);

fn get_pixel(uv:vec2<f32>,sun_color: vec3<f32>,moon_color: vec3<f32>,ambient:vec3<f32>)->vec3<f32>{
    let dir = calculate_dir(uv);
    var sky_color = atmosphere_scattering(dir,sun_color,sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0));
    sky_color += draw_sun(dir,sun_color);
    //if(dir.y > 0.1){return sky_color;}
    let res = heightMapTracing(POS(), dir);
    // if(res.hit_sky){ 
    //     return sky_color;
    // }

    let dist = res.position - POS();
    let n = getNormal(res.position,dot(dir, dist) * EPSILON_NRM());
    var sea_color = getSeaColor(res.position, n, sun_dir, dir, dist,sun_color,moon_color,ambient);
    let alp = 1 - smoothstep(-0.02,0.0,dir.y);

    let horizon_dir = normalize(vec3(dir.xz,min(dir.y,-0.1)).xzy);
    let horizon_color = atmosphere_scattering(horizon_dir,sun_color,sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0));
    let horizon_fac = linear_step(0.1,1.0,exp(-75 * sqr(sun_dir.y + 0.0496)));
    let fog_color = mix(sky_color,horizon_color,sqr(horizon_fac));
    let fog = border_fog(dist,dir);
    
    sea_color = mix(fog_color,sea_color,fog);
    sky_color = mix(fog_color,sky_color,fog);


    return mix(sky_color,sea_color,alp);
}


@fragment
fn fragment_main(@location(0) uv: vec2<f32>,@location(1) color: vec4<f32>,@location(2) sun_color:vec3<f32>,@location(3) moon_color:vec3<f32>,@location(4) ambient:vec3<f32>) -> @location(0) vec4<f32> {
    return vec4<f32>(pow(get_pixel(uv,sun_color,moon_color,ambient),vec3<f32>(0.65)),1.0);
}

//FRAGMENT