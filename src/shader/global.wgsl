#ifndef GLOBAL_WGSL
#define GLOBAL_WGSL

#define HighQualitySea

struct Uniforms {
    iTime: f32,
    @align(8) iResolution: vec2<f32>,
    @align(16) sun_dir: vec3<f32>,
    @align(32) ball_data: vec2<f32>,
};

struct Pos_data {
    ball_pos: vec3<f32>,
    camera_pos: vec3<f32>,
    light_pos: vec3<f32>,
};

struct Material {
    albedo: vec3<f32>,
    f0: vec3<f32>,
    roughness: f32,
    metallic: f32,
};

struct pbr_shading_res {
    diffuse: vec3<f32>,
    specular: vec3<f32>
};

struct TracingResult {
    distance: f32,
    position: vec3<f32>,
    hit_type: u32, // 0 SEA 1 SKY 2 BALL
}; 

// 常量定义
const SUN_I: f32 = 1.0;
const MOON_I: f32 = 1.0;
const MOON_R: f32 = 0.75;
const MOON_G: f32 = 0.83;
const MOON_B: f32 = 1.0;
const PI: f32 = 3.141592;
const TAU: f32 = 6.283185307179586;
const HALF_PI: f32 = PI / 2.0;
const RCP_PI: f32 = 0.3183098861837907;
const EPSILON: f32 = 1e-3;
const DEGREE: f32 = PI / 180.0;
const TX:f32 = 500.0;

#if defined(HighQualitySea)
    const NUM_STEPS: i32 = 8;
#else
    const NUM_STEPS: i32 = 4;
#endif
const ITER_GEOMETRY: i32 = 4;
const ITER_FRAGMENT: i32 = 8;
const SEA_HEIGHT: f32 = 0.6;
const SEA_CHOPPY: f32 = 4.0;
const SEA_SPEED: f32 = 0.8;
const SEA_FREQ: f32 = 0.16;
const FOV:f32 = tan(radians(60.0));
const SEA_BASE: vec3<f32> = vec3<f32>(0.0, 0.09, 0.18);
const SEA_WATER_COLOR: vec3<f32> = vec3<f32>(0.7, 0.8, 0.9) * 0.6;
const octave_m: mat2x2<f32> = mat2x2<f32>(1.6, 1.2, -1.2, 1.6);
const ball_pos : vec3<f32> = vec3<f32>(0.0, 6.0, -10.0);
const ball_radius : f32 = 1.0;
const light_dist : f32 = 2.0;
const light_speed : f32 = 5.0;
const light_strength : f32 = 1.0;
const ball_material:Material = Material(vec3(0.8,0.5,0.5),vec3(0.04),0.1,0.0);
const env_sample_size:f32 = 128;
const KERNEL_5: array<f32, 5> = array<f32, 5>(
    0.0625, 
    0.25, 
    0.375, 
    0.25, 
    0.0625 
);

// 绑定资源
@group(0) @binding(0) var<uniform> uniforms: Uniforms;
@group(0) @binding(1) var atmosphere_scattering_lut: texture_3d<f32>;
@group(0) @binding(2) var atmosphere_scattering_sampler: sampler;
@group(0) @binding(3) var env_texture: texture_2d<f32>;
@group(0) @binding(4) var env_sampler: sampler;
@group(0) @binding(5) var env_last_frame: texture_2d<f32>;
@group(0) @binding(6) var env_last_frame_sampler: sampler;

// ================================
// 数学工具函数
// ================================
// fn linear_step(edge0: f32, edge1: f32, x: f32) -> f32 {
//     return clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
// }
#define linear_step(edge0, edge1, x) (clamp(((x) - (edge0)) / ((edge1) - (edge0)), 0.0, 1.0))

// fn pulse(a: f32, b: f32, x: f32) -> f32 {
//     return step(a, x) - step(b, x);
// }

#define pulse(a, b, x) (step((a), (x)) - step((b), (x)))

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

#define cube(x) ((x) * (x) * (x))

#define sqr(x) ((x) * (x))

#define pow5(x) ((x) * (x) * (x) * (x) * (x))

// fn max30(a: vec3<f32>)->vec3<f32>{
//     return vec3<f32>(max(a.x,0.0),max(a.y,0.0),max(a.z,0.0));
// }
#define max30(a) (vec3<f32>(max(a.x,0.0),max(a.y,0.0),max(a.z,0.0)))

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
    return pow(max(dot(reflect(e, n), l), 0.0), s);
}

fn f0_to_ior(f0: f32) -> f32 {
    return (1.0 + sqrt(f0)) / (1.0 - sqrt(f0));
}

// ================================
// 时间/方向相关辅助函数
// ================================
//fn SEA_TIME() -> f32 { return 1.0 + uniforms.iTime * SEA_SPEED; }
#define SEA_TIME() (1.0 + uniforms.iTime * SEA_SPEED)
fn EPSILON_NRM() -> f32 { return 0.1 / uniforms.iResolution.x; }
fn TIME() -> f32 {return uniforms.iTime * 0.1;}
fn ANGLE() -> vec3<f32> {return vec3<f32>(sin(TIME()*3.0)*0.1,0.0,TIME());}
fn POS() -> vec3<f32> {return vec3<f32>(0.0,5.0,0.0);}

fn FADE()-> f32{
    if(uniforms.sun_dir.y < 0.18){
        return 0.37 + 1.2 * max(0.0, -uniforms.sun_dir.y);
    }
    return 0.17;
}

fn SUN_WEIGHT()-> f32{
    return pow(clamp(1.0 - FADE() * abs(uniforms.sun_dir.y - 0.18) ,0.0,1.0),2.0);
}

fn SUN_RAISE()-> f32{
    return step(0.0001,uniforms.sun_dir.x) * SUN_WEIGHT();
}

fn SUN_SET()-> f32{
    return step(0.0001,-uniforms.sun_dir.x) * SUN_WEIGHT();
}
#endif