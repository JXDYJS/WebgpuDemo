const n=`
struct VertexOutput {
    @builtin(position) position: vec4<f32>,
    @location(0) uv: vec2<f32>,
    @location(1) color: vec4<f32>
};

@vertex
fn vertex_main(@location(0) position: vec2<f32>,@location(1) color:vec4<f32>) -> VertexOutput{
    var output: VertexOutput;
    output.position = vec4<f32>(position * 2.0 - 1.0,0.0,1.0);
    output.uv = vec2<f32>(position.x,position.y);
    output.color = color;
    return output;
}

const PI: f32 = 3.141592;
const EPSILON: f32 = 1e-3;
const NUM_STEPS: i32 = 8;
const ITER_GEOMETRY: i32 = 4;
const ITER_FRAGMENT: i32 = 6;
const SEA_HEIGHT: f32 = 0.6;
const SEA_CHOPPY: f32 = 4.0;
const SEA_SPEED: f32 = 0.8;
const SEA_FREQ: f32 = 0.16;
const SEA_BASE: vec3<f32> = vec3<f32>(0.0, 0.09, 0.18);
const SEA_WATER_COLOR: vec3<f32> = vec3<f32>(0.8, 0.9, 0.6) * 0.6;
const octave_m: mat2x2<f32> = mat2x2<f32>(1.6, 1.2, -1.2, 1.6);

struct Uniforms {
    iTime: f32,
    @align(8) iResolution: vec2<f32>,
};

@group(0) @binding(0) var<uniform> uniforms: Uniforms;

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
    let nrm = (s + 8.0) / (PI * 8.0); 
    return pow(max(dot(reflect(e, n), l), 0.0), s) * nrm;
}

fn getSkyColor(e: vec3<f32>) -> vec3<f32> {
    let ey = (max(e.y, 0.0) * 0.8 + 0.2) * 0.8;
    return vec3<f32>(
        pow(1.0 - ey, 2.0),
        1.0 - ey,
        0.6 + (1.0 - ey) * 0.4
    ) * 1.1;
}

fn getSeaColor(p: vec3<f32>, n: vec3<f32>, l: vec3<f32>, eye: vec3<f32>, dist: vec3<f32>) -> vec3<f32> {  
    // 菲涅尔项计算
    var fresnel: f32 = clamp(1.0 - dot(n, -eye), 0.0, 1.0);
    fresnel = min(fresnel * fresnel * fresnel, 0.5);
    
    // 反射和折射成分
    let reflected: vec3<f32> = getSkyColor(reflect(eye, n));    
    let refracted: vec3<f32> = SEA_BASE + (diffuse(n, l, 80.0) * SEA_WATER_COLOR) * 0.12; 
    
    // 基础颜色混合
    var color: vec3<f32> = mix(refracted, reflected, fresnel);
    
    // 距离衰减效果
    let atten: f32 = max(1.0 - dot(dist, dist) * 0.001, 0.0);
    color += (SEA_WATER_COLOR * (p.y - SEA_HEIGHT)) * 0.18 * atten;
    
    // 高光添加
    color += specular(n, l, eye, 60.0);
    
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
    var tx: f32 = 1000.0;
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
    var dir = normalize(vec3<f32>(ndc_uv.x,ndc_uv.y,-2.0));
    dir.z += length(dir.xy) * 0.14;
    //return normalize(dir);
     return normalize(dir) * fromEuler(ANGLE());
} 
//const pos:vec3<f32> = vec3<f32>(3.0,3.0,3.0);

fn get_pixel(uv:vec2<f32>)->vec3<f32>{
    let dir = calculate_dir(uv);
    if(dir.y > 0.1){return getSkyColor(dir);}
    let res = heightMapTracing(POS(), dir);
    let sky_color = getSkyColor(dir);
    if(res.hit_sky){ 
        return sky_color;
    }

    let dist = res.position - POS();
    let n = getNormal(res.position,dot(dir, dist) * EPSILON_NRM());
    let light = normalize(vec3<f32>(0.0, 1.0, 0.8));
    let sea_color = getSeaColor(res.position, n, light, dir, dist);
    let alp = 1 - smoothstep(-0.02,0.0,dir.y);
    return mix(sky_color,sea_color,alp);
}


@fragment
fn fragment_main(@location(0) uv: vec2<f32>,@location(1) color: vec4<f32>)-> @location(0) vec4<f32> {
    return vec4<f32>(pow(get_pixel(uv),vec3<f32>(0.65)),1.0);
}`;export{n as default};
