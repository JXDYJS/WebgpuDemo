const n=`/*\r
参考的代码和资料\r
https://www.shadertoy.com/view/Ms2SD1\r
https://zhuanlan.zhihu.com/p/127026136\r
https://zhuanlan.zhihu.com/p/477489052\r
https://www.gdcvault.com/play/1024478/PBR-Diffuse-Lighting-for-GGX\r
*/\r
\r
\r
\r
\r
struct Uniforms {\r
    iTime: f32,\r
    @align(8) iResolution: vec2<f32>,\r
    \r
};\r
@group(0) @binding(0) var<uniform> uniforms: Uniforms;\r
@group(0) @binding(1) var atmosphere_scattering_lut: texture_3d<f32>;\r
@group(0) @binding(2) var atmosphere_scattering_sampler: sampler;\r
\r
\r
//大气部分函数\r
const sun_dir = normalize(vec3<f32>(0.0, 0.2,1.0));\r
const moon_dir = normalize(vec3<f32>(0.0, -0.3,1.0));//目前不考虑夜晚，所以这行没用\r
const SUN_I: f32 = 1.0;\r
const MOON_I: f32 = 1.0;\r
const MOON_R: f32 = 0.75;\r
const MOON_G: f32 = 0.83;\r
const MOON_B: f32 = 1.0;\r
const PI: f32 = 3.141592;\r
const TAU: f32 = 6.283185307179586;\r
const HALF_PI: f32 = PI / 2.0;\r
const EPSILON: f32 = 1e-3;\r
const DEGREE: f32 = PI / 180.0;\r
const TX:f32 = 500.0;\r
\r
fn FADE()-> f32{\r
    if(sun_dir.y < 0.18){\r
        return 0.37 + 1.2 * max(0.0, -sun_dir.y);\r
    }\r
    return 0.17;\r
}\r
\r
fn SUN_WEIGHT()-> f32{\r
    return pow(clamp(1.0 - FADE() * abs(sun_dir.y - 0.18) ,0.0,1.0),2.0);\r
}\r
\r
fn SUN_RAISE()-> f32{\r
    return step(0.0001,sun_dir.x) * SUN_WEIGHT();\r
}\r
\r
fn SUN_SET()-> f32{\r
    return step(0.0001,-sun_dir.x) * SUN_WEIGHT();\r
}\r
\r
// 天体参数\r
const SUN_ANGULAR_RADIUS: f32 = 0.275 * DEGREE;\r
const MOON_ANGULAR_RADIUS: f32 = 0.278 * DEGREE;\r
const SUNLIGHT_COLOR: vec3<f32> = vec3(1.051, 0.985, 0.940);\r
const MOONLIGHT_COLOR: vec3<f32> = vec3(0.05, 0.1, 0.15);\r
\r
// 大气参数\r
const PLANET_RADIUS: f32 = 6371e3;\r
const ATMOSPHERE_INNER_RADIUS: f32 = PLANET_RADIUS - 1e3;\r
const ATMOSPHERE_OUTER_RADIUS: f32 = PLANET_RADIUS + 110e3;\r
const ATMOSPHERE_THICKNESS: f32 = ATMOSPHERE_OUTER_RADIUS - ATMOSPHERE_INNER_RADIUS;\r
\r
// LUT参数\r
const SCATTERING_RES: vec3<f32> = vec3(16.0, 64.0, 32.0);\r
const MIN_MU_S: f32 = -0.35;\r
\r
fn linear_step(edge0: f32, edge1: f32, x: f32) -> f32 {\r
    return clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);\r
}\r
\r
fn pulse(a: f32, b: f32, x: f32) -> f32 {\r
    return step(a, x) - step(b, x);\r
}\r
\r
fn from_srgb(srgb: vec3<f32>) -> vec3<f32> {\r
    let linear = select(\r
        srgb / 12.92,\r
        pow((srgb + 0.055) / 1.055, vec3<f32>(2.4)),\r
        srgb > vec3(0.04045)\r
    );\r
    return linear;\r
}\r
\r
fn smoothstep(low: f32, high: f32, x: f32) -> f32 {\r
    let t = clamp((x - low) / (high - low), 0.0, 1.0);\r
    return t * t * (3.0 - 2.0 * t);\r
}\r
\r
fn cube(a: f32) -> f32 {\r
    return a * a * a;\r
}\r
\r
fn sqr(a:f32)->f32{\r
    return a * a;\r
}\r
\r
fn pow5(x: f32) -> f32 {\r
    return x * x * x * x * x;\r
}\r
\r
fn fast_acos(x:f32) -> f32 {\r
	let C0 = 1.57018;\r
	let C1 = -0.201877;\r
	let C2 = 0.0464619;\r
\r
	var res = (C2 * abs(x) + C1) * abs(x) + C0; // p(x)\r
	res *= sqrt(1.0 - abs(x));\r
\r
	if(x>=0.0){\r
		return res;\r
    }\r
    else{\r
        return PI - res;\r
    }\r
}\r
\r
fn cubic_length(v:vec2<f32>)->f32 {\r
	return pow(cube(abs(v.x)) + cube(abs(v.y)), 0.333333);\r
}\r
\r
fn dampen(x:f32)->f32{\r
	let t = clamp(x,0.0,1.0);\r
	return t * (2.0 - t);\r
}\r
\r
fn ACESToneMapping(color: vec3<f32>, adapted_lum: f32) -> vec3<f32> {\r
    const A: f32 = 2.51;\r
    const B: f32 = 0.03;\r
    const C: f32 = 2.43;\r
    const D: f32 = 0.59;\r
    const E: f32 = 0.14;\r
\r
    let scaled_color = color * adapted_lum;\r
    return (scaled_color * (A * scaled_color + B)) / \r
           (scaled_color * (C * scaled_color + D) + E);\r
}\r
\r
fn get_uv_from_unit_range(val: f32, resolution: f32) -> f32 {\r
    return clamp(0.5 / resolution + val * (1.0 - 1.0 / resolution), 0.0, 1.0);\r
}\r
\r
fn intersect_sphere(mu: f32, radius: f32, outer_radius: f32) -> vec2<f32> {\r
    let discriminant = radius * radius * (mu * mu - 1.0) + outer_radius * outer_radius;\r
    let sqrt_d = sqrt(max(discriminant, 0.0));\r
    return vec2(-radius * mu - sqrt_d, -radius * mu + sqrt_d);\r
}\r
\r
fn atmosphere_mie_phase(nu: f32) -> f32 {\r
    const G: f32 = 0.76;\r
    let gg = G * G;\r
    return (3.0 * (1.0 - gg)) / (8.0 * PI * (2.0 + gg)) * (1.0 + nu * nu) / pow(1.0 + gg - 2.0 * G * nu, 1.5);\r
}\r
\r
fn get_sun_exposure() -> f32 {\r
    let time_sunset = SUN_SET();\r
    let time_sunrise = SUN_RAISE();\r
    let base_scale = 7.0 * SUN_I;\r
    \r
    // 蓝调时刻计算\r
    let blue_hour = linear_step(\r
        0.05, \r
        1.0, \r
        exp(-190.0 * pow(sun_dir.y + 0.09604, 2.0))\r
    );\r
    \r
    let daytime_mul = 1.0 + 0.5 * (time_sunset + time_sunrise) + 40.0 * blue_hour;\r
    \r
    return base_scale * daytime_mul;\r
}\r
\r
fn get_sun_tint() -> vec3<f32> {\r
    let time_sunset = SUN_SET();\r
    let time_sunrise = SUN_RAISE();\r
    // 蓝调时刻系数\r
    let blue_hour = linear_step(\r
        0.05, \r
        1.0, \r
        exp(-190.0 * pow(sun_dir.y + 0.09604, 2.0))\r
    );\r
    \r
    // 早晚色调\r
    var morning_evening_tint = vec3<f32>(1.05, 0.84, 0.93) * 1.2;\r
    morning_evening_tint = mix(\r
        vec3<f32>(1.0), \r
        morning_evening_tint, \r
        pow(pulse(0.17, 0.40, sun_dir.y), 2.0)\r
    );\r
    \r
    // 蓝调时刻色调\r
    var blue_hour_tint = vec3<f32>(0.95, 0.80, 1.0);\r
    blue_hour_tint = mix(vec3<f32>(1.0), blue_hour_tint, blue_hour);\r
    \r
    // 自定义色调(未来可能在设置内传入)\r
    let tint_morning = from_srgb(vec3<f32>(1.0, 1.0, 1.0));\r
    let tint_noon    = from_srgb(vec3<f32>(1.0, 1.0, 1.0));\r
    let tint_evening = from_srgb(vec3<f32>(1.0, 1.0, 1.0));\r
    \r
    var user_tint = mix(tint_noon, tint_morning, time_sunrise);\r
    user_tint = mix(user_tint, tint_evening, time_sunset);\r
    \r
    return morning_evening_tint * blue_hour_tint * user_tint;\r
}\r
\r
// 月光计算 -------------------------------------------------\r
fn get_moon_exposure() -> f32 {\r
    return 0.66 * MOON_I * 1.0;\r
}\r
\r
fn get_moon_tint() -> vec3<f32> {\r
    return from_srgb(vec3<f32>(MOON_R, MOON_G, MOON_B));\r
}\r
\r
// 主散射函数 ---------------------------------------------------------------\r
fn atmosphere_scattering(\r
    ray_dir: vec3<f32>,\r
    sun_color: vec3<f32>,\r
    sun_dir: vec3<f32>,\r
    moon_color: vec3<f32>,\r
    moon_dir: vec3<f32>\r
) -> vec3<f32> {\r
    // 基础参数\r
    let mu = ray_dir.y;\r
    let nu_sun = dot(ray_dir, sun_dir);\r
    let nu_moon = dot(ray_dir, moon_dir);\r
    let mu_sun = sun_dir.y;\r
    let mu_moon = moon_dir.y;\r
\r
    // 地平线优化\r
    var adjusted_mu = mu;\r
    let sun_step = smoothstep(-0.05, 0.1, mu_sun);\r
    let moon_step = smoothstep(0.05, 0.1, mu_moon);\r
    let horizon_mu = mix(-0.01, 0.03, clamp(sun_step + moon_step, 0.0, 1.0));\r
    adjusted_mu = max(mu, horizon_mu);\r
\r
    // 太阳散射计算 ---------------------------------------------------------\r
    // nu映射\r
    let half_range_nu_sun = sqrt((1.0 - adjusted_mu * adjusted_mu) * (1.0 - mu_sun * mu_sun));\r
    let nu_min_sun = adjusted_mu * mu_sun - half_range_nu_sun;\r
    let nu_max_sun = adjusted_mu * mu_sun + half_range_nu_sun;\r
    let u_nu_sun = select(\r
        (nu_sun - nu_min_sun) / (nu_max_sun - nu_min_sun),\r
        nu_min_sun,\r
        nu_min_sun == nu_max_sun\r
    );\r
    let u_nu_sun_uv = get_uv_from_unit_range(u_nu_sun, SCATTERING_RES.x);\r
\r
    // mu映射\r
    let r = PLANET_RADIUS;\r
    let H = sqrt(ATMOSPHERE_OUTER_RADIUS * ATMOSPHERE_OUTER_RADIUS - ATMOSPHERE_INNER_RADIUS * ATMOSPHERE_INNER_RADIUS);\r
    let rho = sqrt(max(r * r - ATMOSPHERE_INNER_RADIUS * ATMOSPHERE_INNER_RADIUS, 0.0));\r
    let rmu = r * adjusted_mu;\r
    let discriminant = rmu * rmu - r * r + ATMOSPHERE_INNER_RADIUS * ATMOSPHERE_INNER_RADIUS;\r
\r
    var u_mu: f32;\r
    if (adjusted_mu < 0.0 && discriminant >= 0.0) {\r
        let d = -rmu - sqrt(discriminant);\r
        let d_min = r - ATMOSPHERE_INNER_RADIUS;\r
        let d_max = rho;\r
        u_mu = (d - d_min) / (d_max - d_min);\r
        u_mu = get_uv_from_unit_range(u_mu, SCATTERING_RES.y / 2.0);\r
        u_mu = 0.5 - 0.5 * u_mu;\r
    } else {\r
        let d = -rmu + sqrt(discriminant + H * H);\r
        let d_min = ATMOSPHERE_OUTER_RADIUS - r;\r
        let d_max = rho + H;\r
        u_mu = (d - d_min) / (d_max - d_min);\r
        u_mu = get_uv_from_unit_range(u_mu, SCATTERING_RES.y / 2.0);\r
        u_mu = 0.5 + 0.5 * u_mu;\r
    }\r
\r
    // mu_s映射\r
    let intersect_sun = intersect_sphere(MIN_MU_S, ATMOSPHERE_INNER_RADIUS, ATMOSPHERE_OUTER_RADIUS);\r
    let D_sun = intersect_sun.y;\r
    let A_sun = (D_sun - ATMOSPHERE_THICKNESS) / (H - ATMOSPHERE_THICKNESS);\r
    let d_sun = intersect_sphere(mu_sun, ATMOSPHERE_INNER_RADIUS, ATMOSPHERE_OUTER_RADIUS).y;\r
    let a_sun = (d_sun - ATMOSPHERE_THICKNESS) / (H - ATMOSPHERE_THICKNESS);\r
    let u_mu_s_sun = max(1.0 - a_sun / A_sun, 0.0) / (1.0 + a_sun);\r
    let u_mu_s_sun_uv = get_uv_from_unit_range(u_mu_s_sun, SCATTERING_RES.z);\r
\r
    // 月光散射计算 ---------------------------------------------------------\r
    // (结构相同，替换变量)\r
    let half_range_nu_moon = sqrt((1.0 - adjusted_mu * adjusted_mu) * (1.0 - mu_moon * mu_moon));\r
    let nu_min_moon = adjusted_mu * mu_moon - half_range_nu_moon;\r
    let nu_max_moon = adjusted_mu * mu_moon + half_range_nu_moon;\r
    let u_nu_moon = select(\r
        (nu_moon - nu_min_moon) / (nu_max_moon - nu_min_moon),\r
        nu_min_moon,\r
        nu_min_moon == nu_max_moon\r
    );\r
    let u_nu_moon_uv = get_uv_from_unit_range(u_nu_moon, SCATTERING_RES.x);\r
\r
    // mu_s映射\r
    let intersect_moon = intersect_sphere(MIN_MU_S, ATMOSPHERE_INNER_RADIUS, ATMOSPHERE_OUTER_RADIUS);\r
    let D_moon = intersect_moon.y;\r
    let A_moon = (D_moon - ATMOSPHERE_THICKNESS) / (H - ATMOSPHERE_THICKNESS);\r
    let d_moon = intersect_sphere(mu_moon, ATMOSPHERE_INNER_RADIUS, ATMOSPHERE_OUTER_RADIUS).y;\r
    let a_moon = (d_moon - ATMOSPHERE_THICKNESS) / (H - ATMOSPHERE_THICKNESS);\r
    let u_mu_s_moon = max(1.0 - a_moon / A_moon, 0.0) / (1.0 + a_moon);\r
    let u_mu_s_moon_uv = get_uv_from_unit_range(u_mu_s_moon, SCATTERING_RES.z);\r
\r
    //纹理采样 ------------------------------------------------------------\r
    let uv_sc_sun = vec3<f32>(u_nu_sun_uv * 0.5, u_mu, u_mu_s_sun_uv);\r
    let uv_sm_sun = vec3<f32>(u_nu_sun_uv * 0.5 + 0.5, u_mu, u_mu_s_sun_uv);\r
    let scattering_sc_sun = textureSample(atmosphere_scattering_lut, atmosphere_scattering_sampler, uv_sc_sun).rgb;\r
    let scattering_sm_sun = textureSample(atmosphere_scattering_lut, atmosphere_scattering_sampler, uv_sm_sun).rgb;\r
\r
    let uv_sc_moon = vec3<f32>(u_nu_moon_uv * 0.5, u_mu, u_mu_s_moon_uv);\r
    let uv_sm_moon = vec3<f32>(u_nu_moon_uv * 0.5 + 0.5, u_mu, u_mu_s_moon_uv);\r
    let scattering_sc_moon = textureSample(atmosphere_scattering_lut, atmosphere_scattering_sampler, uv_sc_moon).rgb;\r
    let scattering_sm_moon = textureSample(atmosphere_scattering_lut, atmosphere_scattering_sampler, uv_sm_moon).rgb;\r
\r
    // // 相位函数计算 ---------------------------------------------------------\r
    let mie_phase_sun = atmosphere_mie_phase(nu_sun);\r
    let mie_phase_moon = atmosphere_mie_phase(nu_moon);\r
\r
    // 最终合成 ------------------------------------------------------------\r
    let sun_contribution = (scattering_sc_sun + scattering_sm_sun * mie_phase_sun) * sun_color;\r
    let moon_contribution = (scattering_sc_moon + scattering_sm_moon * mie_phase_moon) * moon_color;\r
    var atmosphere = sun_contribution + moon_contribution;\r
\r
    return  atmosphere;\r
}\r
\r
fn get_ambient_color(sun_color: vec3<f32>, moon_color: vec3<f32>) -> vec3<f32> {\r
    var sky_dir = normalize(vec3(0.0, 1.0, -0.8));\r
    var sky_color = atmosphere_scattering(sky_dir, sun_color, sun_dir, moon_color, moon_dir);\r
    sky_color = TAU * mix(sky_color,vec3(sky_color.b) * sqrt(2.0),1.0 / PI);\r
    return sky_color;\r
}\r
\r
fn border_fog(scene_pos:vec3<f32>,world_dir:vec3<f32>)->f32 {\r
    var fog = cubic_length(scene_pos.xz) / TX;\r
    fog = exp2(-8.0 * pow(fog,8.0));\r
    fog = mix(fog, 1.0, 0.75 * dampen(linear_step(0.0, 0.2, world_dir.y)));\r
    return fog;\r
}\r
\r
fn draw_sun(ray_dir:vec3<f32>,sun_color:vec3<f32>)->vec3<f32> {\r
	let nu = dot(ray_dir, sun_dir);\r
	let alpha = vec3(0.429, 0.522, 0.614);\r
	let center_to_edge = max(2 * (TAU / 360.0) - fast_acos(nu),0.0);\r
	let limb_darkening = pow(vec3(1.0 - sqr(1.0 - center_to_edge)), 0.5 * alpha);\r
\r
	return 40.0 * sun_color * step(0.0, center_to_edge) * limb_darkening;\r
}\r
//ATOMOSPHERE\r
\r
//\r
//PBR部分\r
//\r
\r
\r
// GGX法线分布函数\r
fn distribution_ggx(NoH: f32, alpha_sq: f32) -> f32 {\r
    let denom = (1.0 - NoH * NoH) + NoH * NoH * alpha_sq;\r
    return alpha_sq / (PI * denom * denom);\r
}\r
\r
// Smith单方向遮挡函数\r
fn visibility_smith_ggx_single(cos_theta: f32, alpha_sq: f32) -> f32 {\r
    let term1 = (-cos_theta * alpha_sq + cos_theta) * cos_theta;\r
    return 1.0 / (cos_theta + sqrt(term1 + alpha_sq));\r
}\r
\r
// Smith双方向遮挡函数（高度相关）\r
// 正确的高度相关Smith-Joint遮蔽函数\r
fn visibility_smith_ggx_joint(NoL: f32, NoV: f32, alpha_sq: f32) -> f32 {\r
    let lambda_v = NoV * sqrt( (-NoL * alpha_sq + NoL) * NoL + alpha_sq );\r
    let lambda_l = NoL * sqrt( (-NoV * alpha_sq + NoV) * NoV + alpha_sq );\r
    return 1.0 / (1.0 + lambda_v + lambda_l + 1e-5); // 添加防除零保护\r
}\r
\r
// 菲涅尔近似（Schlick）\r
fn fresnel_schlick(cos_theta: f32, f0: vec3<f32>) -> vec3<f32> {\r
    return f0 + (vec3<f32>(1.0) - f0) * pow5(1.0 - cos_theta);\r
}\r
\r
// 介质菲涅尔计算\r
fn fresnel_dielectric(cos_theta: f32, ior: f32) -> f32 {\r
    let g_sq = sqr(ior) + sqr(cos_theta) - 1.0;\r
    if (g_sq < 0.0){return 1.0;}\r
    let g = sqrt(g_sq);\r
    let a = (g - cos_theta) / max((g + cos_theta),EPSILON);\r
    let b = (cos_theta * (g + cos_theta) - 1.0) / max((cos_theta * (g - cos_theta) + 1.0),EPSILON);\r
    return 0.5 * a * a * (1.0 + b * b);\r
}\r
\r
// fn diffuse_hammon(\r
//     albedo: vec3<f32>,\r
//     roughness: f32,\r
//     f0: f32,\r
//     NoL: f32,\r
//     NoV: f32,\r
//     VoH: f32,\r
//     LoV: f32\r
// ) -> vec3<f32> {\r
//     // 参数预处理\r
//     let alpha = roughness * roughness;\r
//     let alpha_sq = alpha * alpha;\r
//     let sqrt_f0 = sqrt(f0) * 0.99999;\r
//     let ior = (1 + sqrt_f0) / (1 - sqrt_f0);\r
    \r
//     // 能量守恒因子\r
//     let energy_factor = 1.0 - (4.0 * sqrt(f0) + 5.0 * f0 * f0) / 9.0;\r
    \r
//     // 单次散射项\r
//     let facing = 0.5 * LoV + 0.5;\r
//     let fresnel_l = 1.0 - fresnel_dielectric(max(NoL, 0.01), ior);\r
//     let fresnel_v = 1.0 - fresnel_dielectric(max(NoV, 0.01), ior);\r
//     let single_smooth = (fresnel_l * fresnel_v) / max(energy_factor,EPSILON);\r
    \r
//     // 粗糙表面修正项（论文启发式公式）\r
//     let single_rough = max(facing, 0.0) * (-0.2 * facing + 0.45) * (1.0 / VoH + 2.0);\r
    \r
//     // 混合单次散射项\r
//     let single = mix(single_smooth, single_rough, roughness) / PI;\r
    \r
//     // 多次散射项（论文数值拟合）\r
//     let multi = 0.1159 * roughness;\r
    \r
//     // 最终组合\r
//     return albedo * (multi +single);\r
// }\r
\r
fn diffuse_water(\r
    albedo: vec3<f32>,\r
    roughness: f32,\r
    f0: f32,\r
    L: vec3<f32>,\r
    V: vec3<f32>,\r
    H: vec3<f32>\r
) -> vec3<f32> {\r
    // 微表面法线相关\r
    let LoH = max(dot(L, H), 0.0);\r
    let VoH = max(dot(V, H), 0.0);\r
    let sqrt_f0 = sqrt(f0) * 0.99999;\r
    let ior = (1 + sqrt_f0) / (1 - sqrt_f0);\r
    \r
    // 双向菲涅尔透射\r
    let F_in = fresnel_dielectric(LoH,ior);\r
    let F_out = fresnel_dielectric(VoH,ior);\r
    let transmittance = (1.0 - F_in) * (1.0 - F_out);\r
    \r
    // 补偿因子（适应水的IOR）\r
    let compensation = 1.0 / (1.0 - 0.28 * roughness);\r
    \r
    // 单次散射项改造\r
    let single = transmittance * compensation / PI;\r
    \r
    // 多次散射项（水体的吸收效应）\r
    let absorption = exp(-albedo * (1.0 - roughness) * 2.0);\r
    let multi = 0.1159 * roughness * absorption;\r
    \r
    return (single + multi) * albedo;\r
}\r
\r
struct pbr_shading_res{\r
    diffuse: vec3<f32>,\r
    specular: vec3<f32>\r
};\r
\r
fn pbr_shading(\r
    albedo: vec3<f32>,\r
    roughness: f32,\r
    f0:f32,\r
    metallic: f32,\r
    N: vec3<f32>,\r
    L: vec3<f32>,\r
    V: vec3<f32>,\r
    H: vec3<f32>\r
) -> pbr_shading_res {\r
    // 几何参数\r
    var res: pbr_shading_res;\r
    let NoL = max(dot(N, L), 0.0);\r
    let NoV = max(dot(N, V), 0.0);\r
    let NoH = max(dot(N, H), 0.0);\r
    let VoH = max(dot(V, H), 0.0);\r
    let LoV = max(dot(L, V), 0.0);\r
    \r
    // 高光项计算\r
    let alpha_sq = (roughness * roughness);\r
    let D = distribution_ggx(NoH, alpha_sq);\r
    let G = visibility_smith_ggx_joint(NoL, NoV, alpha_sq);\r
    let F = fresnel_schlick(VoH, mix(vec3<f32>(0.02), albedo, metallic));\r
    \r
    // 高光BRDF\r
    let specular_brdf = (D * G * F) / (4.0 * NoL * NoV + 1e-5);\r
    \r
    // 漫反射项\r
    let diffuse_brdf = diffuse_water(albedo, roughness, f0, L,V,H);\r
    \r
    let kD = (vec3<f32>(1.0) - F) * (1.0 - metallic);\r
    res.diffuse = diffuse_brdf * kD * NoL;\r
    res.specular = specular_brdf * NoL;\r
    return res;\r
}\r
\r
//PBR\r
\r
\r
//VERTEX\r
\r
struct VertexOutput {\r
    @builtin(position) position: vec4<f32>,\r
    @location(0) uv: vec2<f32>,\r
    @location(1) color: vec4<f32>,\r
    @location(2) sun_color:vec3<f32>,\r
    @location(3) moon_color:vec3<f32>,\r
    @location(4) ambient:vec3<f32>\r
};\r
\r
@vertex\r
fn vertex_main(@location(0) position: vec2<f32>,@location(1) color:vec4<f32>) -> VertexOutput{\r
    var output: VertexOutput;\r
    output.position = vec4<f32>(position * 2.0 - 1.0,0.0,1.0);\r
    output.uv = vec2<f32>(position.x,position.y);\r
    output.color = color;\r
    output.sun_color = get_sun_exposure() * get_sun_tint();\r
    output.moon_color = get_moon_exposure() * get_moon_tint();\r
    output.ambient = vec3(0.15,0.2,0.1);\r
    return output;\r
}\r
\r
//VERTEX\r
\r
//FRAGMENT\r
const NUM_STEPS: i32 = 8;\r
const ITER_GEOMETRY: i32 = 4;\r
const ITER_FRAGMENT: i32 = 8;\r
const SEA_HEIGHT: f32 = 0.6;\r
const SEA_CHOPPY: f32 = 4.0;\r
const SEA_SPEED: f32 = 0.8;\r
const SEA_FREQ: f32 = 0.16;\r
const FOV:f32 = tan(radians(60.0));\r
const SEA_BASE: vec3<f32> = vec3<f32>(0.0, 0.09, 0.18);\r
const SEA_WATER_COLOR: vec3<f32> = vec3<f32>(0.7, 0.8, 0.9) * 0.6;\r
const octave_m: mat2x2<f32> = mat2x2<f32>(1.6, 1.2, -1.2, 1.6);\r
\r
fn SEA_TIME() -> f32 { return 1.0 + uniforms.iTime * SEA_SPEED; }\r
fn EPSILON_NRM() -> f32 { return 0.05 / uniforms.iResolution.x; }\r
fn TIME() -> f32 {return uniforms.iTime * 0.1;}\r
fn ANGLE() -> vec3<f32> {return vec3<f32>(sin(TIME()*3.0)*0.1,sin(TIME())*0.2+0.3,TIME());}\r
fn POS() -> vec3<f32> {return vec3<f32>(0.0,5.0,TIME() * 5.0);}\r
\r
fn fromEuler(ang: vec3<f32>) -> mat3x3<f32> {\r
    let a1 = vec2<f32>(sin(ang.x), cos(ang.x));\r
    let a2 = vec2<f32>(sin(ang.y), cos(ang.y));\r
    let a3 = vec2<f32>(sin(ang.z), cos(ang.z));\r
    return mat3x3<f32>(\r
        vec3<f32>(a1.y*a3.y + a1.x*a2.x*a3.x, a1.y*a2.x*a3.x + a3.y*a1.x, -a2.y*a3.x),\r
        vec3<f32>(-a2.y*a1.x, a1.y*a2.y, a2.x),\r
        vec3<f32>(a3.y*a1.x*a2.x + a1.y*a3.x, a1.x*a3.x - a1.y*a3.y*a2.x, a2.y*a3.y)\r
    );\r
}\r
\r
fn hash(p: vec2<f32>) -> f32 {\r
    let h = dot(p, vec2<f32>(127.1, 311.7));    \r
    return fract(sin(h) * 43758.5453123);\r
}\r
\r
\r
fn noise(p: vec2<f32>) -> f32 {\r
    let i = floor(p);\r
    let f = fract(p);\r
    let u = f * f * (3.0 - 2.0 * f);\r
    \r
    return -1.0 + 2.0 * mix(\r
        mix(hash(i + vec2<f32>(0.0, 0.0)), hash(i + vec2<f32>(1.0, 0.0)), u.x),\r
        mix(hash(i + vec2<f32>(0.0, 1.0)), hash(i + vec2<f32>(1.0, 1.0)), u.x),\r
        u.y\r
    );\r
}\r
\r
fn diffuse(n: vec3<f32>, l: vec3<f32>, p: f32) -> f32 {\r
    return pow(dot(n, l) * 0.4 + 0.6, p);\r
}\r
\r
fn specular(n: vec3<f32>, l: vec3<f32>, e: vec3<f32>, s: f32) -> f32 {\r
    //let nrm = (s + 8.0) / (PI * 8.0); \r
    return pow(max(dot(reflect(e, n), l), 0.0), s);\r
}\r
\r
fn getSkyColor(e: vec3<f32>) -> vec3<f32> {\r
    let ey = (max(e.y, 0.0) * 0.8 + 0.2) * 0.8;\r
    return vec3<f32>(\r
        pow(1.0 - ey, 2.0),\r
        1.0 - ey,\r
        0.6 + (1.0 - ey) * 0.4\r
    ) * 1.1;\r
}\r
\r
fn getSeaColor(p: vec3<f32>, n: vec3<f32>, l: vec3<f32>, eye: vec3<f32>, dist: vec3<f32>,\r
                sun_color:vec3<f32>,moon_color:vec3<f32>,\r
                ambient:vec3<f32>) -> vec3<f32> {  \r
    // // //菲涅尔项计算\r
    // var fresnel: f32 = clamp(1.0 - dot(n, -eye), 0.0, 1.0);\r
    // fresnel = min(fresnel * fresnel * fresnel, 0.5);\r
    \r
    //反射和折射成分\r
    //let reflected: vec3<f32> = getSkyColor(reflect(eye, n));\r
    // let rl = reflect(eye, n);\r
    // let h = normalize(-eye + l);\r
    // let pbr_res = pbr_shading(SEA_WATER_COLOR * sun_color * 0.12,0.002,0.02,0.0,n,l,-eye,h);\r
    // let reflected: vec3<f32> = atmosphere_scattering(rl,sun_color,sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0)) + draw_sun(rl,sun_color);\r
    // let refracted: vec3<f32> = SEA_BASE  + (diffuse(n, l, 80.0) * sun_color * SEA_WATER_COLOR) * 0.12; \r
    \r
    // // 基础颜色混合\r
    // var color: vec3<f32> = mix(refracted, reflected, fresnel);\r
    \r
    // // 距离衰减效果\r
    // let atten: f32 = max(1.0 - dot(dist, dist) * 0.001, 0.0);\r
    // color += (SEA_WATER_COLOR * (p.y - SEA_HEIGHT)) * 0.18 * atten;\r
    \r
    // //高光添加\r
    // color += ACESToneMapping(vec3(specular(n, l, eye, 60.0)) * sun_color,0.25);\r
    \r
    //return ACESToneMapping(refracted,0.72);\r
    //\r
    //\r
    //\r
\r
    let fresnel = fresnel_schlick(dot(n, -eye), vec3(0.02));\r
    \r
    //反射和折射成分\r
    //let reflected: vec3<f32> = getSkyColor(reflect(eye, n));\r
    let rl = reflect(eye, n);\r
    let h = normalize(-eye + l);\r
    let pbr_res = pbr_shading(SEA_WATER_COLOR * sun_color * 0.12,0.002,0.02,0.0,n,l,-eye,h);\r
    let reflected: vec3<f32> = atmosphere_scattering(rl,sun_color,sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0)) + draw_sun(rl,sun_color);\r
    let refracted: vec3<f32> = SEA_BASE  + (pbr_res.diffuse * sun_color) * 0.05; \r
    \r
    // 基础颜色混合\r
    var color: vec3<f32> = mix(refracted, reflected, fresnel);\r
    \r
    // 距离衰减效果\r
    let atten: f32 = max(1.0 - dot(dist, dist) * 0.001, 0.0);\r
    color += (SEA_WATER_COLOR * (p.y - SEA_HEIGHT)) * 0.18 * atten;\r
    \r
    //高光添加\r
    color += ACESToneMapping(pbr_res.specular * sun_color,0.25);\r
    \r
    return ACESToneMapping(color,0.72);\r
}\r
\r
fn sea_octave(uv: vec2<f32>, choppy: f32) -> f32 {\r
    var uv_mut = uv;\r
    uv_mut += noise(uv);\r
    let wv = 1.0 - abs(sin(uv_mut));\r
    let swv = abs(cos(uv_mut));\r
    return pow(1.0 - pow(mix(wv, swv, wv).x * mix(wv, swv, wv).y, 0.65), choppy);\r
}\r
\r
fn map(p: vec3<f32>) -> f32 {\r
    var freq: f32 = SEA_FREQ;\r
    var amp: f32 = SEA_HEIGHT;\r
    var choppy: f32 = SEA_CHOPPY;\r
    var uv: vec2<f32> = p.xz;\r
    uv.x *= 0.75;\r
    \r
    var h: f32 = 0.0;\r
    for(var i: i32 = 0; i < ITER_GEOMETRY; i++) { \r
        let d = sea_octave((uv + SEA_TIME()) * freq, choppy) +\r
                sea_octave((uv - SEA_TIME()) * freq, choppy);\r
        h += d * amp;\r
        uv = uv * octave_m;\r
        freq *= 1.9;\r
        amp *= 0.22;\r
        choppy = mix(choppy, 1.0, 0.2);\r
    }\r
    return p.y - h;\r
}\r
\r
fn map_detail(p: vec3<f32>) -> f32 {\r
    var freq: f32 = SEA_FREQ;\r
    var amp: f32 = SEA_HEIGHT;\r
    var choppy: f32 = SEA_CHOPPY;\r
    var uv: vec2<f32> = p.xz;\r
    uv.x *= 0.75;\r
    \r
    var h: f32 = 0.0;\r
    for(var i: i32 = 0; i < ITER_FRAGMENT; i++) { \r
        let d = sea_octave((uv + SEA_TIME()) * freq, choppy) +\r
                sea_octave((uv - SEA_TIME()) * freq, choppy);\r
        h += d * amp;\r
        uv = uv * octave_m;\r
        freq *= 1.9;\r
        amp *= 0.22;\r
        choppy = mix(choppy, 1.0, 0.2);\r
    }\r
    return p.y - h;\r
}\r
\r
fn getNormal(p:vec3<f32>,eps:f32) -> vec3<f32> {\r
    var n:vec3<f32>;\r
    n.y = map_detail(p);\r
    n.x = map_detail(vec3<f32>(p.x + eps, p.y, p.z)) - n.y;\r
    n.z = map_detail(vec3<f32>(p.x, p.y, p.z + eps)) - n.y;\r
    n.y = eps;\r
    return normalize(n);\r
}\r
\r
struct TracingResult {\r
    distance: f32,\r
    position: vec3<f32>,\r
    hit_sky: bool\r
};\r
\r
fn heightMapTracing(ori: vec3<f32>, dir: vec3<f32>) -> TracingResult {\r
    var tm: f32 = 0.0;\r
    var tx: f32 = TX;\r
    var p: vec3<f32>;\r
    \r
    var hx = map(ori + dir * tx);\r
    if(hx > 0.0) {\r
        return TracingResult(tx, ori + dir * tx, true);\r
    }\r
    \r
    //var hm: f32 = map(ori);\r
    tx = min(POS().y / dir.y,tx);\r
    tm = max(POS().y / dir.y - 5.0,0.0);\r
    var hm = map(ori + dir * tm);\r
    hx = map(ori + dir * tx);\r
    for(var i: i32 = 0; i < NUM_STEPS; i++) {\r
        let tmid: f32 = mix(tm, tx, hm / (hm - hx));\r
        p = ori + dir * tmid;\r
        let hmid = map(p);\r
        \r
        if(hmid < 0.0) {\r
            tx = tmid;\r
            hx = hmid;\r
        } else {\r
            tm = tmid;\r
            hm = hmid;\r
        }\r
        if(abs(hmid) < EPSILON) { break; }\r
    }\r
    return TracingResult(\r
        mix(tm, tx, hm / (hm - hx)),  \r
        p,\r
        false                             \r
    );\r
}\r
\r
fn calculate_dir(uv: vec2<f32>) -> vec3<f32>{\r
    let ndc_uv = vec2<f32>(uv.x * 2.0 - 1.0,uv.y * 2.0 - 1.0);\r
    var dir = normalize(vec3<f32>(ndc_uv.x,ndc_uv.y,-FOV));\r
    //dir.z += length(dir.xy) * 0.14;\r
    //dir = normalize(dir + vec3(0.0,0.2,0.0));\r
    //return normalize(dir);\r
     return normalize(dir) * fromEuler(ANGLE());\r
} \r
//const pos:vec3<f32> = vec3<f32>(3.0,3.0,3.0);\r
\r
fn get_pixel(uv:vec2<f32>,sun_color: vec3<f32>,moon_color: vec3<f32>,ambient:vec3<f32>)->vec3<f32>{\r
    let dir = calculate_dir(uv);\r
    var sky_color = atmosphere_scattering(dir,sun_color,sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0));\r
    sky_color += draw_sun(dir,sun_color);\r
    sky_color = ACESToneMapping(sky_color,0.5);\r
    let res = heightMapTracing(POS(), dir);\r
\r
    let dist = res.position - POS();\r
    let n = getNormal(res.position,dot(dir, dist) * EPSILON_NRM());\r
    var sea_color = getSeaColor(res.position, n, sun_dir, dir, dist,sun_color,moon_color,ambient);\r
    let alp = 1 - smoothstep(-0.02,0.0,dir.y);\r
\r
    let horizon_dir = normalize(vec3(dir.xz,min(dir.y,-0.1)).xzy);\r
    let horizon_color = atmosphere_scattering(horizon_dir,sun_color,sun_dir,vec3(0.0,0.0,0.0),vec3(0.0,1.0,0.0));\r
    let horizon_fac = linear_step(0.1,1.0,exp(-75 * sqr(sun_dir.y + 0.0496)));\r
    let fog_color = mix(sky_color,horizon_color,sqr(horizon_fac));\r
    let fog = border_fog(dist,dir);\r
    \r
    sea_color = mix(fog_color,sea_color,fog);\r
    sky_color = mix(fog_color,sky_color,fog);\r
\r
\r
    return mix(sky_color,sea_color,alp);\r
}\r
\r
\r
@fragment\r
fn fragment_main(@location(0) uv: vec2<f32>,@location(1) color: vec4<f32>,@location(2) sun_color:vec3<f32>,@location(3) moon_color:vec3<f32>,@location(4) ambient:vec3<f32>) -> @location(0) vec4<f32> {\r
    return vec4<f32>(pow(get_pixel(uv,sun_color,moon_color,ambient),vec3<f32>(0.65)),1.0);\r
}\r
\r
//FRAGMENT`;export{n as default};
