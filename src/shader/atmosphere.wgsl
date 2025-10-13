#ifndef ATMOSPHERE_WGSL
#define ATMOSPHERE_WGSL

#include "global.wgsl"

const sun_dir = normalize(vec3<f32>(1.0, 0.2,0.0));
const moon_dir = normalize(vec3<f32>(0.0, -0.3,1.0));
const SUN_ANGULAR_RADIUS: f32 = 0.275 * DEGREE;
const MOON_ANGULAR_RADIUS: f32 = 0.278 * DEGREE;
const SUNLIGHT_COLOR: vec3<f32> = vec3(1.051, 0.985, 0.940);
const MOONLIGHT_COLOR: vec3<f32> = vec3(0.05, 0.1, 0.15);
const PLANET_RADIUS: f32 = 6371e3;
const ATMOSPHERE_INNER_RADIUS: f32 = PLANET_RADIUS - 1e3;
const ATMOSPHERE_OUTER_RADIUS: f32 = PLANET_RADIUS + 110e3;
const ATMOSPHERE_THICKNESS: f32 = ATMOSPHERE_OUTER_RADIUS - ATMOSPHERE_INNER_RADIUS;
const SCATTERING_RES: vec3<f32> = vec3(16.0, 64.0, 32.0);
const MIN_MU_S: f32 = -0.35;

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
        exp(-190.0 * pow(uniforms.sun_dir.y + 0.09604, 2.0))
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
        exp(-190.0 * pow(uniforms.sun_dir.y + 0.09604, 2.0))
    );
    
    // 早晚色调
    var morning_evening_tint = vec3<f32>(1.05, 0.84, 0.93) * 1.2;
    morning_evening_tint = mix(
        vec3<f32>(1.0), 
        morning_evening_tint, 
        pow(pulse(0.17, 0.40, uniforms.sun_dir.y), 2.0)
    );
    
    // 蓝调时刻色调
    var blue_hour_tint = vec3<f32>(0.95, 0.80, 1.0);
    blue_hour_tint = mix(vec3<f32>(1.0), blue_hour_tint, blue_hour);
    
    // 自定义色调(未来可能在设置内传入)
    let tint_morning = from_srgb(vec3<f32>(1.0, 1.0, 1.0));
    let tint_noon    = from_srgb(vec3<f32>(1.0, 1.0, 1.0));
    let tint_evening = from_srgb(vec3<f32>(1.0, 1.0, 1.0));
    
    var user_tint = mix(tint_noon, tint_morning, time_sunrise);
    user_tint = mix(user_tint, tint_evening, time_sunset);
    
    return morning_evening_tint * blue_hour_tint * user_tint;
}

// 月光计算
fn get_moon_exposure() -> f32 {
    return 0.66 * MOON_I * 1.0;
}

fn get_moon_tint() -> vec3<f32> {
    return from_srgb(vec3<f32>(MOON_R, MOON_G, MOON_B));
}

// 主散射函数
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

    // 相位函数计算 ---------------------------------------------------------
    let mie_phase_sun = atmosphere_mie_phase(nu_sun);
    let mie_phase_moon = atmosphere_mie_phase(nu_moon);

    // 最终合成 ------------------------------------------------------------
    let sun_contribution = (scattering_sc_sun + scattering_sm_sun * mie_phase_sun) * sun_color;
    let moon_contribution = (scattering_sc_moon + scattering_sm_moon * mie_phase_moon) * moon_color;
    var atmosphere = sun_contribution + moon_contribution;

    return  atmosphere;
}

// ================================
// 天空/环境映射函数
// ================================
fn project_sky(direction: vec3<f32>) -> vec2<f32> {
    var projected_dir: vec2<f32> = normalize(direction.xz);
    
    let azimuth_angle: f32 = PI + atan2(projected_dir.x, -projected_dir.y);
    let altitude_angle: f32 = HALF_PI - fast_acos(clamp(direction.y, -1.0, 1.0));
    
    // 计算基础坐标
    var coord: vec2<f32>;
    coord.x = azimuth_angle * (1.0 / TAU);
    let signed_alt: f32 = sign(altitude_angle);
    coord.y = 0.5 + 0.5 * signed_alt * sqrt(2.0 * RCP_PI * abs(altitude_angle));
    
    // 边缘填充计算
    let sky_map_pixel_size: vec2<f32> = 1.0 / uniforms.iResolution;
    let pad_amount: f32 = 2.0;
    let mul: f32 = 1.0 - 2.0 * pad_amount * sky_map_pixel_size.x;
    let add: f32 = pad_amount * sky_map_pixel_size.x;
    coord.x = coord.x * mul + add;
    
    return coord;
}

fn unproject_sky(coord: vec2<f32>) -> vec3<f32> {
    // 从Uniform获取分辨率
    let sky_map_pixel_size: vec2<f32> = 1.0 / uniforms.iResolution;
    let pad_amount: f32 = 2.0;
    
    // 反向填充计算
    let inverse_mul: f32 = 1.0 / (1.0 - 2.0 * pad_amount * sky_map_pixel_size.x);
    let inverse_add: f32 = -pad_amount * sky_map_pixel_size.x * inverse_mul;
    let unpacked_x: f32 = coord.x * inverse_mul + inverse_add;
    
    // 保证坐标在[0,1)范围内
    var new_coord: vec2<f32> = vec2(fract(unpacked_x), coord.y);
    
    // 高度角非线性映射
    new_coord.y = select(
        -pow(1.0 - 2.0 * new_coord.y, 2.0),
        pow(2.0 * new_coord.y - 1.0, 2.0),
        new_coord.y >= 0.5
    );
    
    // 角度重建
    let azimuth_angle = new_coord.x * TAU - PI;
    let altitude_angle = new_coord.y * HALF_PI;
    
    // 三维坐标重建
    let altitude_cos = cos(altitude_angle);
    let altitude_sin = sin(altitude_angle);
    let azimuth_cos = cos(azimuth_angle);
    let azimuth_sin = sin(azimuth_angle);
    
    return vec3<f32>(
        altitude_cos * azimuth_sin,
        altitude_sin,
        -altitude_cos * azimuth_cos
    );
}

fn get_ambient_color(sun_color: vec3<f32>, moon_color: vec3<f32>) -> vec3<f32> {
    var sky_dir = normalize(vec3(0.0, 1.0, -0.8));
    var sky_color = atmosphere_scattering(sky_dir, sun_color, uniforms.sun_dir, moon_color, moon_dir);
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
	let nu = dot(ray_dir, uniforms.sun_dir);
	let alpha = vec3(0.429, 0.522, 0.614);
	let center_to_edge = max(nu - 0.998,0.0);
	let limb_darkening = pow(vec3(1.0 - sqr(1.0 - center_to_edge)), 0.5 * alpha);

	return 40.0 * sun_color * step(1e-5, center_to_edge);
}
#endif