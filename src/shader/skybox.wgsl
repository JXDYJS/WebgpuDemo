// 常量定义 -----------------------------------------------------------------
const PI: f32 = 3.141592653589793;
const TAU: f32 = 6.283185307179586;
const DEGREE: f32 = PI / 180.0;

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

// 纹理绑定
@group(0) @binding(0) var atmosphere_scattering_lut: texture_3d<f32>;
@group(0) @binding(1) var atmosphere_scattering_sampler: sampler;

// 工具函数 ----------------------------------------------------------------
fn smoothstep(low: f32, high: f32, x: f32) -> f32 {
    let t = clamp((x - low) / (high - low), 0.0, 1.0);
    return t * t * (3.0 - 2.0 * t);
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

    // 纹理采样 ------------------------------------------------------------
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

    // 后处理饱和度增强
    let luma_weights = vec3<f32>(0.2627, 0.6780, 0.0593); // Rec.2020
    let luma = dot(atmosphere, luma_weights);
    let saturation_factor = 1.5; 
    atmosphere = mix(vec3(luma), atmosphere, saturation_factor);

    return atmosphere;
}