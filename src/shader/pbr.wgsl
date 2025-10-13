#ifndef PBR_WGSL
#define PBR_WGSL

#include "global.wgsl"
// GGX法线分布函数
fn distribution_ggx(NoH: f32, alpha_sq: f32) -> f32 {
    let denom = (1.0 - NoH * NoH) + NoH * NoH * alpha_sq;
    return alpha_sq / (PI * denom * denom);
}

// Smith单方向遮挡函数
fn visibility_smith_ggx_single(cos_theta: f32, alpha_sq: f32) -> f32 {
    let term1 = (-cos_theta * alpha_sq + cos_theta) * cos_theta;
    return 1.0 / (cos_theta + sqrt(term1 + alpha_sq));
}

// Smith双方向遮挡函数
fn visibility_smith_ggx_joint(NoL: f32, NoV: f32, alpha_sq: f32) -> f32 {
    let lambda_v = NoV * sqrt( (-NoL * alpha_sq + NoL) * NoL + alpha_sq );
    let lambda_l = NoL * sqrt( (-NoV * alpha_sq + NoV) * NoV + alpha_sq );
    return 1.0 / (1.0 + lambda_v + lambda_l + 1e-5); // 添加防除零保护
}

// 菲涅尔近似（Schlick）
fn fresnel_schlick(cos_theta: f32, f0: vec3<f32>) -> vec3<f32> {
    return f0 + (vec3<f32>(1.0) - f0) * pow5(1.0 - cos_theta);
}

// 介质菲涅尔计算
fn fresnel_dielectric(cos_theta: f32, ior: f32) -> f32 {
    let g_sq = sqr(ior) + sqr(cos_theta) - 1.0;
    if (g_sq < 0.0){return 1.0;}
    let g = sqrt(g_sq);
    let a = (g - cos_theta) / max((g + cos_theta),EPSILON);
    let b = (cos_theta * (g + cos_theta) - 1.0) / max((cos_theta * (g - cos_theta) + 1.0),EPSILON);
    return 0.5 * a * a * (1.0 + b * b);
}

fn diffuse_hammon(
    albedo: vec3<f32>,
    roughness: f32,
    f0: f32,
    NoL: f32,
    NoV: f32,
    VoH: f32,
    LoV: f32
) -> vec3<f32> {
    // 参数预处理
    let alpha = roughness * roughness;
    let alpha_sq = alpha * alpha;
    let sqrt_f0 = sqrt(f0) * 0.99999;
    let ior = (1 + sqrt_f0) / (1 - sqrt_f0);
    
    // 能量守恒因子
    let energy_factor = 1.0 - (4.0 * sqrt(f0) + 5.0 * f0 * f0) / 9.0;
    
    // 单次散射项
    let facing = 0.5 * LoV + 0.5;
    let fresnel_l = 1.0 - fresnel_dielectric(max(NoL, 0.01), ior);
    let fresnel_v = 1.0 - fresnel_dielectric(max(NoV, 0.01), ior);
    let single_smooth = (fresnel_l * fresnel_v) / max(energy_factor,EPSILON);
    
    // 粗糙表面修正项（论文启发式公式）
    let single_rough = max(facing, 0.0) * (-0.2 * facing + 0.45) * (1.0 / max(VoH,0.01) + 2.0);
    
    // 混合单次散射项
    let single = mix(single_smooth, single_rough, roughness) / PI;
    
    // 多次散射项（论文数值拟合）
    let multi = 0.1159 * roughness;
    
    // 最终组合
    return albedo * (multi +single);
}

fn diffuse_water(
    albedo: vec3<f32>,
    roughness: f32,
    f0: f32,
    L: vec3<f32>,
    V: vec3<f32>,
    H: vec3<f32>
) -> vec3<f32> {
    // 微表面法线相关
    let LoH = max(dot(L, H), 0.0);
    let VoH = max(dot(V, H), 0.0);
    let sqrt_f0 = sqrt(f0) * 0.99999;
    let ior = (1 + sqrt_f0) / (1 - sqrt_f0);
    
    // 双向菲涅尔透射
    let F_in = fresnel_dielectric(LoH,ior);
    let F_out = fresnel_dielectric(VoH,ior);
    let transmittance = (1.0 - F_in) * (1.0 - F_out);
    
    // 补偿因子（适应水的IOR）
    let compensation = 1.0 / (1.0 - 0.28 * roughness);
    
    // 单次散射项改造
    let single = transmittance * compensation / PI;
    
    // 多次散射项（水体的吸收效应）
    let absorption = exp(-albedo * (1.0 - roughness) * 2.0);
    let multi = 0.1159 * roughness * absorption;
    
    return (single + multi) * albedo;
}

fn pbr_shading_water(
    albedo: vec3<f32>,
    roughness: f32,
    f0:f32,
    metallic: f32,
    N: vec3<f32>,
    L: vec3<f32>,
    V: vec3<f32>,
    H: vec3<f32>
) -> pbr_shading_res {
    // 几何参数
    var res: pbr_shading_res;
    let NoL = max(dot(N, L), 0.0);
    let NoV = max(dot(N, V), 0.0);
    let NoH = max(dot(N, H), 0.0);
    let VoH = max(dot(V, H), 0.0);
    let LoV = max(dot(L, V), 0.0);
    
    // 高光项计算
    let alpha_sq = (roughness * roughness);
    let D = distribution_ggx(NoH, alpha_sq);
    let G = visibility_smith_ggx_joint(NoL, NoV, alpha_sq);
    let F = fresnel_schlick(VoH, mix(vec3<f32>(0.02), albedo, metallic));
    
    // 高光BRDF
    let specular_brdf = (D * G * F) / (4.0 * NoL * NoV + 1e-5);
    
    // 漫反射项
    let diffuse_brdf = diffuse_water(albedo, roughness, f0, L,V,H);
    
    let kD = (vec3<f32>(1.0) - F) * (1.0 - metallic);
    res.diffuse = diffuse_brdf * kD * NoL;
    res.specular = specular_brdf * NoL;
    return res;
}

fn pbr_shading(
    albedo: vec3<f32>,
    roughness: f32,
    f0:f32,
    metallic: f32,
    N: vec3<f32>,
    L: vec3<f32>,
    V: vec3<f32>,
    H: vec3<f32>
) -> pbr_shading_res {
    // 几何参数
    var res: pbr_shading_res;
    let NoL = max(dot(N, L), 0.0);
    let NoV = max(dot(N, V), 0.0);
    let NoH = max(dot(N, H), 0.0);
    let VoH = max(dot(V, H), 0.0);
    let LoV = max(dot(L, V), 0.0);
    
    // 高光项计算
    let alpha_sq = (roughness * roughness);
    let D = distribution_ggx(NoH, alpha_sq);
    let G = visibility_smith_ggx_joint(NoL, NoV, alpha_sq);
    let F = fresnel_schlick(VoH, mix(vec3<f32>(0.02), albedo, metallic));
    
    // 高光BRDF
    let specular_brdf = (D * G * F) / (4.0 * NoL * NoV + 1e-5);
    
    // 漫反射项
    let diffuse_brdf = diffuse_hammon(albedo, roughness, f0,NoL,NoV, VoH,LoV);
    
    let kD = (vec3<f32>(1.0) - F) * (1.0 - metallic);
    res.diffuse = diffuse_brdf * kD * NoL;
    res.specular = specular_brdf * NoL;
    return res;
}

fn lambertNoTangent(normal: vec3<f32>, uv: vec2<f32>) -> vec3<f32> {
    let theta = 6.283185 * uv.x;
    
    // 将v坐标从[0,1]映射到[-1,1]
    let adjusted_vy = 2.0 * uv.y - 1.0;
    
    // 生成球面坐标点
    let sphere_radius_xy = sqrt(1.0 - adjusted_vy * adjusted_vy);
    let sphere_point = vec3<f32>(
        sphere_radius_xy * cos(theta),
        sphere_radius_xy * sin(theta),
        adjusted_vy
    );
    
    // 合成并归一化方向向量
    return normalize(normal + sphere_point);
}

fn sample_ggx_vndf(viewer_dir: vec3<f32>, alpha: vec2<f32>, hash: vec2<f32>) -> vec3<f32> {
    // Improved GGX importance sampling function by zombye
    // https://ggx-research.github.io/publication/2023/06/09/publication-ggx.html

    // Transform viewer direction to the hemisphere configuration
    var v_dir = normalize(vec3<f32>(alpha * viewer_dir.xy, viewer_dir.z));

    // Sample a reflection direction off the hemisphere
    let tau: f32 = 6.2831853; // 2 * PI
    let phi: f32 = tau * hash.x;
    let cos_theta: f32 = fma(1.0 - hash.y, 1.0 + v_dir.z, -v_dir.z);
    let sin_theta: f32 = sqrt(clamp(1.0 - cos_theta * cos_theta, 0.0, 1.0));
    let reflected: vec3<f32> = vec3<f32>(
        vec2<f32>(cos(phi), sin(phi)) * sin_theta,
        cos_theta
    );

    // Evaluate halfway direction
    let halfway: vec3<f32> = reflected + v_dir;

    // Transform back to hemiellipsoid configuration
    return normalize(vec3<f32>(
        alpha * halfway.xy,
        halfway.z
    ));
}

fn get_tbn(n:vec3<f32>)->mat3x3<f32>{
    var tangent : vec3<f32>;
    if(n.y > 0.99){
        tangent = vec3<f32>(1.0,0.0,0.0);
    }
    else{
        tangent = normalize(cross(vec3<f32>(0.0,1.0,0.0),n));
    }
    let bitangent = normalize(cross(n,tangent));
    return mat3x3<f32>(tangent,bitangent,n);
}

fn pbr_env(    
    albedo: vec3<f32>,
    roughness: f32,
    f0:vec3<f32>,
    metallic: f32,
    N: vec3<f32>,
    L: vec3<f32>,
    V: vec3<f32>,
    H: vec3<f32>) -> vec3<f32>{

    let NoL = dot(N, L);
    let NoV = max(dot(N, V), 0.0);
    let NoH = max(dot(N, H), 0.0);
    let VoH = max(dot(V, H), 0.0);
    let LoV = max(dot(L, V), 0.0);
    var cout = vec3<f32>(0.0);

    if(NoL > 0.0){
        let alpha_sq = (roughness * roughness);
        //let D = distribution_ggx(NoH, alpha_sq);
        let v1 = visibility_smith_ggx_single(NoV, alpha_sq);
        let v2 = visibility_smith_ggx_joint(NoL, NoV, alpha_sq);
        var F: vec3<f32>;
        if(metallic > 0.001){
            F = fresnel_schlick(NoL,f0);
        }
        else{
            F = vec3(fresnel_dielectric(NoL,f0_to_ior(f0.x)));
        }
        cout += F * (2.0 * NoL * v2 / v1);
    }
    return cout;
}
#endif