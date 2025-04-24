const PI:f32 = 3.14159265359;
const AtmosphereHeight = 80000.0;
const PlaneRadius = 6371000.0;
const SCATTER_R:vec3<f32> = vec3<f32>(5.8e-6, 1.35e-5, 3.3e-5);
const SCATTER_M:vec3<f32> = vec3<f32>(3.996e-6, 3.996e-6, 3.996e-6);
const EXTINCTION_M:vec3<f32> = vec3<f32>(4.4e-6, 4.4e-6, 4.4e-6);
const DensityScaleHeight:vec2<f32> = vec2<f32>(7994.0, 8192.0);

struct Uniforms {
    @group(0) @binding(0) var DensityTexture: texture_2d<f32>;
    @group(0) @binding(1) var LinearSampler: sampler;
}

//求交函数，返回.x为第一个交点，.y为第二个交点
fn RaysphereIntersect(ori:vec3<f32>, dir:vec3<f32>, center:vec3<f32>, radius:f32) -> vec2<f32> {
    let rayOrigin = ori - center;
    let a = dot(dir, dir);
    let b = dot(rayOrigin, dir);
    let c = dot(rayOrigin, rayOrigin) - radius * radius;
    var d = b * b - 4 * a * c;

    if(d < 0.0){
        //无交点
        return vec2<f32>(-1.0, -1.0);
    }
    else{
        d = sqrt(d);
        return vec2<f32>((-b - d) / (2.0 * a), (-b + d) / (2.0 * a));
    }
}

fn ComputeDensity(ori:vec3<f32>,dir:vec3<f32>)->vec2<f32>{
    let planeCenter = vec3<f32>(0.0,-PlaneRadius,0.0);

    let stepCount = 250.0;
    var intersect = RaysphereIntersect(ori,dir,planeCenter,PlaneRadius);
    if(intersect.x > 0.0){
        //打到了地面
        return vec2(1e+20);
    }
    intersect = RaysphereIntersect(ori,dir,planeCenter,AtmosphereHeight + PlaneRadius);
    let end = ori + intersect.y * dir;
    let step = (end - ori) / stepCount;
    let step_size = length(step);
    var density:vec2<f32> = vec2<f32>(0.0,0.0);
    for(var i = 0.5;i < stepCount;i++){
        var p = ori + step * i;
        var h = abs(length(p - planeCenter) - PlaneRadius);
        var hh = vec2(h,h);
        density += exp(-hh / DensityScaleHeight) * step_size;
    }
    return density;
}

struct GetDensityResult{
    localDensity:vec2<f32>,
    atomTopDensity:vec2<f32>
};

fn GetDensity(ori:vec3<f32>,planeCenter:vec3<f32>,light_dir:vec3<f32>)->GetDensityResult{
    var result:GetDensityResult;
    var h = length(ori - planeCenter) - PlaneRadius;
    result.localDensity = exp(-vec2(h,h) / DensityScaleHeight);

    let cosAngle = dot(normalize(ori - planeCenter), -light_dir);
    result.atomTopDensity = 
}


