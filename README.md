图形渲染代码功能文档



一、概述



本代码主要实现了一个结合大气散射效果与物理渲染（PBR）的图形渲染功能，可用于模拟天空、海洋等自然场景。代码中包含了多种功能模块，如大气参数定义、PBR 光照计算、海面高度场生成、光线追踪等，通过一系列自定义函数实现复杂的渲染效果。


二、代码结构与模块功能



### 2.1 全局常量与结构体定义&#xA;



1.  **Uniforms 结构体**：用于存储全局统一变量，包括时间（`iTime`）、分辨率（`iResolution`）、太阳方向（`sun_dir`）以及球体相关数据（`ball_data`）。这些变量在渲染过程中可动态更新，影响场景表现。


2.  **Pos\_data 结构体**：存储位置相关数据，如球体位置（`ball_pos`）、相机位置（`camera_pos`）和光源位置（`light_pos`），但在现有代码中未被充分使用。


3.  **Material 结构体**：定义材质属性，包括反照率（`albedo`）、菲涅尔反射系数（`f0`）、粗糙度（`roughness`）和金属度（`metallic`），用于 PBR 光照计算。


4.  **TracingResult 结构体**：用于存储光线追踪结果，包含距离（`distance`）、位置（`position`）和命中类型（`hit_type`，0 代表海面，1 代表天空，2 代表球体）。


### 2.2 大气散射相关函数&#xA;



1.  **大气参数定义**：


*   定义了太阳、月亮相关参数，如方向（`sun_dir`、`moon_dir`）、强度（`SUN_I`、`MOON_I`）、颜色（`SUNLIGHT_COLOR`、`MOONLIGHT_COLOR`）等。


*   定义了行星和大气相关参数，如行星半径（`PLANET_RADIUS`）、大气内半径（`ATMOSPHERE_INNER_RADIUS`）、外半径（`ATMOSPHERE_OUTER_RADIUS`）等。


*   定义了查找表（LUT）相关参数，如散射分辨率（`SCATTERING_RES`）、最小太阳天顶角余弦值（`MIN_MU_S`）等。


1.  **大气散射计算函数**：


*   `FADE`、`SUN_WEIGHT`、`SUN_RAISE`、`SUN_SET`等函数用于计算太阳相关的权重和状态，如日出日落的过渡效果。


*   `get_sun_exposure`和`get_sun_tint`函数分别计算太阳的曝光强度和色调，考虑了日出日落以及蓝调时刻的影响。


*   `get_moon_exposure`和`get_moon_tint`函数计算月亮的曝光强度和色调。


*   `atmosphere_scattering`函数是核心大气散射计算函数，通过对太阳和月亮方向与光线方向的计算，结合查找表采样，计算大气散射的颜色贡献。


### 2.3 PBR 光照计算相关函数&#xA;



1.  **PBR 基础函数**：


*   实现了 GGX 法线分布函数（`distribution_ggx`）、Smith 单方向和双方向遮挡函数（`visibility_smith_ggx_single`、`visibility_smith_ggx_joint`）、菲涅尔近似函数（`fresnel_schlick`）以及介质菲涅尔计算函数（`fresnel_dielectric`）等基础 PBR 计算函数。


*   `f0_to_ior`函数将菲涅尔反射系数转换为折射率。


1.  **漫反射计算函数**：


*   `diffuse_hammon`函数基于 Hammon 模型计算漫反射，考虑了粗糙度、菲涅尔反射系数等因素。


*   `diffuse_water`函数专门针对水体的漫反射计算，考虑了双向菲涅尔透射和水体的吸收效应。


1.  **PBR 着色函数**：


*   `pbr_shading`和`pbr_shading_water`函数分别用于普通材质和水体材质的 PBR 着色计算，结合漫反射和高光反射计算最终的着色结果。


### 2.4 海面相关函数&#xA;



1.  **海面高度场生成**：


*   `sea_octave`函数通过噪声函数生成海面的一个八度高度场。


*   `map`和`map_detail`函数通过多次调用`sea_octave`函数，生成不同细节程度的海面高度场。


1.  **海面相关计算**：


*   `getNormal`函数根据海面高度场计算海面法线。


*   `getSeaColor`和`getSeaColor_Ph`函数计算海面颜色，考虑了菲涅尔反射、反射和折射成分、距离衰减以及高光等因素。


### 2.5 光线追踪相关函数&#xA;



1.  **光线追踪函数**：


*   `heightMapTracing`函数用于基于高度场的光线追踪，判断光线是否与海面相交，并返回相交距离和位置。


*   `is_hit_ball`函数用于判断光线是否与球体相交，并返回相交结果。


1.  **颜色计算函数**：


*   `get_color`和`get_color_without_sun`函数根据光线追踪结果计算最终颜色，考虑了天空颜色、海面颜色、大气散射、雾效等因素。


### 2.6 其他辅助函数&#xA;



1.  **坐标转换函数**：


*   `project_sky`和`unproject_sky`函数实现了三维方向向量与二维天空纹理坐标之间的转换。


1.  **噪声与哈希函数**：


*   `hash`函数生成哈希值，`noise`函数生成二维噪声。


1.  **环境纹理采样函数**：


*   `sample_env_texture`函数根据方向向量对环境纹理进行采样，`fast_gaussian_5x5`函数实现快速的 5x5 高斯滤波采样。


1.  **其他辅助计算函数**：


*   包含各种数学辅助函数，如`linear_step`、`pulse`、`from_srgb`、`smoothstep`等，用于数值计算和颜色空间转换。


三、使用说明





1.  **全局变量设置**：通过`Uniforms`结构体设置时间、分辨率、太阳方向等全局变量，这些变量将影响整个场景的渲染效果。


2.  **材质设置**：通过`Material`结构体设置物体材质属性，可根据需求调整反照率、粗糙度、金属度等参数。


3.  **光线追踪与颜色计算**：在主渲染逻辑中，调用`heightMapTracing`和`is_hit_ball`函数进行光线追踪，然后使用`get_color`或`get_color_without_sun`函数计算最终颜色。


4.  **环境纹理采样**：如需使用环境纹理，可通过`sample_env_texture`函数结合`pbr_env`函数进行环境光的采样和计算。


四、注意事项





1.  代码中的一些常量和参数（如大气参数、海面参数等）是根据特定场景需求设定的，如需调整渲染效果，可能需要修改这些参数。


2.  光线追踪和 PBR 计算涉及大量数学运算，可能对性能有一定要求，在实际应用中可根据硬件性能进行优化，如减少采样次数等。


3.  部分函数（如`get_ball_material`）预留了扩展接口，未来可根据需求进一步完善物体材质的动态设置功能。

