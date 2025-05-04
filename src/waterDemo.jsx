import { useEffect, useRef, useState } from 'react';

function angleToSunDirection(angle) {

  const t = (angle - 20) / 160;

  const y = 1.0 / (5 * t + 1.2); 
  
  const x = Math.sin(t * Math.PI * 1.8) * 0.6;
  const z = Math.cos(t * Math.PI * 1.3) * 1.2;

  const length = Math.sqrt(x*x + y*y + z*z);
  return [x/length, y/length, z/length];
}

export default function WaterDemo({ resolutionScale,SunAngle}) { // 修改1: 参数名改为resolutionScale
  const canvasRef = useRef(null);
  const [error, setError] = useState(null);
  const [shaderCode, setShaderCode] = useState(null);
  const animationFrameId = useRef(null);
  const deviceRef = useRef(null);
  const contextRef = useRef(null);
  const uniformBufferRef = useRef(null);
  const bindGroupRef = useRef(null);
  const atmosphereScatteringLUT = useRef(null);
  const lastFrameRef = useRef(0);

  const deviceGeneration = useRef(0);
const isUnmounted = useRef(false);

useEffect(() => {
  isUnmounted.current = false;
  return () => { isUnmounted.current = true; };
}, []);

  const initializationLock = { current: false };


  const TARGET_FPS = 30;
  const FRAME_INTERVAL = 1000 / TARGET_FPS;

  const errorStyle = {
    color: '#ff4444',
    backgroundColor: '#1a1a1a',
    padding: '1rem',
    borderRadius: '4px',
    margin: '1rem 0',
    whiteSpace: 'pre-wrap'
  };

  const checkShaderCompilation = async (shaderModule) => {
    if (!shaderModule.compilationInfo) {
      console.warn("浏览器不支持 compilationInfo()，跳过着色器编译检查");
      return;
    }
    const compilationInfo = await shaderModule.compilationInfo();
    if (compilationInfo.messages.some(msg => msg.type === 'error')) {
      const errors = compilationInfo.messages
        .filter(msg => msg.type === 'error')
        .map(msg => `[${msg.lineNum}:${msg.linePos}] ${msg.message}`)
        .join('\n');
      throw new Error(`Shader编译错误:\n${errors}`);
    }
  };

  useEffect(() => {
      if (initializationLock.current) return;
      initializationLock.current = true;
    const initWebGPU = async () => {
      try {
        const generationBeforeInit = ++deviceGeneration.current;
        const init = async () => {
          if (animationFrameId.current) {
            cancelAnimationFrame(animationFrameId.current);
            animationFrameId.current = null;
          }
          
          // 异步等待设备销毁
          if (deviceRef.current) {
            const oldDevice = deviceRef.current;
            deviceRef.current = null;
            await oldDevice.destroy(); 
          }
    
          // 解除上下文绑定
          if (contextRef.current) {
            contextRef.current.unconfigure();
            contextRef.current = null;
          }
    
          await new Promise(requestAnimationFrame); // 确保GPU进程完成
        };
        
        // 1. 加载着色器
        const code = await import('./shader/waterdemo.wgsl?raw')
          .then(m => m.default)
          .catch(e => { throw new Error(`无法加载着色器: ${e.message}`) });
        setShaderCode(code);

        // 2. 检查WebGPU支持
        if (!navigator.gpu) throw new Error("请使用Chrome 113+或Edge 113+");
        
        // 3. 初始化设备
        const canvas = canvasRef.current;
        if (!canvas) throw new Error("找不到canvas元素");
        
        const adapter = await navigator.gpu.requestAdapter();
        if (!adapter) throw new Error("无法获取WebGPU适配器");

        if (contextRef.current) contextRef.current.unconfigure();
        
        const device = await adapter.requestDevice();
        deviceRef.current = device;
        device.lost.then(info => {
          throw new Error(`设备丢失: ${info.message}`);
        });

        // 修改2: 分辨率计算增加取整
        const dpr = window.devicePixelRatio || 1;
        const screenWidth = Math.floor(window.innerWidth * dpr * resolutionScale); 
        const screenHeight = Math.floor(window.innerHeight * dpr * resolutionScale);

        canvas.width = screenWidth;
        canvas.height = screenHeight;

        canvas.style.width = '100vw';
        canvas.style.height = '100vh';

        // 4. 配置上下文（修改3: 每次分辨率变化时重新配置）
        const context = canvas.getContext('webgpu');
        if (!context) throw new Error("无法获取WebGPU上下文");
        contextRef.current = context;
        
        const format = navigator.gpu.getPreferredCanvasFormat();
        context.configure({ 
          device, 
          format,
          size: [canvas.width, canvas.height], // 使用新尺寸
          alphaMode: 'opaque' 
        });


        const load3DTexture = async (device) => {
          try {
            // 1. 加载二进制文件
            const response = await fetch('./image/scattering.dat');
            if (!response.ok) throw new Error(`HTTP error! status: ${response.status}`);
            const arrayBuffer = await response.arrayBuffer();
            
            // 2. 验证数据尺寸 (32x64x32 RGB16F)
            const expectedSize = 32 * 64 * 32 * 6; // 3 channels × 2 bytes each
            if (arrayBuffer.byteLength !== expectedSize) {
              throw new Error(`Invalid texture data size. Expected ${expectedSize}, got ${arrayBuffer.byteLength}`);
            }
            // 3. 创建GPU纹理
            const texture = device.createTexture({
              size: [32, 64, 32],  // width, height, depth
              format: 'rgba16float', // WebGPU不支持rgb16float，需要转换为RGBA
              usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST,
              dimension: '3d'
            });
            atmosphereScatteringLUT.current = texture;
        
            // 4. 将RGB数据转换为RGBA（添加Alpha通道）
            const srcData = new Uint16Array(arrayBuffer);

            console.log('数据示例:', {
              R: srcData[0].toString(16), // 应为半精度浮点
              G: srcData[1].toString(16),
              B: srcData[2].toString(16)
            });

            const dstData = new Uint16Array(32 * 64 * 32 * 4); // RGBA
            
            for (let i = 0, j = 0; i < srcData.length; i += 3, j += 4) {
              dstData[j] = srcData[i];     // R
              dstData[j + 1] = srcData[i + 1]; // G
              dstData[j + 2] = srcData[i + 2]; // B
              dstData[j + 3] = 0x3C00;     // Alpha=1.0 in half float (0x3C00)
            }
        
            // 5. 上传数据到GPU
            device.queue.writeTexture(
              { texture },
              dstData.buffer,
              {
                offset: 0,
                bytesPerRow: 32 * 4 * 2,   // 32 pixels × 4 channels × 2 bytes
                rowsPerImage: 64
              },
              [32, 64, 32] // 完整3D尺寸
            );
        
            return texture;
          } catch (e) {
            throw new Error(`Failed to load 3D texture: ${e.message}`);
          }
        };


        // 5. 创建顶点缓冲区
        const vertices = new Float32Array([
          0, 0,   1,0,0,1,  1,0, 0,1,0,1,
          0,1, 0,0,1,1,  1,1, 1,1,0,1
        ]);

        const vertexBuffer = device.createBuffer({
          size: vertices.byteLength,
          usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST,
          mappedAtCreation: true
        });
        new Float32Array(vertexBuffer.getMappedRange()).set(vertices);
        vertexBuffer.unmap();

        // 创建Uniform缓冲区
        const uniformBuffer = device.createBuffer({
          size: 256,
          usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
        });
        uniformBufferRef.current = uniformBuffer;

        // 6. 创建绑定组布局
        const bindGroupLayout = device.createBindGroupLayout({
          entries: [
            {
              binding: 0,
              visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT,
              buffer: { type: 'uniform' }
            },
            {
              binding: 1,
              visibility: GPUShaderStage.FRAGMENT,
              texture: {
                sampleType: 'float',
                viewDimension: '3d'
              }
            },
            {
              binding: 2,
              visibility: GPUShaderStage.FRAGMENT,
              sampler: { type: 'filtering' }
            }
          ]
        });

        const sampler = device.createSampler({
          addressModeU: 'clamp-to-edge',
          addressModeV: 'clamp-to-edge',
          addressModeW: 'clamp-to-edge',
          magFilter: 'linear',
          minFilter: 'linear'
        });
    
        // 7. 创建绑定组
        const scatteringTexture = await load3DTexture(device);
        bindGroupRef.current = device.createBindGroup({
          layout: bindGroupLayout,
          entries: [
            {
              binding: 0,
              resource: { buffer: uniformBuffer }
            },
            {
              binding: 1,
              resource: scatteringTexture.createView()
            },
            {
              binding: 2,
              resource: sampler
            }
          ]
        });

        // 创建渲染管线
        const vertexShader = device.createShaderModule({ code });
        const fragmentShader = device.createShaderModule({ code });

        await checkShaderCompilation(vertexShader)
          .catch(e => { throw new Error(`顶点着色器错误: ${e.message}`) });
          
        await checkShaderCompilation(fragmentShader)
          .catch(e => { throw new Error(`片段着色器错误: ${e.message}`) });

        const pipeline = await device.createRenderPipelineAsync({
          layout: device.createPipelineLayout({
            bindGroupLayouts: [bindGroupLayout]
          }),
          vertex: {
            module: vertexShader,
            entryPoint: 'vertex_main',
            buffers: [{
              arrayStride: 6 * 4,
              attributes: [
                { shaderLocation: 0, offset: 0, format: 'float32x2' },
                { shaderLocation: 1, offset: 2 * 4, format: 'float32x4' }
              ]
            }]
          },
          fragment: {
            module: fragmentShader,
            entryPoint: 'fragment_main',
            targets: [{ format }]
          },
          primitive: {
            topology: 'triangle-strip'
          }
        }).catch(e => { throw new Error(`管线创建失败: ${e.message}`) });

        // 启动渲染循环
        const startTime = performance.now();
        let lastFrameTime = 0;
        
        const render = (timestamp) => {
          if (!deviceRef.current || deviceRef.current.destroyed) {
            return;
          }
          const currentDevice = deviceRef.current;
          try {
            const now = performance.now();
            const elapsed = now - lastFrameRef.current;

            const encoder = device.createCommandEncoder();
            const texture = context.getCurrentTexture();
            if (texture.device?.id !== deviceRef.current?.id) {
              throw new Error("纹理与当前设备不匹配！");
            }
            texture.autoReleaseAfterRender = true;
            const uniformData = new Float32Array(32);
            const dataView = new DataView(uniformData.buffer);

            // 修改4: 使用当前canvas的实际尺寸
            const time = (now - startTime) / 1000;
            const sun_dir = angleToSunDirection(SunAngle);
            dataView.setFloat32(0, time, true);
            dataView.setFloat32(8, canvas.width, true); 
            dataView.setFloat32(12, canvas.height, true);
            dataView.setFloat32(16, sun_dir[0], true);          // sun_dir.x (16)
            dataView.setFloat32(20, sun_dir[1], true);          // sun_dir.y (20)
            dataView.setFloat32(24, sun_dir[2], true);          // sun_dir.z (24)
            dataView.setFloat32(28, 0, true);                   // vec3填充 (28-31)
            
            device.queue.writeBuffer(
              uniformBuffer,
              0,
              uniformData.buffer,
              uniformData.byteOffset,
              uniformData.byteLength
            );

            lastFrameRef.current = now; 
            
            if (!texture) throw new Error("无法获取当前纹理");
            
            const pass = encoder.beginRenderPass({
              colorAttachments: [{
                view: texture.createView(),
                loadOp: 'clear',
                storeOp: 'store'
              }]
            });

            pass.setPipeline(pipeline);
            pass.setVertexBuffer(0, vertexBuffer);
            
            if (vertexBuffer.size < 4 * 6 * 4) {
              throw new Error("顶点数据不完整");
            }
            
            pass.setBindGroup(0, bindGroupRef.current);
            pass.draw(4);
            pass.end();

            const commandBuffer = encoder.finish(); 
            device.queue.submit([commandBuffer]);
            //encoder.destroy();
            animationFrameId.current = requestAnimationFrame(render);
          } catch (e) {
            setError(`渲染错误: ${e.message}`);
            cancelAnimationFrame(animationFrameId.current);
          }
        };

        // 初始清除
        const initialEncoder = device.createCommandEncoder();
        const initialPass = initialEncoder.beginRenderPass({
          colorAttachments: [{
            view: context.getCurrentTexture().createView(),
            clearValue: [0, 1, 0, 1],
            loadOp: 'clear',
            storeOp: 'store'
          }]
        });
        initialPass.end();
        device.queue.submit([initialEncoder.finish()]);

        render();

      } catch (error) {
        setError(error.message);
        console.error("初始化失败:", error);
      }
    };

    initWebGPU();

    return () => {
      if (animationFrameId.current) {
        cancelAnimationFrame(animationFrameId.current);
      }
      // 修改5: 仅当组件卸载时销毁设备
      if (deviceRef.current) {
        deviceRef.current.destroy();
        deviceRef.current = null;
      }
    };
  }, [resolutionScale]); // 修改6: 添加依赖项

  return (
    <div className="gpu-demo">
      {error && (
        <div style={errorStyle}>
          <strong>⚠️ 发生错误:</strong>
          <div>{error}</div>
        </div>
      )}

      <canvas 
        ref={canvasRef}
        style={{ 
          position: 'fixed',
          top: 0,
          left: 0,
          width: '100%',
          height: '100%',
          imageRendering: 'crisp-edges',
          border: 'none'
        }}
      />
    </div>
  );
}