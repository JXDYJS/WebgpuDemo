
import { useEffect, useRef, useState } from 'react';

function angleToSunDirection(angle) {
  let c = Math.cos(angle / 180 * Math.PI);
  let s = Math.sin(angle / 180 * Math.PI);
  var dir = new Float32Array([c, s, 0.1]);
  let l = Math.sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
  dir = [dir[0] / l, dir[1] / l, dir[2] / l];
  return dir;
}
async function loadCombinedShaders({file1,file2}) {
  try {
    const [waterdemo, main] = await Promise.all([
      fetch(`./shader/${file1}`).then(r => r.text()),
      fetch(`./shader/${file2}`).then(r => r.text())
    ]);
    // 确保main内容在最后且只出现一次
    const combinedCode =`${waterdemo}\n//# MAIN_SHADER_START\n${main}`;
    return combinedCode
  } catch (error) {
    throw new Error(`Shader加载失败: ${error.message}`);
  }
}

export default function WaterDemo({ resolutionScale, SunAngle,roughness,metallic}) {
  const canvasRef = useRef(null);
  const [error, setError] = useState(null);
  const deviceRef = useRef(null);
  const contextRef = useRef(null);
  const uniformBufferRef = useRef(null);
  const bindGroupRef = useRef(null);
  const animationFrameId = useRef(null);
  const vertexBufferRef = useRef(null);
  const pipelineRef = useRef(null);
  const envTextureRef = useRef(null);
  const stagingTextureRef = useRef(null);
  const envSamplerRef = useRef(null);
  const prePipelineRef = useRef(null);
  const sunAngleRef = useRef(SunAngle);
  const roughnessRef = useRef(roughness);
  const metallicRef = useRef(metallic);
  const initializationLock = useRef(false);
  const startTimeRef = useRef(0);
  const env = true;

  useEffect(() => {
    sunAngleRef.current = SunAngle;
    roughnessRef.current = roughness;
    metallicRef.current = metallic;
  }, [SunAngle, roughness, metallic]);

  useEffect(() => {
    if (initializationLock.current) return;
    initializationLock.current = true;

    const initWebGPU = async () => {
      try {
        if (deviceRef.current) {
          await deviceRef.current.destroy();
          deviceRef.current = null;
        }

        const adapter = await navigator.gpu.requestAdapter();
        const device = await adapter.requestDevice();
        deviceRef.current = device;

        const canvas = canvasRef.current;
        const dpr = window.devicePixelRatio || 1;
        canvas.width = Math.floor(window.innerWidth * dpr * resolutionScale);
        canvas.height = Math.floor(window.innerHeight * dpr * resolutionScale);

        const context = canvas.getContext('webgpu');
        const format = navigator.gpu.getPreferredCanvasFormat();
        context.configure({ device, format, alphaMode: 'opaque' });
        contextRef.current = context;

        const vertices = new Float32Array([
          0, 0, 1,0,0,1, 1,0, 0,1,0,1,
          0,1, 0,0,1,1, 1,1, 1,1,0,1
        ]);
        vertexBufferRef.current = device.createBuffer({
          size: vertices.byteLength,
          usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST,
          mappedAtCreation: true
        });
        new Float32Array(vertexBufferRef.current.getMappedRange()).set(vertices);
        vertexBufferRef.current.unmap();

        uniformBufferRef.current = device.createBuffer({
          size: 256,
          usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
        });
        var shaderCode;
        const pre_shaderCode = await loadCombinedShaders({ file1: "waterdemo.wgsl", file2: "prepra.wgsl" });
        if(!env){
          shaderCode = await loadCombinedShaders({ file1: "waterdemo.wgsl", file2: "main.wgsl" });
        }
        else{
          shaderCode = await loadCombinedShaders({ file1: "waterdemo.wgsl", file2: "main_env.wgsl" });
        }

        
        const bindGroupLayout = device.createBindGroupLayout({
            entries: env? [
                { binding: 0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
                { binding: 1, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float', viewDimension: '3d' } },
                { binding: 2, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } },
                { binding: 3, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float' } },
                { binding: 4, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } }
            ]:[
                { binding:0, visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } },
                { binding:1, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float', viewDimension: '3d' } },
                { binding:2, visibility: GPUShaderStage.FRAGMENT, sampler: { type: 'filtering' } }
            ]
        });

        const sampler = device.createSampler({
          addressModeU: 'clamp-to-edge',
          magFilter: 'linear',
          minFilter: 'linear'
        });

        const scatteringTexture = await (async () => {
          const response = await fetch('./image/scattering.dat');
          const arrayBuffer = await response.arrayBuffer();
          const texture = device.createTexture({
            size: [32, 64, 32],
            format: 'rgba16float',
            usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST,
            dimension: '3d'
          });
          const srcData = new Uint16Array(arrayBuffer);
          const dstData = new Uint16Array(32 * 64 * 32 * 4);
          for (let i = 0, j = 0; i < srcData.length; i += 3, j += 4) {
            dstData[j] = srcData[i];
            dstData[j + 1] = srcData[i + 1];
            dstData[j + 2] = srcData[i + 2];
            dstData[j + 3] = 0x3C00;
          }
          device.queue.writeTexture(
            { texture },
            dstData.buffer,
            { offset: 0, bytesPerRow: 32 * 4 * 2, rowsPerImage: 64 },
            [32, 64, 32]
          );
          return texture;
        })();

        if (env) {
          // 创建半分辨率纹理
          const envTexture = device.createTexture({
              size: [
                  Math.floor(canvas.width * 0.5),
                  Math.floor(canvas.height * 0.5),
                  1
              ],
              format: 'rgba8unorm', 
              usage: GPUTextureUsage.RENDER_ATTACHMENT | 
                    GPUTextureUsage.TEXTURE_BINDING |
                    GPUTextureUsage.COPY_SRC
          });
          envTextureRef.current = envTexture;

          const stagingTexture = device.createTexture({
              size: [
                  Math.floor(canvas.width * 0.5),
                  Math.floor(canvas.height * 0.5),
                  1
              ],
              format: 'rgba8unorm', 
              usage: GPUTextureUsage.COPY_DST | 
                    GPUTextureUsage.TEXTURE_BINDING
          });
          stagingTextureRef.current = stagingTexture;

          envSamplerRef.current = device.createSampler({
            magFilter: 'linear',
            minFilter: 'linear'
          });

          if (!stagingTextureRef.current || !envSamplerRef.current) {
            throw new Error('环境纹理或采样器未初始化');
          }

        }
        const entries = [
            { binding: 0, resource: { buffer: uniformBufferRef.current } },
            { binding: 1, resource: scatteringTexture.createView() },
            { binding: 2, resource: sampler }
        ];
        if (env) {
            entries.push(
                { binding: 3, resource: stagingTextureRef.current.createView() },
                { binding: 4, resource: envSamplerRef.current }
            );
        }
        bindGroupRef.current = device.createBindGroup({ layout: bindGroupLayout, entries });

        if (env) {
          prePipelineRef.current = await device.createRenderPipelineAsync({
              layout: device.createPipelineLayout({ bindGroupLayouts: [bindGroupLayout] }),
              vertex: {
                  module: device.createShaderModule({ code: pre_shaderCode }),
                  entryPoint: 'vertex_main', // 预处理着色器的入口点
                  buffers: [{
                  arrayStride: 6 * 4,
                  attributes: [
                    { shaderLocation: 0, offset: 0, format: 'float32x2' },
                    { shaderLocation: 1, offset: 2 * 4, format: 'float32x4' }
                  ]
                }]
              },
              fragment: {
                  module: device.createShaderModule({ code: pre_shaderCode }),
                  entryPoint: 'fragment_main', // 预处理着色器的入口点
                  targets: [{ format: 'rgba8unorm' }] // 匹配预处理纹理格式
              },
              primitive: { topology: 'triangle-strip' }
          });
      }

        pipelineRef.current = await device.createRenderPipelineAsync({
          layout: device.createPipelineLayout({ bindGroupLayouts: [bindGroupLayout] }),
          vertex: {
            module: device.createShaderModule({ code: shaderCode }),
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
            module: device.createShaderModule({ code: shaderCode }),
            entryPoint: 'fragment_main',
            targets: [{ format }]
          },
          primitive: { topology: 'triangle-strip' }
        });

        startTimeRef.current = performance.now();
        startRenderLoop(device);

      } catch (error) {
        setError(error.message);
      }
    };

    initWebGPU();
      return () => {
    initializationLock.current = false;
    
    if (contextRef.current) {
      contextRef.current.unconfigure();
      contextRef.current = null;
    }
    if (deviceRef.current) {
      deviceRef.current.destroy();
      deviceRef.current = null;
    }
    if (animationFrameId.current) {
      cancelAnimationFrame(animationFrameId.current);
    }
    if (envTextureRef.current) {
    envTextureRef.current.destroy();
    envTextureRef.current = null;
    }
    if (envSamplerRef.current) {
        // 注意：Sampler 没有 destroy 方法，只需置null
        envSamplerRef.current = null;
    }
  };
}, [resolutionScale]);

  const startRenderLoop = (device) => {
    const render = () => {
      if (!deviceRef.current) return;

      try {
        const now = performance.now();
        const time = (now - startTimeRef.current) / 1000;
        const Roughness = roughnessRef.current;
        const Metallic = metallicRef.current;
        const sun_dir = angleToSunDirection(sunAngleRef.current);
        console.log(Roughness);

        const uniformData = new Float32Array(64);
        const dataView = new DataView(uniformData.buffer);
        dataView.setFloat32(0, time, true);
        dataView.setFloat32(8, canvasRef.current.width, true);
        dataView.setFloat32(12, canvasRef.current.height, true);
        dataView.setFloat32(16, sun_dir[0], true);
        dataView.setFloat32(20, sun_dir[1], true);
        dataView.setFloat32(24, sun_dir[2], true);
        dataView.setFloat32(32, Roughness, true);
        dataView.setFloat32(36, Metallic, true);


        device.queue.writeBuffer(uniformBufferRef.current, 0, uniformData);

        const encoder = device.createCommandEncoder();
        if (env) {
            // 执行预处理Pass
            const prePassEncoder = encoder.beginRenderPass({
                colorAttachments: [{
                    view: envTextureRef.current.createView(),
                    loadOp: 'clear',
                    clearValue: [0,0,0,1],
                    storeOp: 'store'
                }]
            });
            prePassEncoder.setPipeline(prePipelineRef.current);
            prePassEncoder.setVertexBuffer(0, vertexBufferRef.current);
            prePassEncoder.setBindGroup(0, bindGroupRef.current);
            prePassEncoder.draw(4);
            prePassEncoder.end();

            encoder.copyTextureToTexture(
            { texture: envTextureRef.current },
            { texture: stagingTextureRef.current }, // 需要创建中间过渡纹理
            [envTextureRef.current.width, envTextureRef.current.height]
    );
        }
        const pass = encoder.beginRenderPass({
          colorAttachments: [{
            view: contextRef.current.getCurrentTexture().createView(),
            loadOp: 'clear',
            storeOp: 'store'
          }]
        });

        pass.setPipeline(pipelineRef.current);
        pass.setVertexBuffer(0, vertexBufferRef.current);
        pass.setBindGroup(0, bindGroupRef.current);
        pass.draw(4);
        pass.end();

        device.queue.submit([encoder.finish()]);
        animationFrameId.current = requestAnimationFrame(render);
      } catch (e) {
        setError(`渲染错误: ${e.message}`);
      }
    };

    animationFrameId.current = requestAnimationFrame(render);
  };

  return (
    <div className="gpu-demo">
      {error && (
        <div style={{
          color: '#ff4444',
          backgroundColor: '#1a1a1a',
          padding: '1rem',
          borderRadius: '4px',
          margin: '1rem 0',
          whiteSpace: 'pre-wrap'
        }}>
          <strong>⚠️ 发生错误:</strong>
          <div>{error}</div>
        </div>
      )}
      <canvas ref={canvasRef} style={{
        position: 'fixed',
        top: 0,
        left: 0,
        width: '100%',
        height: '100%',
        imageRendering: 'crisp-edges'
      }}/>
    </div>
  );
}