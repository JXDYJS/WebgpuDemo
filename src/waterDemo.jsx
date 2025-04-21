import { useEffect, useRef, useState } from 'react';

export default function WaterDemo() {
  const canvasRef = useRef(null);
  const [error, setError] = useState(null);
  const [shaderCode, setShaderCode] = useState(null);
  const animationFrameId = useRef(null);
  const deviceRef = useRef(null); // 新增设备引用
  const uniformBufferRef = useRef(null); // 新增uniform缓冲区引用
  const bindGroupRef = useRef(null); // 需要添加这行

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
    const initWebGPU = async () => {
      try {
        setError(null);
        
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
        
        const device = await adapter.requestDevice();
        deviceRef.current = device; // 保存设备引用
        device.lost.then(info => {
          throw new Error(`设备丢失: ${info.message}`);
        });
        //配置画布
        const dpr = window.devicePixelRatio || 1;
        canvas.width = 800 * dpr;
        canvas.height = 600 * dpr;
        canvas.style.width = '800px';
        canvas.style.height = '600px';

        // 4. 配置上下文
        const context = canvas.getContext('webgpu');
        if (!context) throw new Error("无法获取WebGPU上下文");
        
        const format = navigator.gpu.getPreferredCanvasFormat();
        context.configure({ device, format, alphaMode: 'opaque' });

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

    //  创建Uniform缓冲区
        const uniformBuffer = device.createBuffer({
            size: 256, // 对齐到256字节
            usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
        });
        uniformBufferRef.current = uniformBuffer;

        // 6. 创建绑定组布局
        const bindGroupLayout = device.createBindGroupLayout({
            entries: [{
                binding: 0,
                visibility: GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT,
                buffer: { type: 'uniform' }
            }]
            });
    
            // 7. 创建绑定组
            bindGroupRef.current = device.createBindGroup({
            layout: bindGroupLayout,
            entries: [{
                binding: 0,
                resource: { buffer: uniformBuffer }
            }]
            });

        // 6. 创建渲染管线
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

        // 7. 启动渲染循环
        const render = () => {
          try {
            const encoder = device.createCommandEncoder();
            const texture = context.getCurrentTexture();
            const uniformData = new Float32Array(16);
            const dataView = new DataView(uniformData.buffer);
            const time = (performance.now() - startTime) / 1000;
            dataView.setFloat32(0, time, true); // iTime offset 0
            dataView.setFloat32(8, canvas.width, true); // iResolution.x offset 8
            dataView.setFloat32(12, canvas.height, true); // iResolution.y offset 12
            device.queue.writeBuffer(
                uniformBuffer,
                0,
                uniformData.buffer,
                uniformData.byteOffset,
                uniformData.byteLength
              );
            
            if (!texture) throw new Error("无法获取当前纹理");
            
            const pass = encoder.beginRenderPass({
              colorAttachments: [{
                view: texture.createView(),
                clearValue: [0,0,0,1],
                loadOp: 'clear',
                storeOp: 'store'
              }]
            });

            pass.setPipeline(pipeline);
            pass.setVertexBuffer(0, vertexBuffer);
            
            if (vertexBuffer.size < 4 * 6 * 4) {
              throw new Error("顶点数据不完整");
            }
            
            pass.draw(4);
            pass.end();

            device.queue.submit([encoder.finish()]);
            animationFrameId.current = requestAnimationFrame(render);
          } catch (e) {
            setError(`渲染错误: ${e.message}`);
            cancelAnimationFrame(animationFrameId.current);
          }
        };

        render(); // 开始渲染循环

      } catch (error) { // 外层catch捕获所有初始化错误
        setError(error.message);
        console.error("初始化失败:", error);
      }
    };

    initWebGPU();

    // 清理函数
    return () => {
      if (animationFrameId.current) {
        cancelAnimationFrame(animationFrameId.current);
      }
      if (deviceRef.current) {
        deviceRef.current.destroy(); // 销毁设备
        deviceRef.current = null;
      }
    };
  }, []);

  return (
    <div className="gpu-demo">
      <h2>WebGPU Demo</h2>
      
      {error && (
        <div style={errorStyle}>
          <strong>⚠️ 发生错误:</strong>
          <div>{error}</div>
        </div>
      )}

      <canvas 
        ref={canvasRef}
        width={800}
        height={600}
        style={{ 
          border: `2px solid ${error ? '#ff4444' : '#333'}`,
          display: error ? 'none' : 'block',
          width: "1600px",   // 显示尺寸：1600×1200（放大 200%）
          height: "1200px",
          imageRendering: "crisp-edges"
        }}
      />
    </div>
  );
}