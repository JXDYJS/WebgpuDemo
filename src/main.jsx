import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import './index.css'
//import App from './App.jsx'
import GpuDemo from './gpuDemo.jsx'

createRoot(document.getElementById('root')).render(
  <StrictMode>
    <GpuDemo />
  </StrictMode>,
)
