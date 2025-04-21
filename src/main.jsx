import { StrictMode } from 'react'
import { createRoot } from 'react-dom/client'
import './index.css'
//import App from './App.jsx'
import WaterDemo from './waterDemo.jsx'

createRoot(document.getElementById('root')).render(
  <StrictMode>
    < WaterDemo/>
  </StrictMode>,
)
