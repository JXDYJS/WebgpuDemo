import { useState } from 'react'
import reactLogo from './assets/react.svg'
import viteLogo from '/vite.svg'
import './App.css'
import WaterDemoControl from './waterDemoControl'
import WaterDemo from './waterDemo'

// function App() {
//   const [count, setCount] = useState(0)

//   return (
//     <>
//       <div>
//         <a href="https://vite.dev" target="_blank">
//           <img src={viteLogo} className="logo" alt="Vite logo" />
//         </a>
//         <a href="https://react.dev" target="_blank">
//           <img src={reactLogo} className="logo react" alt="React logo" />
//         </a>
//       </div>
//       <h1>Vite + React</h1>
//       <div className="card">
//         <button onClick={() => setCount((count) => count + 1)}>
//           count is {count}
//         </button>
//         <p>
//           Edit <code>src/App.jsx</code> and save to test HMR
//         </p>
//       </div>
//       <p className="read-the-docs">
//         Click on the Vite and React logos to learn more
//       </p>
//     </>
//   )
// }

function App() {
  const [resolutionScale, setResolutionScale] = useState(0.8);
  const [sunAngle, setSunAngle] = useState(45);
  const [roughness, setRoughness] = useState(0.25);
  const [metallic, setMetallic] = useState(0);

  return (
    <div>
      <WaterDemo resolutionScale={resolutionScale} SunAngle={sunAngle} roughness={roughness} metallic={metallic}/>
      <WaterDemoControl 
        onResolutionChange={(scale) => setResolutionScale(scale)} onSunAngleChange={(angle) => setSunAngle(angle)}
        onRoughnessChange={(roughness) => setRoughness(roughness)} onMetallicChange={(metallic) => setMetallic(metallic)}
      />
    </div>
  );
}

export default App
