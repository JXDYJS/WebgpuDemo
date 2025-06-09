import { useState } from 'react';
import './WaterDemoControl.css';

export default function WaterDemoControl({ 
  onResolutionChange, 
  onSunAngleChange,
  onRoughnessChange,
  onMetallicChange
}) {
  const [resolution, setResolution] = useState(80);
  const [sunAngle, setSunAngle] = useState(45);
  const [roughness, setRoughness] = useState(25);
  const [metallic, setMetallic] = useState(0);

  const handleResolutionChange = (e) => {
    const value = parseInt(e.target.value);
    setResolution(value);
    onResolutionChange(value / 100);
  };

  const handleSunAngleChange = (e) => {
    const value = parseInt(e.target.value);
    setSunAngle(value);
    onSunAngleChange(value);
  };

  const handleRoughnessChange = (e) => {
    const value = parseInt(e.target.value);
    setRoughness(value);
    onRoughnessChange(value / 100);
  };

  const handleMetallicChange = (e) => {
    const value = parseInt(e.target.value);
    setMetallic(value);
    onMetallicChange(value / 100);
  };

  return (
    <div className="control-panel">
      {/* 精细度和表面粗糙度在同一行 */}
      <div className="control-row">
        <div className="slider-container">
          <label>分辨率: {resolution}%</label>
          <input
            type="range"
            min="25"
            max="125"
            step="5"
            value={resolution}
            onChange={handleResolutionChange}
            className="resolution-slider"
          />
        </div>
        
        <div className="slider-container">
          <label>粗糙度: {roughness}%</label>
          <input
            type="range"
            min="5"
            max="95"
            step="1"
            value={roughness}
            onChange={handleRoughnessChange}
            className="roughness-slider"
          />
        </div>
      </div>
      
      {/* 太阳角度和金属压痕在同一行 */}
      <div className="control-row">
        <div className="slider-container">
          <label>太阳角度: {sunAngle}°</label>
          <input
            type="range"
            min="10"
            max="170"
            step="1"
            value={sunAngle}
            onChange={handleSunAngleChange}
            className="sun-angle-slider"
          />
        </div>
        
        <div className="slider-container">
          <label>金属度: {metallic}%</label>
          <input
            type="range"
            min="0"
            max="100"
            step="1"
            value={metallic}
            onChange={handleMetallicChange}
            className="metallic-slider"
          />
        </div>
      </div>
    </div>
  );
}