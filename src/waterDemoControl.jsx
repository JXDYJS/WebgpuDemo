import { useState } from 'react';
import './WaterDemoControl.css';

export default function WaterDemoControl({ onResolutionChange, onSunAngleChange }) {
  const [resolution, setResolution] = useState(80);
  const [sunAngle, setSunAngle] = useState(45); // 默认正午角度

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

  return (
    <div className="control-panel">
      <div className="slider-container">
        <label>分辨率比例: {resolution}%</label>
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
        <label>太阳角度: {sunAngle}°</label>
        <input
          type="range"
          min="0"
          max="80"
          step="1"
          value={sunAngle}
          onChange={handleSunAngleChange}
          className="sun-angle-slider"
        />
      </div>
    </div>
  );
}