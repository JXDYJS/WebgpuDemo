import { useState } from 'react';
import './WaterDemoControl.css';

export default function WaterDemoControl({ onResolutionChange }) {
  const [resolution, setResolution] = useState(80);

  const handleChange = (e) => {
    const value = parseInt(e.target.value);
    setResolution(value);
    onResolutionChange(value / 100);
  };

  return (
    <div className="control-panel">
      <div className="slider-container">
        <label>分辨率比例: {resolution}%</label>
        <input
          type="range"
          min="25"
          max="100"
          step="5"
          value={resolution}
          onChange={handleChange}
          className="resolution-slider"
        />
      </div>
    </div>
  );
}