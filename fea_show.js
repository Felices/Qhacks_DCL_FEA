for(i = 0; i < pointCount; i++) 
{
   r = 10 * Math.cos(i / 10);
   x.push(i);
   y.push(i);
   z.push(i);
   c.push(i)
}

Plotly.plot('graph', [{
  type: 'scatter3d',
  mode: 'lines',
  x: x,
  y: y,
  z: z,
  line: {
    width: 6,
    color: c,
    colorscale: "Viridis"},
  marker: {
    size: 3.5,
    color: c,
    colorscale: "Greens",
    cmin: -20,
    cmax: 50
  }},                  
], {}, {showSendToCloud: true});