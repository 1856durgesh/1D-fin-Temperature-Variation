
function relode(){
  setTimeout(() => {
    location.reload();
  }, 100);
}
// JavaScript logic for solving the heat conduction problem

function longFin() {

  // Get user inputs
 
  const rho = parseFloat(document.getElementById('rho').value);
  const Cp = parseFloat(document.getElementById('Cp').value);
  const k = parseFloat(document.getElementById('k').value);
  const h = parseFloat(document.getElementById('h').value);
  const r = parseFloat(document.getElementById('r').value);
  const T_inf = parseFloat(document.getElementById('T_inf').value);
  const L = parseFloat(document.getElementById('L').value);
  const dx = parseFloat(document.getElementById('dx').value);
  const dt = parseFloat(document.getElementById('dt').value);
  const Nt = parseInt(document.getElementById('Nt').value);
  const T2 = parseInt(document.getElementById('T2').value);

  const P = 2 * Math.PI * r; // Perimeter of the circular fin
  const A = Math.PI * r * r; // Area of the circular fin
  const del_p = (h * P * T_inf) / (rho * Cp * A);
  const N = Math.floor(L / dx); // Number of nodes
  const a_w = -(k * dt) / (rho * Cp * dx * dx);
  const a_e = a_w;
  const a_p = 1.0 - 2 * a_w + (h * P) / (rho * Cp * A);

  // Initialize temperature distribution
  let T = new Array(N).fill(T_inf);
  T[0] = T2; // Base temperature

  // TDMA algorithm
  function TDMA(a, b, c, d) {
    let n = b.length;
    let Tn1 = new Array(n);
    let cp = new Array(n - 1);
    let dp = new Array(n);

    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for (let i = 1; i < n; i++) {
      let m = 1.0 / (b[i] - a[i - 1] * cp[i - 1]);
      if (i < n - 1) cp[i] = c[i] * m;
      dp[i] = (d[i] - a[i - 1] * dp[i - 1]) * m;
    }

    Tn1[n - 1] = dp[n - 1];
    for (let i = n - 2; i >= 0; i--) {
      Tn1[i] = dp[i] - cp[i] * Tn1[i + 1];
    }
    return Tn1;
  }

  // Time stepping loop
  for (let n = 0; n < Nt; n++) {
    let a = new Array(N - 1).fill(a_w);
    let b = new Array(N).fill(a_p);
    let c = new Array(N - 1).fill(a_e);
    let d = T.map((t) => t + del_p);

    b[0] -= a_w;
    b[N - 1] = a_p - a_e;
    d[0] -= 2 * a_w * 100;
    let T2=T_inf;
    d[N - 1] = T[N-1] + del_p - 2 * a_e * T2;


    // Solve TDMA
    T = TDMA(a, b, c, d);
  }

  // Generate the data for the graph
  const x = [];
  const temperature = [];
  x.push(0);
  temperature.push(T2);
  for (let i = 0; i < N; i++) {
    x.push(i * dx +dx/2);
    temperature.push(T[i]);
  }

  // Plot the results using Chart.js
  plotGraph(x, temperature);

}


function insulatedFin() {
  // Get user inputs
  const rho = parseFloat(document.getElementById('rho').value);
  const Cp = parseFloat(document.getElementById('Cp').value);
  const k = parseFloat(document.getElementById('k').value);
  const h = parseFloat(document.getElementById('h').value);
  const r = parseFloat(document.getElementById('r').value);
  const T_inf = parseFloat(document.getElementById('T_inf').value);
  const L = parseFloat(document.getElementById('L').value);
  const dx = parseFloat(document.getElementById('dx').value);
  const dt = parseFloat(document.getElementById('dt').value);
  const Nt = parseInt(document.getElementById('Nt').value);
  const T2 = parseInt(document.getElementById('T2').value);


  const P = 2 * Math.PI * r; // Perimeter of the circular fin
  const A = Math.PI * r * r; // Area of the circular fin
  const del_p = (h * P * T_inf) / (rho * Cp * A);
  const N = Math.floor(L / dx); // Number of nodes
  const a_w = -(k * dt) / (rho * Cp * dx * dx);
  const a_e = a_w;
  const a_p = 1.0 - 2 * a_w + (h * P) / (rho * Cp * A);

  // Initialize temperature distribution
  let T = new Array(N).fill(T_inf);
  T[0] = T2; // Base temperature

  // TDMA algorithm
  function TDMA(a, b, c, d) {
    let n = b.length;
    let Tn1 = new Array(n);
    let cp = new Array(n - 1);
    let dp = new Array(n);

    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for (let i = 1; i < n; i++) {
      let m = 1.0 / (b[i] - a[i - 1] * cp[i - 1]);
      if (i < n - 1) cp[i] = c[i] * m;
      dp[i] = (d[i] - a[i - 1] * dp[i - 1]) * m;
    }

    Tn1[n - 1] = dp[n - 1];
    for (let i = n - 2; i >= 0; i--) {
      Tn1[i] = dp[i] - cp[i] * Tn1[i + 1];
    }
    return Tn1;
  }

  // Time stepping loop
  for (let n = 0; n < Nt; n++) {
    let a = new Array(N - 1).fill(a_w);
    let b = new Array(N).fill(a_p);
    let c = new Array(N - 1).fill(a_e);
    let d = T.map((t) => t + del_p);

    b[0] -= a_w;
    b[N - 1] = a_p - a_e;
    d[0] -= 2 * a_w * 100;
    let T2=T[N-1];
    d[N - 1] = T[N-1] + del_p - 2 * a_e * T2;


    // Solve TDMA
    T = TDMA(a, b, c, d);
  }

  // Generate the data for the graph
  const x = [];
  const temperature = [];
  x.push(0);
  temperature.push(T2);
  for (let i = 0; i < N; i++) {
    x.push(i * dx +dx/2);
    temperature.push(T[i]);
  }

  // Plot the results using Chart.js
  plotGraph(x, temperature);
}

function convectionFin() {
  // Get user inputs
  const rho = parseFloat(document.getElementById('rho').value);
  const Cp = parseFloat(document.getElementById('Cp').value);
  const k = parseFloat(document.getElementById('k').value);
  const h = parseFloat(document.getElementById('h').value);
  const r = parseFloat(document.getElementById('r').value);
  const T_inf = parseFloat(document.getElementById('T_inf').value);
  const L = parseFloat(document.getElementById('L').value);
  const dx = parseFloat(document.getElementById('dx').value);
  const dt = parseFloat(document.getElementById('dt').value);
  const Nt = parseInt(document.getElementById('Nt').value);
  const T2 = parseInt(document.getElementById('T2').value);


  const P = 2 * Math.PI * r; // Perimeter of the circular fin
  const A = Math.PI * r * r; // Area of the circular fin
  const del_p = (h * P * T_inf) / (rho * Cp * A);
  const N = Math.floor(L / dx); // Number of nodes
  const a_w = -(k * dt) / (rho * Cp * dx * dx);
  const a_e = a_w;
  const a_p = 1.0 - 2 * a_w + (h * P) / (rho * Cp * A);

  // Initialize temperature distribution
  let T = new Array(N).fill(T_inf);
  T[0] = T2; // Base temperature

  // TDMA algorithm
  function TDMA(a, b, c, d) {
    let n = b.length;
    let Tn1 = new Array(n);
    let cp = new Array(n - 1);
    let dp = new Array(n);

    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    for (let i = 1; i < n; i++) {
      let m = 1.0 / (b[i] - a[i - 1] * cp[i - 1]);
      if (i < n - 1) cp[i] = c[i] * m;
      dp[i] = (d[i] - a[i - 1] * dp[i - 1]) * m;
    }

    Tn1[n - 1] = dp[n - 1];
    for (let i = n - 2; i >= 0; i--) {
      Tn1[i] = dp[i] - cp[i] * Tn1[i + 1];
    }
    return Tn1;
  }

  // Time stepping loop
  for (let n = 0; n < Nt; n++) {
    let a = new Array(N - 1).fill(a_w);
    let b = new Array(N).fill(a_p);
    let c = new Array(N - 1).fill(a_e);
    let d = T.map((t) => t + del_p);

    b[0] -= a_w;
    b[N - 1] = a_p - a_e;
    d[0] -= 2 * a_w * 100;
    let  T2=(T_inf-(2*k*T[N-1])/(h*dx))/(1-2*k/(h*dx));
    d[N - 1] = T[N-1] + del_p - 2 * a_e * T2;

    // Solve TDMA
    T = TDMA(a, b, c, d);
  }

  // Generate the data for the graph
  const x = [];
  const temperature = [];
  x.push(0);
  temperature.push(T2);
  for (let i = 0; i < N; i++) {
    x.push(i * dx + dx/2);
    temperature.push(T[i]);
  }

  // Plot the results using Chart.js
  plotGraph(x, temperature);
}



// Function to plot the temperature distribution
function plotGraph(x, temperature) {
  const ctx = document.getElementById('temperatureChart').getContext('2d');
  new Chart(ctx, {
    type: 'line',
    data: {
      labels: x,
      datasets: [
        {
          label: 'Temperature vs. Position',
          data: temperature,
          borderColor: 'rgba(75, 192, 192, 1)',
          borderWidth: 2,
          fill: false,
        },
      ],
    },
    options: {
      scales: {
        x: {
          title: {
            display: true,
            text: 'Position along the fin (m)',
          },
        },
        y: {
          title: {
            display: true,
            text: 'Temperature (°C)',
          },
        },
      },
    },
  });
}
