# Quantum Tunneling Simulator  

A Python-based Quantum Tunneling Simulator that models quantum wave packet evolution and quantum phase estimation using **Qiskit** and **Matplotlib**.

## Features  
- Simulates **quantum tunneling** of a wave packet using numerical methods.  
- Implements **Quantum Phase Estimation (QPE)** using Qiskit.  
- **Matplotlib animation** for visualizing wave packet evolution over time.  
- Probability visualization for tunneling transmission and reflection.  

## Installation  

Make sure you have **Python 3.8+** installed, then install the required dependencies:  

```bash
pip install numpy matplotlib qiskit qiskit-aer scipy
```

## Usage  

Run the main script to start the simulation:  

```bash
python app.py
```

This will:  
1. Set up a quantum barrier.  
2. Simulate wave packet evolution using numerical methods.  
3. Run QPE-based quantum simulation via Qiskit.  
4. Visualize the results using **Matplotlib animations** and **histograms**.  

## Example Output  

The simulation produces:  
- **Animated probability distribution** of the wave packet.  
- **Histogram of quantum phase estimation results.**  
- **Graph showing probability distribution on both sides of the barrier.**  

## Dependencies  
- `numpy`  
- `matplotlib`  
- `qiskit`  
- `qiskit-aer`  
- `scipy`  

## Requirement
- `numpy>=1.23.5`
- `matplotlib>=3.5.1`
- `qiskit==1.4.2`
- `qiskit-aer>=0.16.1`
- `scipy>=1.10.1`
