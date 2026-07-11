import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from qiskit import QuantumCircuit, transpile
from qiskit_aer import Aer
from qiskit.visualization import plot_histogram
from scipy.linalg import expm
from qiskit.circuit.library import QFT, UnitaryGate


class QuantumTunnelingSimulator:
    def __init__(self, num_qubits=5, dt=0.1, total_time=10.0):
        self.num_qubits = num_qubits
        self.dim = 2**num_qubits
        self.dt = dt
        self.total_time = total_time
        self.steps = int(total_time / dt)
        self.x_range = np.linspace(0, 1, self.dim)
        self.probabilities = np.zeros((self.steps + 1, self.dim))

    
    def setup_barrier(self, barrier_start=0.4, barrier_end=0.6, barrier_height=1.0):
        self.potential = np.zeros(self.dim)
        barrier_start_idx = int(barrier_start * self.dim)
        barrier_end_idx = int(barrier_end * self.dim)
        self.potential[barrier_start_idx:barrier_end_idx] = barrier_height
        self.barrier_start = barrier_start
        self.barrier_end = barrier_end
        self.barrier_height = barrier_height

    
    def create_hamiltonian(self, wavepacket_energy=0.5):
        kinetic = np.zeros((self.dim, self.dim))
        for i in range(self.dim):
            kinetic[i, i] = -2 * wavepacket_energy
            if i > 0:
                kinetic[i, i-1] = wavepacket_energy
            if i < self.dim - 1:
                kinetic[i, i+1] = wavepacket_energy
        potential = np.diag(self.potential)
        hamiltonian = kinetic + potential
        self.evolution_operator = expm(-1j * hamiltonian * self.dt)

    
    def create_initial_wavepacket(self, center=0.2, width=0.05, k0=30):
        x = self.x_range
        psi = np.exp(-(x - center)**2 / (2 * width**2)) * np.exp(1j * k0 * x)
        self.initial_state = psi / np.sqrt(np.sum(np.abs(psi)**2))

    
    def evolve_state(self):
        psi = self.initial_state.copy()
        self.probabilities[0] = np.abs(psi)**2
        for i in range(self.steps):
            psi = self.evolution_operator @ psi
            self.probabilities[i+1] = np.abs(psi)**2

    
    def simulate_with_qiskit(self):
        self.evolve_state()
        fig, ax = plt.subplots(figsize=(10, 6))
        barrier_height_viz = self.barrier_height / 2
        barrier_x = np.linspace(self.barrier_start, self.barrier_end, 100)
        barrier_y = np.ones(100) * barrier_height_viz
        line, = ax.plot(self.x_range, self.probabilities[0], 'b-', lw=2)
        barrier_plot, = ax.plot(barrier_x, barrier_y, 'r-', lw=3)
        ax.fill_between(barrier_x, 0, barrier_y, color='red', alpha=0.3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, max(np.max(self.probabilities), barrier_height_viz) * 1.1)
        ax.set_xlabel('Position')
        ax.set_ylabel('Probability / Potential')
        ax.set_title('Quantum Tunneling Simulation')
        ax.legend(['Probability', 'Potential Barrier'], loc='upper right')
        time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes)

        
        def update(frame):
            line.set_ydata(self.probabilities[frame])
            time_text.set_text(f'Time: {frame * self.dt:.1f}')
            return line, time_text

        
        ani = FuncAnimation(fig, update, frames=self.steps+1, interval=50, blit=True)
        self.qiskit_qpe_implementation()
        return ani

    
    def qiskit_qpe_implementation(self):
        u_gate = UnitaryGate(self.evolution_operator)
        n_counting = 3
        qpe = QuantumCircuit(n_counting + self.num_qubits)
        initial_state_prep = self.create_qiskit_state_preparation()
        qpe = qpe.compose(initial_state_prep, list(range(self.num_qubits)))
        
        for qubit in range(n_counting):
            qpe.h(qubit)
            
        for i in range(n_counting):
            power = 2**i
            
            for _ in range(power):
                qpe.append(u_gate.control(1), [i] + list(range(n_counting, n_counting + self.num_qubits)))
                
        qpe.append(QFT(n_counting).inverse(), range(n_counting))
        qpe.measure_all()
        simulator = Aer.get_backend('qasm_simulator')
        compiled_circuit = transpile(qpe, simulator)
        job = simulator.run(compiled_circuit, shots=1024)
        result = job.result() # default return by qiskit
        counts = result.get_counts()
        plt.figure(figsize=(10, 6))
        plot_histogram(counts)
        plt.title('QPE Results: Eigenvalues of Time Evolution Operator')
        plt.xlabel('Measured Phase')
        plt.ylabel('Counts')
        plt.tight_layout()

    def create_qiskit_state_preparation(self):
        qc = QuantumCircuit(self.num_qubits)
        center_qubit = int(self.num_qubits * 0.2)
        
        for i in range(center_qubit - 1, center_qubit + 2):
            if 0 <= i < self.num_qubits:
                qc.h(i)
                
        for i in range(self.num_qubits):
            qc.p(i * np.pi / 8, i)
            
        return qc

    def visualize_barrier_transmission(self):
        left_of_barrier = np.where(self.x_range < self.barrier_start)[0]
        right_of_barrier = np.where(self.x_range > self.barrier_end)[0]
        prob_left = np.sum(self.probabilities[:, left_of_barrier], axis=1)
        prob_right = np.sum(self.probabilities[:, right_of_barrier], axis=1)
        time_points = np.linspace(0, self.total_time, self.steps + 1)
        plt.figure(figsize=(10, 6))
        plt.plot(time_points, prob_left, 'b-', label='Left of Barrier')
        plt.plot(time_points, prob_right, 'g-', label='Right of Barrier')
        plt.xlabel('Time')
        plt.ylabel('Probability')
        plt.title('Tunneling Probability wrt Time')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()

def run_simulation():
    simulator = QuantumTunnelingSimulator(num_qubits=5, dt=0.1, total_time=40.0)
    barrier_width = 0.5
    barrier_start = 0.2
    barrier_end = barrier_start + barrier_width
    barrier_height = 1.0
    simulator.setup_barrier(barrier_start, barrier_end, barrier_height)
    simulator.create_hamiltonian(wavepacket_energy=0.5)
    simulator.create_initial_wavepacket(center=0.0, width=0.05, k0=30)
    animation = simulator.simulate_with_qiskit()
    simulator.visualize_barrier_transmission()
    plt.show()
    return animation

if __name__ == "__main__":
    animation = run_simulation()
