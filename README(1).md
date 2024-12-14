# THz Pulse Propagation Simulation

This project provides a Python-based framework for simulating and analyzing the propagation of terahertz (THz) electromagnetic pulses through various media. The simulations include both linear effects (diffraction and group velocity dispersion, GVD) and nonlinear effects (Kerr effect) using the Unidirectional Pulse Propagation Equation (UPPE).

## Project Structure

The project consists of the following key files:

### 1. `main.py`
- **Description**: Defines the core simulation class `THzProp`, which models the propagation of a THz pulse.
- **Features**:
  - Implements UPPE for linear and nonlinear effects.
  - Allows parameterized propagation of THz pulses.
  - Includes visualization capabilities for the evolution of the pulse in both space and time domains.
- **Usage**: This file can be run independently to test the propagation for a single set of parameters.

### 2. `calculation.py`
- **Description**: Contains scripts to run parameter sweeps and save results.
- **Features**:
  - Varies GVD coefficient (`beta2`), nonlinear index (`n2`), and wavelength (`lambda`).
  - Saves results (e.g., spatial and temporal maps, phase shifts) as `.npy` files.
- **Usage**: Run this file to generate datasets for further analysis.

### 3. `analysis.py`
- **Description**: Provides tools for visualizing and analyzing the simulation results.
- **Features**:
  - Loads saved data from `calculation.py`.
  - Produces detailed plots for time/frequency domain evolution, spatial and temporal propagation, and parameter sweeps.
- **Usage**: Run this file to visualize and analyze the generated datasets.

## Getting Started

### Prerequisites
- Python 3.7+
- Required Python packages:
  - `numpy`
  - `matplotlib`
  - `tqdm`

Install dependencies using:
```bash
pip install numpy matplotlib tqdm
```

### Running the Project

1. **Single Simulation**:
   - Modify the parameters in `main.py` as needed.
   - Run `main.py`:
     ```bash
     python main.py
     ```

2. **Parameter Sweeps**:
   - Edit `calculation.py` to adjust the ranges of `beta2`, `n2`, and `lambda`.
   - Run `calculation.py`:
     ```bash
     python calculation.py
     ```

3. **Visualization**:
   - Use `analysis.py` to load and visualize results.
   - Run `analysis.py`:
     ```bash
     python analysis.py
     ```

## Results
The project saves results as `.npy` files in the `result/` directory. Subdirectories are created for each parameter sweep, storing data for spatial maps, temporal maps, and phase shifts. 

### Example Visualizations
- **Temporal Evolution**: Tracks the pulse in the time domain as it propagates.
- **Spatial Evolution**: Visualizes the beam profile changes in the spatial domain.
- **Parameter Sweeps**: Compares dynamics for varying `beta2`, `n2`, and `lambda`.

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

## Acknowledgements
This project was developed to model advanced optical pulse propagation phenomena, with applications in THz spectroscopy and nonlinear optics research.

---

For further questions or contributions, feel free to create an issue or submit a pull request on the GitHub repository.
