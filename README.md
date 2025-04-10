
# Artificial Elbow – Kinematic and Dynamic Analysis

This repository contains the full analysis and simulation of an artificial elbow joint, developed as part of the course *Applied Mechanics to Machines for Biomedical Engineering* (Academic Year 2020/2021). The project includes kinematic and dynamic modeling, motor control strategies, and vibration analysis, all implemented in MATLAB.

## 📑 Project Overview

The goal of this project is to model and simulate the movement of an artificial elbow using a multi-link mechanical system. Key areas of focus include:

- Kinematic modeling
- Dynamic force analysis
- Motor torque calculations
- Motion law simulation
- Vibration analysis and system response

## 📌 Contents

- **System Geometry & Parameters**: Description of the mechanical arm’s configuration and reference data.
- **Kinematic Scheme**: Link and joint positions, including velocity and acceleration vectors.
- **Kinematic Analysis**: Angular velocities and accelerations of arm components.
- **Motion Law**: Time-based simulation of motion for 0.25s.
- **Dynamic Analysis**: Torque requirements, constraint reactions, and external force effects.
- **Motor Drive**: Motor characteristic curve (torque vs angular speed) and ideal operating conditions.
- **Vibrations**: Study of natural frequency, damping, and system transfer function.

## 🧮 Technologies Used

- MATLAB: For simulations, numerical solving (`fsolve`), animations, and plotting.
- Mathematical modeling: For dynamic equilibrium and kinematic closure equations.

## 🔧 Setup & Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/artificial-elbow.git
   cd artificial-elbow
   ```

2. Open MATLAB and run the main simulation script:
   ```matlab
   run('artificial_elbow_matlab.m')
   ```

3. Explore the plots and animation generated for system motion and vibration response.

## 📁 Files

- `artificial_elbow_matlab.m` – matlab file with all simulations
- `FINAL REPORT ARTIFICIAL ELBOW.pdf` – Report of the project (in Italian)

## 👩‍🎓 Author

**Chiara Raineri**  

## 📜 License

This project is for academic and educational purposes only. 
For licensing terms, see the [LICENSE](LICENSE) file.
