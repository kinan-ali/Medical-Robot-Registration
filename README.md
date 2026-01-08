## Project Description
This repository contains the source code for the "Medical Robot Registration" project, part of the Healthtech Master track (2025-2026).

The goal of this project is to simulate and solve the robotic registration problem for a computer-assisted surgical intervention. The specific medical task is to automatically position a Cartesian robot holding a straight needle to reach a target inside a patient's abdominal cavity, under the guidance of an endoscopic camera.

## System Architecture
The simulation models a realistic surgical setup involving:
* **Cartesian Robot:** Moves the instrument (needle) with a passive wrist.
* **Endoscopic Camera:** Manually positioned, providing visual feedback.
* **Optical Navigation System:** Tracks rigid bodies attached to the instrument and camera.
* **Geometric Constraints:** The instrument is constrained to pass through a fixed incision point (trocart) on the abdominal wall.

## Project Structure

### 1. Registration with Optical Navigation
Implementation of the registration procedure using an external optical tracker. This module computes the necessary frame transformations between the robot, the camera, and the patient to guide the needle accurately.

### 2. Registration without Navigation (Vision-Based)
An alternative approach that removes the dependency on the optical tracker. This module estimates the instrument pose relative to the camera using markers on the instrument body, exploring the trade-offs between system complexity and registration accuracy.

### 3. Error Propagation Analysis
A study of how physiological movements affect precision.
* **Scenario:** Modeling uncertainty in the trocart (incision point) position (standard deviation of 1mm).
* **Methods:** Evaluation using both variance propagation (analytical) and numerical simulation (Monte Carlo).
* **Output:** Computation and visualization of error ellipsoids at the distal end of the instrument.

### 4. Visual Servoing
Implementation of a closed-loop control system to compensate for calibration errors and patient movement.
* **Technique:** Position-Based Visual Servoing (PBVS).
* **Control Law:** Gauss-Newton optimization to minimize the error between the current and desired instrument position using 3D reconstruction feedback.

## Setup and Usage

### Prerequisites
* **MATLAB** (Required for simulation and matrix manipulation).
* **Lab Files:** Ensure the provided simulation state-machine functions (e.g., `InitConfig`, `MoveCamera`) are in the MATLAB path.

### How to Run
1.  Clone this repository.
2.  Open MATLAB and navigate to the project folder.
3.  Initialize the simulation environment:
    ```matlab
    InitConfig();
    ```
4.  Run the main registration scripts or specific section modules.
5.  Use `DisplayConfig()` to visualize the current state of the robot and camera.

## Authors
* **Kinan ALI**
* **Filippo MORINI**

---
*Course: Healthtech Master Track - Robot Registration 2025-2026*
