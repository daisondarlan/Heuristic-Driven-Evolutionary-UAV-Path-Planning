# Evolutionary Algorithm with Domain-Specific Operators for UAV Path Planning

This repository contains the MATLAB source code associated with the paper:  
**Evolutionary Algorithm with Domain-Specific Operators for UAV Path Planning**  
*Authors: Daison Darlan, Oladayo S. Ajani, Rammohan Mallipeddi*

The paper proposes a multi-objective evolutionary framework for UAV path planning that incorporates domain-specific operators tailored to the structure of the problem. In particular, the method integrates A*-guided crossover and mutation operators within the MOEA/D-AWA framework, together with an adaptive polynomial mutation strategy, to generate feasible and efficient UAV trajectories in complex environments.

This repository provides the MATLAB code for:

- **A*-guided crossover operator:** A domain-specific recombination strategy that connects selected parent waypoints using A* search to generate feasible offspring paths.
- **A*-guided mutation operator:** A local refinement operator that replaces inefficient or infeasible path segments through A*-based replanning.
- **Adaptive polynomial mutation:** A generation-dependent mutation mechanism that balances exploration and exploitation during evolution.
- **MOEA/D-AWA-based optimization framework:** The overall multi-objective evolutionary algorithm used to optimize UAV trajectories.

## Repository Structure

- `MOEAD_AWA_Astar/`  
  Contains the main implementation of the proposed algorithm.
  - `MOEADAWA_Astar.m` — main file to run the algorithm
  - `aStarCrossover.m` — implementation of the A*-guided crossover operator
  - `aStarMutation.m` — implementation of the A*-guided mutation operator

- `problems/`  
  Contains the benchmark environments used for the path-planning experiments.

For more information about the benchmark environments and their generation, please refer to the companion repository:  
[UAV-Path-Planning-Benchmark](https://github.com/daisondarlan/UAV-Path-Planning-Benchmark)

## Getting Started

### Prerequisites
- MATLAB R2019b or later
- Required MATLAB toolboxes

### Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/daisondarlan/Heuristic-Driven-Evolutionary-UAV-Path-Planning.git
