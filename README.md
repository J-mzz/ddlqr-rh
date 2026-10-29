# Robust Data-Driven Receding-Horizon Control for LQR With Input Constraints

MATLAB code for the paper

> J. Zheng and M. Sznaier, "Robust data-driven receding-horizon control for LQR with input constraints," *IEEE Control Systems Letters*, vol. 9, pp. 2765–2770, 2025.
> [doi:10.1109/LCSYS.2025.3642051](https://doi.org/10.1109/LCSYS.2025.3642051) · [arXiv:2510.06117](https://arxiv.org/abs/2510.06117)

The method computes a robust LQR state feedback with input constraints for every linear system consistent with noisy data, by solving a semidefinite program obtained through duality. In receding horizon, every new measurement refines the set of consistent systems, which makes the controller less conservative. The scripts compare it with the data-driven LQR without input constraints of Dai and Sznaier, DeePC, and a model-based controller.

## Requirements

- MATLAB (tested with R2024a) with the Control System Toolbox, the Optimization Toolbox, and the Statistics and Machine Learning Toolbox
- [YALMIP](https://yalmip.github.io/)
- [MOSEK](https://www.mosek.com/) (free academic licenses are available); tested with MOSEK 9.1

## Usage

```matlab
ddlqr_2d                     % Example 1 (Fig. 1): 2 states, 2 inputs
ddlqr_4d                     % Example 2 (Fig. 2): 4 states, 4 inputs
```

Each script runs all methods for 10 closed-loop steps from the same initial state and noise sequence. It saves the results in `results/2d/` or `results/4d/` and plots them with `plot_results`.

The results reported in the paper are stored in `example_2d/` and `example_4d/` (with the figures as PNG files). To plot them without rerunning the examples:

```matlab
plot_results('example_2d')
plot_results('example_4d')
```

| Legend | Method | Script section |
|---|---|---|
| ours data-driven | receding-horizon input-constrained DD-LQR (Algorithm 2), consistency set refined at every step | Case 3 |
| ours nonRH data-driven | input-constrained DD-LQR, consistency set fixed | Case 2 |
| ours model-based | input-constrained LQR with known `A` and `B` (Algorithm 1) | Model-based |
| ddlqr | DD-LQR [10], without input constraints, designed once from the initial data | Case 1 |
| deepc | DeePC [11] | Case 0 |

Case 4 solves the receding-horizon problem of Case 3 with Scherer's P-satz instead of Putinar's, for the comparison of the two relaxations in the paper. It is off by default; set `run_scherer = true` in the script to run it. Its result for Example 1 is in `example_2d/ddlqr_sch.mat`.

The figures show the trajectories, the inputs (`||UU*u||_2` for several inputs), the worst-case upper bound on the cost, the Lyapunov function, and (2D only) the consistency sets with and without receding horizon, projected onto their first three principal components.

## Repository layout

```
ddlqr_2d.m, ddlqr_4d.m     scripts for the 2D and 4D examples
plot_results.m             comparison figures from a results folder
example_2d/, example_4d/   results and figures reported in the paper
utils/
  DDLQR_IC.m               input-constrained DD-LQR (Algorithm 2)
  DDLQR_scherer.m, psatz_matrix_refining.m    same with Scherer's P-satz (Case 4)
  DDLQR_D.m                data-driven LQR without input constraints
  MBLQR.m                  model-based input-constrained LQR
  deepc_solve.m            DeePC
  psatz_scalar*.m          Positivstellensatz certificates used by the SDPs
  generate_discrete.m, generate_noise.m    noisy data
  plot_consistency_sets.m  consistency-set figure
  third_party/             external code, see utils/third_party/README.md
```

## Citation

If you use this code, please cite:

```bibtex
@article{zheng2025robustlqr,
  title={Robust data-driven receding-horizon control for {LQR} with input constraints},
  author={Zheng, Jian and Sznaier, Mario},
  journal={IEEE Control Systems Letters},
  volume={9},
  pages={2765--2770},
  year={2025},
  publisher={IEEE},
  doi={10.1109/LCSYS.2025.3642051}
}
```

## License

The code in this repository is released under the [MIT License](LICENSE), except the files in `utils/third_party/`, which keep their original licenses.
