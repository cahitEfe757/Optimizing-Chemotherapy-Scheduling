# Optimizing Chemotherapy Scheduling: An In Silico Comparison of MTD, Metronomic, and Synergistic Regimens Using Gompertzian Dynamics

This repository contains mathematical models simulating tumor growth dynamics to evaluate and compare the efficacy of Maximum Tolerated Dose (MTD) chemotherapy and Metronomic Chemotherapy (MCT). The project explores monotherapies as well as sequential, intermittent, and concurrent combination strategies.

## Overview

Traditional chemotherapy relies on MTD, which can lead to severe toxicity and drug resistance. MCT, characterized by frequent administration of low drug doses, aims to overcome these limitations. This computational study utilizes Ordinary Differential Equations (ODEs) to model tumor volume changes over time, integrating pharmacokinetic (PK) elimination and pharmacodynamic (PD) effects based on the Emax (Hill) equation.

## Mathematical Model

The core of the simulation is based on logistic/Gompertzian tumor growth limited by a carrying capacity ($K$), coupled with the drug's effect. 

The primary differential equation governing the tumor volume ($V$) is:

$$\frac{dV}{dt} = r V \left(1 - \frac{V}{K}\right) - (E \cdot res \cdot V)$$

Where:
* $r$: Intrinsic tumor growth rate
* $K$: Maximum tumor carrying capacity
* $E$: Drug effect (Pharmacodynamics)
* $res$: Resistance factor (Barrier logic for monotherapy vs. combination)

The pharmacodynamic effect ($E$) is modeled using the Hill equation:

$$E = \frac{\alpha C^n}{K_d^n + C^n}$$

Where $C$ is the drug concentration, $\alpha$ is the maximum effect, $K_d$ is the dissociation constant, and $n$ is the Hill coefficient.

## Repository Contents

* **`MTD_vs_MCT.py`**: Simulates and visualizes the failure of monotherapies. It models the 21-day cycle of MTD (relapse cycle) versus the daily continuous low-dose of MCT (dormancy plateau).
* **`Combination.py`**: Explores the synergy of combining both regimens. It models three distinct scenarios to breach the Immunological Cure Threshold ($V_{immune}$):
    1.  Sequential Therapy (MTD followed by MCT maintenance)
    2.  Intermittent Therapy (Chemo-Switch)
    3.  Concurrent Therapy (Simultaneous administration)

## Biological Thresholds Investigated

* **Dormancy Threshold ($V_{dormancy} = 200 \ mm^3$):** A barrier that monotherapies typically fail to breach due to resistance or insufficient sustained dosing.
* **Immunological Cure Threshold ($V_{immune} = 50 \ mm^3$):** The target volume below which the immune system is assumed to clear the remaining tumor cells, resulting in a cure.

## Requirements

To run the simulations, you will need Python installed along with the following libraries:

* `numpy`
* `matplotlib`

You can install the dependencies via pip:

```bash
pip install numpy matplotlib
