# Satcolavoidance
LEO satellite conjunction detection and avoidance simulation in MATLAB
# Satellite Conjunction Avoidance Simulation
**MATLAB-based LEO satellite collision avoidance system with ML risk prediction**

## Project Overview
Simulates a 12-satellite Earth-imaging constellation in LEO alongside 
a synthetic debris field. Detects conjunction risks, plans avoidance 
maneuvers, tracks fuel consumption, and uses machine learning to 
predict collision risk from orbital elements.

## Features
- Custom Two-Body + J2 orbital propagator 
- 12-satellite sun-synchronous constellation (3 planes × 4 sats)
- 150 synthetic debris objects with realistic size distribution
- Dual-trigger conjunction detection: miss distance + Chan Pc estimate
- Automated avoidance maneuver planner with fuel budget tracking
- ML risk classifier (Safe / Warning / Critical) + ΔV recommender
- Interactive app: input any satellite/debris data and get instant prediction

## How to Run
**Full simulation:**
```matlab
run_validation      % run tests first
main_simulation     % runs the full 5-day simulation
```
**ML interactive app:**
```matlab
train_conjunction_ml(3000)   % train model once
ml_conjunction_app            % launch interactive predictor
```

## File Structure
| File | Purpose |
|------|---------|
| `main_simulation.m` | Entry point — full simulation |
| `generate_satellites.m` | 12-sat LEO constellation |
| `generate_debris.m` | Synthetic debris field |
| `propagate_all_orbits.m` | J2 orbital propagator |
| `detect_and_avoid.m` | Conjunction detection + maneuver planner |
| `fuel_analysis.m` | Fuel budget + lifetime estimator |
| `visualize_simulation.m` | All output plots |
| `run_validation.m` | 7-test validation suite |
| `generate_ml_training_data.m` | ML training data generator |
| `train_conjunction_ml.m` | ML model trainer |
| `ml_conjunction_app.m` | Interactive prediction app |
