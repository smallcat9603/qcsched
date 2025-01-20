# QC-HPC Co-Scheduler Simulation

This repository is a reference implementation of Quantum-HPC Co-Scheduler simulation.
The simulation supports four job scheduling algorithms: FCFS (First Come First Serve), SJF (Shortest Job First), Priority, and QPriority (default).

## Requirements

```
streamlit
matplotlib
```

## Usage

```
streamlit run app.py
```

## Input File Format

CSV files can be used for job submission. An example of the input file format is as follows.

```
# HPC Selection, Job Type, HPC Nodes, Elapsed Time, Start Time, Priority
HPC1,HPC,4,5,0,2 
HPC2,QC1,5,3,0,3
HPC3,QC2,3,6,0,1
```
