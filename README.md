This repository showcases different python and c++ projects I have done during my studies.

## Ising model
I wrote a numerical simulation of 2D Ising model using C++. I independently learned and applied the ROOT library for data analysis and visualization. The generated figures are included in the same folder as the main code.
- Implemented Metropolis algorithm with periodic boundary conditions to compute thermodynamical quantities such as energy, magnetization, susceptibility and specific heat
- Used ROOT to visualize the results as a functions of temperature. In addition, I used ROOT to fit a function to find the critical temperature for the model
- Used ROOT to generate snapshots of the spin-lattice near critical temperature to visualize phase transition

## Pure dephasing model
In this project I studied open quantum system dynamics using pure dephasing model and used machine learning to identify environmental spectral types in Jupyter Notebook.
- Generated time-dependent coherence signals for different bath types with different parameters
- Converted the signals to frequency domain using FFT
- Built a PyTorch neural-Network to classify bath types

## Structure factor calculation
An analysis tool used to compute structure factor from a data obtained from cp2k ab initio molecular dynamics simulations.
Developed during a summer internship.
- Performs preprocessing on the data
- Uses numba to speed up calculation
- Compares simulation data to experimental data and visualizes the results.

## Image classification
An image classification project using machine learning. Images are read from urls and preprocessed for an alaysis.
- Ten features that capture color and texture information are extracted from the images
- Data is visualized using pairplot, histogram and PCA.
- Ridge, Random Forest and MLP classifiers were trained and evaluated to classify the images

## Superdense coding project
This project implements the superdense coding quantum computing algorithm using IBM qiskit in python. The project is written in jupyter notebook and includes a theoretical introduction, circuit implementation, simulation and conclusions quantum.
- Constructed quantum circuits for the superdense coding protocol which includes creating a bell state, encoding, decoding, measurement and circuit visualization with Qiskit.
- Showed how an eavesdropper's measurement collapses the entangled state, preventing interception of the message.
- Extended the protocol by including a third qubit to the ciruit based on a theoretical work. Analysed the results and showed that three bits of information can be transferred with two qubits.

## Courses_db
In this Python project I built a database which manages students, courses, teachers, study groups and credits.
- Implemented a relational SQLite Database and developed SQL queries using JOINs, GROUP BY, UNION
- Additional SQL queries to analyze, for example grade distribution

## yatzy-game
A complete yatzy game written in Python. The code is written in finnish.
- The game works on command-line interface
- Writes and reads data from a file
