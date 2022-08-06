# MC_stroke
Computer implementation of a mathematical model of ischemic stroke based on a cellular automaton using the Monte Carlo method

Authors: Ananiev A.V., Fursov A.V., Zinchenko D.I.

Computer program:
Module for Creating the Model of Evolution of States of the Brain Tissue Area during Brain Stroke
Abstract: The program is intended for modeling the evolution of the states of cerebral tissue during ischemic cerebral stroke by Monte Carlo method.  The main opportunities of the program are the following: simulation of states of brain tissue area of given size at different input parameters; real-time visualization of model states; visualization of quantitative relations of lattice nodes in each of the states in the form of diagrams; recording on hard disk evolution of model states as csv files, animation in gif format and video in mp4 format.
Type of computer: IBM PC - compatible PC
Programming language: Python 3.8
OS: Windows 10, Linux, Android
Program size:		19 Kb


Getting Started
1) Install project to your PC: git clone https://github.com/alex625051/MC_stroke.git
2) Install Python 3.8: sudo apt-get python3.8
3) Open project's folder cd MC_stroke
4) Update python pip: python3.8 -m pip install pip
5) Install virtual environment: sudo apt-get install virtualenv
6) Create virtual environment: virtualenv --python="/usr/bin/python3.8" "venv/"
7) Activate it: source venv/bin/activate
8) Install python's packegaes: pip install -r requiements.txt
9) Run modeling process: python mc_stroke.py
