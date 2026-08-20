#!/usr/bin/env python3

from pathlib import Path
from sys import stdout

from openmm import *
from openmm.app import *
from openmm.unit import *

input_dir = Path("data/raw_data/amber_inputs/")
simulation_dir = Path("data/raw_data/simulation")
simulation_dir.mkdir(exist_ok=True)

platform = Platform.getPlatformByName("CUDA")
properties = {"Precision": "mixed"}

for prmtop_file in input_dir.glob("*.prmtop"):
    
    complex_name = prmtop_file.stem
    inpcrd_file = input_dir / f"{complex_name}.inpcrd"

    prmtop = AmberPrmtopFile(prmtop_file)
    inpcrd = AmberInpcrdFile(inpcrd_file)

    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=1.0*nanometer,
        constraints=HBonds,
        hydrogenMass=4*amu
    )

    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 4*femtosecond)
    system.addForce(MonteCarloBarostat(1.0*bar, 300*kelvin, 25))

    simulation = Simulation(
        prmtop.topology, 
        system, 
        integrator, 
        platform, 
        properties
    )
    simulation.context.setPositions(inpcrd.positions)

    print("Minimizando energia")
    simulation.minimizeEnergy(maxIterations=1000)

    print("Equilíbrio NPT")
    npt_reporter = StateDataReporter(
        stdout,
        1000,  # Imprime no terminal a cada 1000 passos (2 ps)
        step=True,
        potentialEnergy=True,
        temperature=True,
        density=True,
        speed=True,
        separator='\t'
    )
    simulation.reporters.append(npt_reporter)
    simulation.context.setVelocitiesToTemperature(300*kelvin)
    simulation.step(15000)
    simulation.reporters.remove(npt_reporter)

    sim_dcd = simulation_dir / f"{complex_name}.dcd"

    simulation.reporters.append(DCDReporter(sim_dcd, 500))
    simulation.reporters.append(
        StateDataReporter(
            stdout,
            500,
            step=True,
            potentialEnergy=True,
            temperature=True,
            density=True,
            speed=True,
            separator='\t'
        )
    )

    simulation.step(100000)
    print("Simulação concluída!")