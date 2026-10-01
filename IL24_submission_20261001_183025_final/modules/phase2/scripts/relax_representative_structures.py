#!/usr/bin/env python3
"""Restrained OpenMM minimization for shortlisted Boltz complexes."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from openmm import CustomExternalForce, LangevinMiddleIntegrator, Platform, unit
from openmm.app import ForceField, HBonds, NoCutoff, PDBFile, Simulation
from pdbfixer import PDBFixer


def relax(source: Path, destination: Path) -> dict[str, object]:
    fixer = PDBFixer(filename=str(source))
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)

    forcefield = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    system = forcefield.createSystem(
        fixer.topology,
        nonbondedMethod=NoCutoff,
        constraints=HBonds,
    )
    restraint = CustomExternalForce(
        "0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)"
    )
    restraint.addGlobalParameter("k", 1000.0 * unit.kilojoule_per_mole / unit.nanometer**2)
    for parameter in ("x0", "y0", "z0"):
        restraint.addPerParticleParameter(parameter)
    backbone = {"N", "CA", "C", "O"}
    positions = fixer.positions
    restrained = 0
    for atom in fixer.topology.atoms():
        if atom.name in backbone:
            xyz = positions[atom.index].value_in_unit(unit.nanometer)
            restraint.addParticle(atom.index, xyz)
            restrained += 1
    system.addForce(restraint)

    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin,
        1 / unit.picosecond,
        0.002 * unit.picoseconds,
    )
    platform = Platform.getPlatformByName(os.environ.get("OPENMM_PLATFORM", "CPU"))
    properties = {"CudaPrecision": "mixed"} if platform.getName() == "CUDA" else {}
    simulation = Simulation(fixer.topology, system, integrator, platform, properties)
    simulation.context.setPositions(positions)
    initial = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    max_iterations = int(os.environ.get("RELAX_MAX_ITERATIONS", "500"))
    simulation.minimizeEnergy(
        tolerance=100.0 * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=max_iterations,
    )
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    final = final_state.getPotentialEnergy()

    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w") as handle:
        PDBFile.writeFile(
            fixer.topology,
            final_state.getPositions(),
            handle,
            keepIds=True,
        )
    return {
        "source": str(source),
        "destination": str(destination),
        "restrained_backbone_atoms": restrained,
        "initial_potential_energy_kj_mol": float(initial.value_in_unit(unit.kilojoule_per_mole)),
        "final_potential_energy_kj_mol": float(final.value_in_unit(unit.kilojoule_per_mole)),
        "platform": platform.getName(),
        "forcefield": "amber14-all.xml + implicit/gbn2.xml",
        "backbone_restraint_k_kj_mol_nm2": 1000.0,
        "minimization_tolerance_kj_mol_nm": 100.0,
        "max_iterations": max_iterations,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = relax(args.source, args.destination)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
