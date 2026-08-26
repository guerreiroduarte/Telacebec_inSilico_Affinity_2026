#!/usr/bin/env python3

from pathlib import Path

import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align

simulation_dir = Path("data/raw_data/simulation")
top_dir = Path("data/raw_data/amber_inputs")

rmsd_results = []
rmsf_results = []
for traj in simulation_dir.glob("*.dcd"):
	terms = traj.stem.split("_")

	organism = terms[0]
	replica = terms[-1]

	top = top_dir / f"{organism}_complex.prmtop"

	u = mda.Universe(top, traj)

	# Cálculo do RMSD
	rmsd_run = rms.RMSD(u, u, ref_frame=0, select="backbone and resname HUU").run()

	rmsd_data = rmsd_run.results.rmsd
	for frame, time, rmsd in rmsd_data:
		rmsd_results.append({
			"Organism": organism,
			"Replica": replica,
			"Frame": str(frame),
			"Time_ns": time,
			"RMSD": rmsd
		})

	# Cálculo do RMSF
	average = align.AverageStructure(u, u, select = "protein and name CA", ref_frame=0).run()
	ref = average.results.universe

	aligner = align.AlignTraj(u, ref, select = "protein and name CA", in_memory=True).run()

	c_alphas = u.select_atoms("protein and name CA")
	rmsf_run = rms.RMSF(c_alphas).run()

	rmsf_data = rmsf_run.results.rmsf
	for resid, rmsf_value in zip(c_alphas.resids, rmsf_data):
		rmsf_results.append({
			"Organism": organism,
			"Replica": replica,
			"Residues": resid,
			"RMSF": rmsf_value
		})

rmsd_df = pd.DataFrame(rmsd_results)
rmsd_df.to_csv(simulation_dir / "rmsd.csv", index=False)

rmsf_df = pd.DataFrame(rmsf_results)
rmsf_df.to_csv(simulation_dir / "rmsf.csv", index=False)
