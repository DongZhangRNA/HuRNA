## HuRNA: A Hybrid framework for universal RNA 3D Structure Prediction
Contributed by Deyin Wang, LinXi Zhang, Ruhong Zhou, and Dong Zhang
Version 1.0, 30/10/2026

********************************************************
HuRNA is a coarse-grained model for de novo prediction RNA 3D structures using motif based template searching. The model has four main features: a template database for RNA motifs, a five-bead coarse-grained representation per nucleotide, a coarse-grained force field extracted through rigorous reference state simulations, and replica-exchange molecular dynamics simulations for conformational sampling. The main program is "Predict_RNA_structure": To fold tertiary structures of RNAs with complicated topologies with or without given secondary structure constraints; 

********************************************************
#### Requirements:
Linux operating systems
>=10 threads
RAM of >= 5G
g++ compiler with C++14 support

********************************************************
#### Dependencies:
HuRNA relies on MPICH platform for parallel REMD simulations (10 replicas) in LAMMPS and AmberTools22 to configure energy minimization in Amber force field.

********************************************************
#### Installation:
1. Install an MPI implementation providing mpirun, such as mpich >=3.2 (http://www.mpich.org/downloads/ or sudo apt-get install mpich if you are using Ubuntu)
2. Install AmberTools for all-atom energy minimization, optional but strongly recommended (https://ambermd.org/GetAmber.php#ambertools)
3. Install HuRNA: unzip HuRNA.zip; change to the project root; run `bash ./install.sh`

The installer checks for cmake, make, gcc, g++, and mpirun on PATH, and verifies that g++ can compile C++14 before starting a build. It exits with an error if a prerequisite is missing or the compilation check fails.

Run the installer with `bash ./install.sh` as shown above. If it is sourced, it
runs in a subshell so its exits, working-directory changes, shell options, and
traps do not change the calling shell. A failed installation still returns a
nonzero status; callers using `set -e` must handle that status to continue.

********************************************************
#### Usage descriptions:

(A) For structure prediction:

`Predict_RNA_structure`:
```bash
  -i, input file contains sequence (required) and 2D structure (optional)
  -name, jobname, default "rna"
  -DIR, work directory, default "${jobname}_prediction"
  -exclude, file contained the excluded PDBids when finding template
  -restr, restraint file contained the user_template, distance_restraints, rigid_body, etc
  -nstep, steps for production CGMD simulation, default "-1" 
  -version, CG force field version, default "v2"
  -pe, percentage of top lowest potential energy structrues for structure prediction, default "0.1"
  -cutoff, rmsd cutoff (in Angstrom) for structure clustering, default "-1.0"
  -pc, percentage of structures in largest cluster, default "0.5" (-cutoff and -pc are mutually exclusive)
  -exe, execute the prediction (1) or not (0), default "1"
  -template_rigid, treat the user_template in restraint file as individual rigid bodies, no value.
  -log, log file, default "${jobname}.log"
  -homedir, path to HuRNA, if not specified, automatically parse
```

(B) For isolated main components of HuRNA:

`Construct_RNA_initial_structure`:
```bash
  -i, input file contains sequence (required) and 2D structure (required)
  -seq, input sequence
  -secstr, input 2D structure
  -name, jobname, default "rna"
  -DIR, work directory, default "prediction"
  -restr, restraint file contained the user_template, distance_restraints, rigid_body, etc
  -nstep, steps for CGMD simulation to assemble the initial structure, default "-1"
  -nmax_template, maximum number of templates for each motif 
  -log, log file, default "${jobname}.log"
  -homedir, path to HuRNA, if not specified, automatically parse
```
  
`Perform_CGMD_simulation`:
```bash
  -i, input file contains sequence (required) and 2D structure (required)
  -if, input initial RNA structure file (required)
  -iformat, the format of input initial RNA structure file ('aa' or 'cg'), default 'cg'
  -restr, restraint file contained the user_template, distance_restraints, rigid_body, etc
  -DIR, work directory, default "prediction"
  -name, jobname, default "rna"
  -log, log file, default "${jobname}.log"
  -ti_fold, refolding the teritary helices related to the nested base-pairing interactions, no value.
  -ss_restr, using secondary structure as restraints, no value
  -version, CG force field version, default "v2"
```

`Obtain_top_structure`:
```bash
  -lt, file contained the lammps trj files list (LAMMPS MD trajectories),should be consistent with the order of log files, default '${DIR}/${jobname}_prod_lammpstrj_list.dat'
  -ll, file contained the lammps log files list (LAMMPS log files for potential energies), should be consistent with the order of trj files, default '${DIR}/${jobname}_prod_lammpslog_list.dat'
  -ref, reference file contained the residue information for lammps trj files, default '${DIR}/${jobname}_cg.pdb'
  -log, log file, default "${jobname}.log"
  -o, output file, default "${DIR}/${jobname}_predicted_cg.pdb"
  -b, step of first frame to read from trajectory, default 0
  -e, step of last frame to read from trajectory, default -1
  -pe, percentage of top lowest potential energy structrues for structure prediction, default "0.1"
  -cutoff, rmsd cutoff (in Angstrom) for structure clustering, default "-1.0"
  -pc, percentage of structures in largest cluster, default "0.5" (-cutoff and -pc are mutually exclusive)
```

`Reconstruct_allatom_structure`:
```bash
  -i, input file contains the coarse-grained RNA structure
  -DIR, work directory, default "prediction"
  -name, jobname, default "rna"
  -o, output file contains predicted all-atom structure, default "${DIR}/${jobname}_predicted.pdb"
  -log, log file, default "${jobname}.log"
  -min, perform all-atom structure energy minimization, no value
  -version, CG force field version, default "v2"
```

********************************************************
### Example usage

The [example](example/) folder contains two PU31 prediction workflows. Complete
the installation above first, including the MPI and VMD dependencies, and
check that `bin/Predict_RNA_structure` has been built. Run each script from
its own example directory: the paths in the scripts are relative to that
directory.

#### End-to-end example

From the project root, run:

```bash
cd example/predict_end-to-end
bash predict_pu31_seq.sh
```

The script uses [`pu31.seq`](example/predict_end-to-end/pu31.seq) as input and
[`pu31_exclude.dat`](example/predict_end-to-end/pu31_exclude.dat) as the
template exclusion list. It creates and enters `pu31_predict/` before starting
the prediction.

#### Example with a supplied secondary structure

From the project root, run:

```bash
cd example/predict_with_2D
bash predict_pu31.sh
```

This script uses [`pu31.2d`](example/predict_with_2D/pu31.2d) and the
[`pu31_exclude.dat`](example/predict_with_2D/pu31_exclude.dat) in the same
directory. It also creates and enters its own `pu31_predict/` directory.

#### Inputs and options

Both supplied input files currently contain the same two lines: the RNA
sequence followed by its secondary structure in dot-bracket notation. Dots
represent unpaired bases; parentheses and square brackets denote base pairs.
Despite the `.seq` name and the `predict_end-to-end` directory name, the
supplied `pu31.seq` already includes the secondary structure. Consequently,
both scripts use the supplied structure rather than predicting it from the
sequence alone.

Both scripts invoke `Predict_RNA_structure` with:

- `-i`: the input file, relative to `pu31_predict/`.
- `-name pu31`: use `pu31` as the job name and output prefix.
- `-exclude ../pu31_exclude.dat`: exclude PDB ID `7mlx` from template searches.
- `-exe 1`: execute the prediction rather than only generate the pipeline script.

The scripts do not set `-DIR`, so the program uses its default working
directory, `pu31_prediction/`, inside `pu31_predict/`.

#### Optional: predict from the sequence alone

To run a genuinely sequence-only prediction without changing the supplied
files, create a new working directory and copy just the sequence line. From
the project root:

```bash
cd example/predict_end-to-end
mkdir pu31_seq_only && \
  head -n 1 pu31.seq > pu31_seq_only/pu31.seq && \
  cd pu31_seq_only && \
  ../../../bin/Predict_RNA_structure \
    -i pu31.seq -name pu31 -exclude ../pu31_exclude.dat -exe 1
```

With no second input line, `Predict_RNA_structure` uses CONTRAfold by default
to predict the secondary structure before running the 3D prediction. The
command above requires a fresh `pu31_seq_only/` directory.

#### Optional: predict with user-specified restraints

To run a genuinely sequence-only prediction without changing the supplied
files, create a new working directory and copy just the sequence line. From
the project root:

```bash
cd example/predict_with_restraints
bash predict_pu02.sh
```

With the `-restr pu02_restr.dat` option, `Predict_RNA_structure` uses specified
restraints to predict the secondary structure before running the 3D prediction.

#### Output and rerunning

For either bundled script, look in its `pu31_predict/` directory for:

- `pu31.log`: the prediction log.
- `predict_pu31.sh`: the generated pipeline script.
- `pu31_prediction/pu31_predicted_cg.pdb`: the selected coarse-grained model.
- `pu31_prediction/pu31_predicted.pdb`: the reconstructed all-atom model.
- `pu31_prediction/pu31_predicted_min.pdb`: the minimized model, if produced.

When a final model is available, the pipeline copies it into
`pu31_predict/`: it prefers `pu31_predicted_min.pdb`, falling back to
`pu31_predicted.pdb` if no minimized model is present. The optional
sequence-only command uses the same layout under `pu31_seq_only/`.

The bundled scripts use `mkdir pu31_predict` and do not resume an existing
directory. Before rerunning an example, move or rename the previous
`pu31_predict/` directory, or use a fresh copy of the example.


********************************************************
Citation:
Wang et al. HuRNA: A Hybrid framework for universal RNA 3D Structure Prediction

Please contact us at zhangd_iqb@zju.edu.cn for more information.
