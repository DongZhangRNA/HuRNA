CXX = g++
CXXFLAGS = -std=c++14 -O3 -Wall -pthread
HOME := $(dir $(patsubst %/,%,$(CURDIR)))
CIF = $(HOME)external/cifparse-obj-v7.105-prod-src
EIGEN = $(HOME)external/eigen-3.4.0
BIN := $(HOME)bin
LDFLAGS = $(CIF)/lib/all.a

all: $(BIN) combined RNA Fetch_PDB Process_PDB Edit_chain Replace_sequence Rotate_molecule Extract_RNA_secondary_structure \
	Update_RNA_template Construct_RNA_initial_structure Perform_CGMD_simulation Obtain_top_structure \
	Reconstruct_allatom_structure Predict_RNA_structure  Calculate_CGMD_rmsd Calculate_heavyatom_rmsd Calculate_TMscore \
	Calculate_RNA_INF Extract_CGMD_snapshot Build_RNA_helix Reset_lammps_data Build_RNA_motif_structure
	
$(BIN):
	mkdir -p $@
	
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
	  
%.o: %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@
	
read_PDB.o: read_PDB.cpp
	$(CXX) $(CXXFLAGS) -I $(CIF)/include -c $< -o $@

superpose.o: superpose.cpp
	$(CXX) $(CXXFLAGS) -I $(EIGEN) -c $< -o $@

Combined = cobmms.o error.o initialization.o atom.o residue.o chain.o molecule.o read_PDB.o superpose.o random.o sort_vector.o cluster.o
combined: $(Combined)
	ld -r -o combined.o $(Combined)

RNA = RNA_molecule.o RNA_basepair.o RNA_motif.o read_RNA_template.o RNA_coarse_grained.o RNA_recovery.o	measure.o
RNA: $(Combined) $(RNA)
	ld -r -o RNA.o $(Combined) $(RNA)

FetchPDB = Fetch_PDB.o	
Fetch_PDB: $(FetchPDB)
	$(CXX) $(CXXFLAGS) $(FetchPDB) -o $(BIN)/$@

ProcessPDB = Process_PDB.o RNA.o
Process_PDB: $(ProcessPDB)
	$(CXX) $(CXXFLAGS) $(ProcessPDB) $(LDFLAGS) -o $(BIN)/$@

EditChain = Edit_chain.o combined.o
Edit_chain: $(EditChain)
	$(CXX) $(CXXFLAGS) $(EditChain) $(LDFLAGS) -o $(BIN)/$@

RelpaceSequence = Replace_sequence.o combined.o RNA_molecule.o RNA_basepair.o
Replace_sequence: $(RelpaceSequence)
	$(CXX) $(CXXFLAGS) $(RelpaceSequence) $(LDFLAGS) -o $(BIN)/$@

RotateMolecule = Rotate_molecule.o combined.o
Rotate_molecule: $(RotateMolecule)
	$(CXX) $(CXXFLAGS) $(RotateMolecule) $(LDFLAGS) -o $(BIN)/$@

ExtractRNASecondaryStructure = Extract_RNA_secondary_structure.o combined.o RNA_molecule.o RNA_basepair.o
Extract_RNA_secondary_structure: $(ExtractRNASecondaryStructure)
	$(CXX) $(CXXFLAGS) $(ExtractRNASecondaryStructure) $(LDFLAGS) -o $(BIN)/$@

UpdateRNATemplate = Update_RNA_template.o RNA.o extract_RNA_template.o	
Update_RNA_template: $(UpdateRNATemplate)
	$(CXX) $(CXXFLAGS) $(UpdateRNATemplate) $(LDFLAGS) -o $(BIN)/$@

ConstructRNAInitialStructure = Construct_RNA_initial_structure.o RNA.o read_restraint.o parse_RNA_motif.o match_RNA_motif.o \
	find_RNA_template.o assemble_RNA_motif.o build_RNA_model.o read_amber.o read_lammps.o read_coeff.o prepare_lammps.o 
Construct_RNA_initial_structure: $(ConstructRNAInitialStructure)
	$(CXX) $(CXXFLAGS) $(ConstructRNAInitialStructure) $(LDFLAGS) -o $(BIN)/$@

PerformCGMDSimulation = Perform_CGMD_simulation.o combined.o measure.o read_restraint.o RNA_molecule.o RNA_coarse_grained.o RNA_basepair.o \
	RNA_motif.o parse_RNA_motif.o read_amber.o read_coeff.o prepare_lammps.o read_lammps.o
Perform_CGMD_simulation: $(PerformCGMDSimulation)
	$(CXX) $(CXXFLAGS) $(PerformCGMDSimulation) $(LDFLAGS) -o $(BIN)/$@
	
ObtainTopStructure = Obtain_top_structure.o combined.o read_lammps.o
Obtain_top_structure: $(ObtainTopStructure)
	$(CXX) $(CXXFLAGS) $(ObtainTopStructure) $(LDFLAGS) -o $(BIN)/$@
	
ReconstructAllatomStructure = Reconstruct_allatom_structure.o combined.o measure.o RNA_molecule.o RNA_basepair.o RNA_recovery.o RNA_coarse_grained.o \
	read_amber.o read_lammps.o read_coeff.o prepare_lammps.o 
Reconstruct_allatom_structure: $(ReconstructAllatomStructure)
	$(CXX) $(CXXFLAGS) $(ReconstructAllatomStructure) $(LDFLAGS) -o $(BIN)/$@

PredictRNAStructure = Predict_RNA_structure.o cobmms.o initialization.o error.o atom.o random.o
Predict_RNA_structure: $(PredictRNAStructure)
	$(CXX) $(CXXFLAGS) $(PredictRNAStructure) -o $(BIN)/$@

CalculateCGMDRMSD = Calculate_CGMD_rmsd.o combined.o RNA_molecule.o RNA_basepair.o RNA_coarse_grained.o read_lammps.o
Calculate_CGMD_rmsd: $(CalculateCGMDRMSD)
	$(CXX) $(CXXFLAGS) $(CalculateCGMDRMSD) $(LDFLAGS) -o $(BIN)/$@

CalculateHeavyatomRMSD = Calculate_heavyatom_rmsd.o combined.o
Calculate_heavyatom_rmsd: $(CalculateHeavyatomRMSD)
	$(CXX) $(CXXFLAGS) $(CalculateHeavyatomRMSD) $(LDFLAGS) -o $(BIN)/$@
	
CalculateTMscore = Calculate_TMscore.o combined.o
Calculate_TMscore: $(CalculateTMscore)
	$(CXX) $(CXXFLAGS) $(CalculateTMscore) $(LDFLAGS) -o $(BIN)/$@
	
CalculateRNAINF = Calculate_RNA_INF.o combined.o measure.o RNA_molecule.o RNA_basepair.o
Calculate_RNA_INF: $(CalculateRNAINF)
	$(CXX) $(CXXFLAGS) $(CalculateRNAINF) $(LDFLAGS) -o $(BIN)/$@

ExtractCGMDSnapshot = Extract_CGMD_snapshot.o combined.o read_lammps.o
Extract_CGMD_snapshot: $(ExtractCGMDSnapshot)
	$(CXX) $(CXXFLAGS) $(ExtractCGMDSnapshot) $(LDFLAGS) -o $(BIN)/$@
	
BuildRNAHelix = Build_RNA_helix.o combined.o RNA_molecule.o RNA_basepair.o
Build_RNA_helix: $(BuildRNAHelix)
	$(CXX) $(CXXFLAGS) $(BuildRNAHelix) $(LDFLAGS) -o $(BIN)/$@
	
ResetLammpsData = Reset_lammps_data.o combined.o read_lammps.o read_amber.o read_coeff.o prepare_lammps.o
Reset_lammps_data: $(ResetLammpsData)
	$(CXX) $(CXXFLAGS) $(ResetLammpsData) $(LDFLAGS) -o $(BIN)/$@
	
BuildRNAMotifStructure = Build_RNA_motif_structure.o RNA.o read_restraint.o parse_RNA_motif.o match_RNA_motif.o \
	find_RNA_template.o assemble_RNA_motif.o build_RNA_model.o read_amber.o read_lammps.o read_coeff.o prepare_lammps.o 
Build_RNA_motif_structure: $(BuildRNAMotifStructure)
	$(CXX) $(CXXFLAGS) $(BuildRNAMotifStructure) $(LDFLAGS) -o $(BIN)/$@

clear:
	rm -f *.o
	
clean:
	rm -f *.o $(BIN)/*
	
.PHONY: all clean clear
