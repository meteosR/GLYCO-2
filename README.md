# GLYCO-2

## GLYcan COverage (A webserver to calcaulte glycan coverage of glysoylated proteins)

### How to install

       conda create -n GLYCO python=3.8
       
       conda activate GLYCO <br />
       
       pip install -r requirements.txt <br />
       
       python3 -m pip install Cython && python3 setup.py build_ext --inplace <br />


### How to run GLYCO-2.0
* #### Step 1: Enter input pdb file(s)
    Prepare your glycosylated protein pdb file and upload it. {# #}
Make sure your pdb file does not contain water molecules if it comes from an MD simulation, for example. 
This only increases the size of the file. 
Also, make sure your pdb follows the standard format as defined here. GLYCO only recognizes lines starting with "ATOM" or "HETA".
Each chain must have a different chain ID (you should not have two Alanine 23 with the same chain, for example).

    * ##### How can I glycosylate proteins?
    You can build glycans on your protein with other servers/programs such as CHARMM-GUI Glycan Reader & Modeler (1,2,3,4) or Glycosylator (5) etc. GLYCO 2 does not build glycans.

* #### Step 2: Select parameters
  - Module type is either "all_atom" or "subset".
  - Module "all_atom" calculates glycan coverage for all surface protein residues, and module "sub" calculates glycan coverage only for user input residues.
  - If "subset" is selected, you should provide a list of protein residues of interest as an input.
  - If "All protein atoms" is selected, you should provide the Minimum RSA value in Å² to determine the surface residues.
  - Carefully check which glycan names are found on your pdb files (Example: BMA, AMA).
  - Distance cutoff parameter is generally determined by a histogram of the longest length of glycans (for example, distance between a carbon (C1) in a glycan bonded to nitrogen (ND2) in a protein) of the glycosylated protein. From histograms of Man-5 HIV-1, Influenza HA, Lassa, SARS-CoV-2 500 ns molecular dynamics simulations, the longest length of Man-5 glycan was determined to be 23 Å and 26 Å for Man-9. Users can freely alter the distance cutoff based on their structures.
  - Cylinder radius is set to 1.4 Å by default.
  - Since atomic radii of carbon, nitrogen, and oxygen are around 0.6 - 0.8 Å (GLYCO only considers non-hydrogen atom), we recommends cylinder radius to be twice the atomic radii Example: 1.2 or 1.6 Å.
  - Surface area cutoff SASA (solvent accessible surface area) is normally set to 30 Å² to only select protein surface residues, but you are free to alter it. Probe radius is set to 1.4 Å by default.

### How to analyze output
GLYCO 2.0 an output folder that includes four files:
  - X_bfactor.pdb: PDB containing the glycosylated structure and the glycan coverage for each atom in the b-factor column.
Download X_bfactor.pdb and open it through PyMOL (6). You can visualize glycan coverage on the protein surface with the command "spectrum b, white_black". Of course you can change the color as you wish.
  - X.rsa: output file after running FreeSASA. This file contains information about which residues were determined to be in the surface.
However, does not directly relate to what you want
  - log.txt: A internal log with the computation steps. If you encounter errors, you can email us this log with questions for additional help.
You can ignore this log in most cases.
  - X.csv: Ths file contains the quantified glycan coverage for each protein residue. The file contains the following columns:  
    - Protein_ID: The unique protein residue identifier used during the computation of glycan coverage.
    - Glycans_atoms: The list of all the glycan atoms covering the given protein residue. The length of this list corresponds to the coverage
    - Glycan_density: The glycan coverage of this protein residue
    - Protein Chain: The chain of the residue (Example: A)
    - Protein residue: The residue name (Example: ALA)
    - Protein residue position: The position of the residue (Example: 123B)

### Example commands

Run glyco on the "glyco_examples/2_ha_man5_frame_10__BGL_BMA_AMA__23.pdb" file, with a distance cutoff of 23 and module type all_atom using 12 CPUs. The results will be saved on the 2_ha_man5_frame_10__BGL_BMA_AMA__23 folder.
Look for BGL,BMA,AMA glycans.
python3 glyco.py -pdb glyco_examples/2_ha_man5_frame_10__BGL_BMA_AMA__23.pdb -cutoff 23 -glycans BGL,BMA,AMA -out_folder 2_ha_man5_frame_10__BGL_BMA_AMA__23 -ncpu 12 -module all_atom

Run glyco on multiple files in an input folder /example/GLYCO-2/struct_input and average results. 
python3 glyco.py -in_folder /example/GLYCO-2/struct_input -cutoff 23 -glycans BGL,BMA,AMA -out_folder debug -ncpu 32 -module all_atom -average 
