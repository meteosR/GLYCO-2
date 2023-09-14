## GLYCO-2 (GLYcan COverage, version 2) <br />
(A tool to calcaulte glycan coverage of glysoylated proteins, version 2)

### How to run GLYCO-2

1. Download GLYCO-2 package 

2. Install GLYCO-2
       
       conda create -n GLYCO python=3.8
       conda activate GLYCO
       pip install -r requirements.txt
       python3 -m pip install Cython && python3 setup.py build_ext --inplace
   
4. Run GLYCO-2.0 <br /><br />
        1) Prepare glycosylated protein pdb file(s) <br /> <br />
              - Make sure your pdb file does not contain water molecules. This only increases the filesize which makes process slow. <br />
              - Make sure your pdb follows the standard format of PDB.<br />
                            - GLYCO only recognizes atoms starting with "ATOM" or "HETA".<br />
                            - Each chain must have a different chain ID <br />

    * ##### How can I glycosylate proteins? (MATEO, PLEASE ADD REFERENCE PROPERLY)
    You can build glycans on your protein with other servers/programs such as CHARMM-GUI Glycan Reader & Modeler (1,2,3,4) or Glycosylator (5) etc. GLYCO 2 does not build glycans.<br />


        2) Select parameters<br />
             
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

### Example commands (MATEO, PLEASE CORRECT BELOW EXAMPLES, THERE IS NO GLYCO_EXAMPLES..ETC ..  NO IDEA WHAT IT MEANS)

Run glyco on the "glyco_examples/2_ha_man5_frame_10__BGL_BMA_AMA__23.pdb" file, with a distance cutoff of 23 and module type all_atom using 12 CPUs. The results will be saved on the 2_ha_man5_frame_10__BGL_BMA_AMA__23 folder.
Look for BGL,BMA,AMA glycans.

python3 glyco.py -pdb glyco_examples/2_ha_man5_frame_10__BGL_BMA_AMA__23.pdb -cutoff 23 -glycans BGL,BMA,AMA -out_folder 2_ha_man5_frame_10__BGL_BMA_AMA__23 -ncpu 12 -module all_atom

Run glyco on multiple files in an input folder /example/GLYCO-2/struct_input and average results. 
python3 glyco.py -in_folder /example/GLYCO-2/struct_input -cutoff 23 -glycans BGL,BMA,AMA -out_folder debug -ncpu 32 -module all_atom -average 
