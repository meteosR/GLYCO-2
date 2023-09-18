## GLYCO-2 (version 2) <br />
(A tool to calcaulte glycan coverage of glysoylated proteins, version 2)

### How to run GLYCO-2

1. Download GLYCO-2 package 

2. Install GLYCO-2
       
       conda create -n GLYCO2 python=3.8
       conda activate GLYCO2
       pip install -r requirements.txt
       python3 -m pip install Cython && python3 setup.py build_ext --inplace
   
4. Run GLYCO-2 <br />
(How to run is pretty much the same as [GLYCO ver.1](https://github.com/myungjinlee/GLYCO/blob/main/README.md)) <br /> <br />
        1) Prepare glycosylated protein pdb file(s) <br /> <br />
              - Make sure your pdb file does not contain water molecules. This only increases the filesize which makes process slow. <br />
              - Make sure your pdb follows the standard format of PDB.<br />
                            - GLYCO only recognizes atoms starting with "ATOM" or "HETA".<br />
                            - Each chain must have a different chain ID <br />

    * ##### How can I glycosylate proteins? 
    You can build glycans on your protein with other servers/programs such as [CHARMM-GUI Glycan Reader & Modeler](https://charmm-gui.org/?doc=input/glycan) or [Glycosylator](https://github.com/tlemmin/glycosylator) etc. GLYCO-2 does not build glycans.<br /><br />
        2) Arguments to enter <br />
        <img width="925" alt="image" src="https://github.com/meteosR/GLYCO-2/assets/32939217/68da2c2c-90aa-4b8b-9b0e-88d2cc103cce"><br />
      <br />  Please make sure you are in GLYCO environment<br />
      
       conda activate GLYCO2

### Example commands 
       Single PDB
       python3 glyco.py -pdb FILENAME.pdb -cutoff 23 -glycans BGL,BMA,AMA -out_folder FILENAME -ncpu 12 -module all_atom

       Multiple PDBs
       python3 glyco.py -in_folder FOLDERNAME -cutoff 23 -glycans BGL,BMA,AMA -out_folder OUT_FOLDERNAME -ncpu 32 -module all_atom -average 

### How to analyze output
GLYCO-2 outputs two main output files and six accessory files:<br /><br />
Main output:<br />
  1) X_bfactor.pdb: input PDB with the glycan coverage for each atom in the b-factor column. Please load it in PyMOL and visualize with a command. Of course you can change the color as you wish.<br />
  
            spectrum b, white_green_black 
            
  2) X.csv: Ths file contains the quantified glycan coverage for each protein residue. The file contains the following columns:<br />
    - Protein_ID: The unique protein residue identifier used during the computation of glycan coverage.<br />
    - Glycans_atoms: The list of all the glycan atoms covering the given protein residue. The length of this list corresponds to the coverage<br />
    - Glycan_density: The glycan coverage of this protein residue<br />
    - Protein Chain: The chain of the residue (Example: A)<br />
    - Protein residue: The residue name (Example: ALA)<br />
    - Protein residue position: The position of the residue (Example: 123B)<br />

Accessory files:<br />

  3) X.csv:output file after running FreeSASA. This file contains information about which residues were determined to be in the surface.
However, does not directly relate to what you want<br />

  4) log.txt: A internal log with the computation steps. If you encounter errors, you can email us this log with questions for additional help. It also shows time spent for calculation at the end of the file.<br />
  
  5) glysums_X.txt: (MATEO, WHAT DOES GLYCANSUMS HAVE? <br />
  
  6) params_in.txt: The file contains input parameters that the user entered.<br />

