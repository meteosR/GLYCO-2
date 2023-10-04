## GLYCO-2 (version 2) <br />
(A tool to calcaulte glycan coverage of glysoylated proteins, version 2)

### How to run GLYCO-2

1. Download GLYCO-2 zip file 

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
        2) FreeSASA(https://freesasa.github.io/) - You can use FreeSASA executable in Dependencies folder, or download it by yourself<br /><br />
        3) Arguments to enter <br /><br />
![diagram](https://github.com/meteosR/GLYCO-2/assets/32939217/19628413-a6da-4a49-bb59-917a2854cc2d)


      <br />
   *Argument in green box is mandatory and white box is optional.<br />
   *ncpu x npar = total number of CPUs allocated<br />
   *Use an argument -verbose to print out detailed log file, if you need<br /> <br />
   
   Please make sure you are in GLYCO2 environment<br />
      
       conda activate GLYCO2

### Example commands 
Single PDB<br />

       python3 glyco.py -pdb FILENAME.pdb -module all_atom -freesasa FREESASA_PATH -cutoff 23 -glycans BGL,BMA,AMA -ncpu NUM_CPU -out_folder FILENAME

Multiple PDBs with module "subset" to calculate predefined residues<br />

       python3 glyco.py -in_folder FOLDERNAME -average -module subset -residue RESLIST.txt -cutoff 23 -glycans BGL,BMA,AMA -ncpu NUM_CPU -out_folder OUT_FOLDERNAME -average


*RESLIST.txt should be in the following format: residue name, chain ID, residue number     <br />
    (residuelist.txt)    <br />
    LEU  A  95    <br />
    VAL  A  97    <br />
    THR  A  98<br />
    
### How to analyze output
GLYCO-2 outputs three main output files and two accessory files:<br /><br />
Main outputs in result folder:<br />

  1) X_bfactor.pdb: input PDB with the glycan coverage for each atom in the b-factor column. Please load it in PyMOL and visualize with a command. Of course you can change the color as you wish.<br />
  
            spectrum b, white_green_black 
            
  2) X.csv: Ths file contains the quantified glycan coverage for each protein residue. The file contains the following columns:<br />
    - Chain: ChainID<br />
    - ResNum: Residue number <br />
    - ResName: Residue name <br />
    - Glycan_density: The glycan coverage of the protein residue<br />
    - SASA ABS: Solvent accessible surface area from FreeSASA of the residue<br />

  3) glysums_X.txt: <br />
    - First line: Sum of glycan density<br />
    - Second line: Sum of glycan density excluding any overlaps<br />
  
Accessory files:<br />

  4) params_in.txt: The file contains input parameters that the user entered.<br />
  
  5) log.txt: A internal log with the computation steps. If you encounter errors, you can email us this log with questions for additional help. It also shows time spent for calculation at the end of the file.<br />


