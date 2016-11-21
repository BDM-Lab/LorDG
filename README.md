
----------

## 3D Genome Structure Modeling by Lorentzian Objective Function

----------

#### Bioinformatics, Data Mining, Machine Learning (BDM) Laboratory, 
#### University of Missouri, Columbia MO 65211

----------

#### Developer: <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Tuan Trieu <br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Department of Computer Science <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; University of Missouri, Columbia <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Email: tuantrieu@mail.missouri.edu <br/>

#### Contact: <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Jianlin Cheng, PhD <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Department of Computer Science <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; University of Missouri, Columbia <br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Email: chengji@missouri.edu <br/>





## 1. Content of folders:
- bin: contains executable files
- example: contains example data and parameter files used to reconstruct chromosome/genome structures
- src: source code of LorDG in java




## 2. Usage ##

To run the tool, type: `java -jar 3DDistanceBaseLorentz.jar parameters.txt`


- Parameters are configured in the 'parameters.txt' file:
	+ **NUM**: number of models to generate
	+ **OUTPUT_FOLDER**: output folder
	+ **INPUT_FILE**: hi-C contact file, each line contains 3 numbers (separated by a space) of a contact, position_1 position_2 interaction_frequencies	
	+ **CONVERT_FACTOR**: the factor used to convert IF to distance, distance = 1/(IF^factor), when not specified, the program will search for it in range [0.1, 3.0], step = 0.1
	+ **CHROMOSOME_LENGTH**: remove it if there is only one chromosome. If there are multiple chromosomes in the input data, specify number of points (or beads) of chromosomes in the input data, separated by a comma. These numbers must be consistent with the input data.	
	+ **VERBOSE**: true or false to output gradient values during optmization
	+ **LEARNING_RATE**: learning rate for the optimization, if optimization fails, try to reduce this value
	+ **MAX_ITERATION**: maximum number of iterations, the optimization may converge before this number	

See in `/examples/` for sample files


## 3. Output ##

**LorDG** produces 4 files:
	
- **\*.pdb**: contains the model and can be visualized by pyMol, Chimera
- **\*\_log\_a\_number.txt**: contains the settings used to build the model and Spearman's correlation of reconstructed distances and input IFs
- **\*\_log.txt**: if NUM > 1, the files contains settings and average correlation of Spearman's correlations of separate models
- **\*\_coordinate_mapping.txt**: contains the mapping of genomic positions to indices in the model. Indices start from 0, while in pyMol or Chimera, id starts from 1


## 4. Disclaimer ##

The executable software and the source code of LorDG is distributed free of 
charge as it is to any non-commercial users. The authors hold no liabilities to 
the performance of the program.

## 5. Citations
**T. Trieu, J. Cheng. 3D Genome Structure Modeling by Lorentzian Objective Function. Nucleic Acids Research**
