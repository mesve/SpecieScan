[README_PY]
![Python Version](https://img.shields.io/badge/Python-3.x-brightgreen.svg)

## Description

This file provides documentation for SpecieScan automated species identification algorithm in Python. 

## Features

- Imports pre-processed mass spectrometry data from MALDI-TOF MS (see README_R on preprocessing).
- Matches processed spectra peaks to a reference database.
- Ranks species based on correlation scores.
- Batch processes multiple files for efficiency.
- Identifies potential contaminants in unknown samples.

## Prerequisites

- Python 3.x
- Required Python packages (install using `pip install package_name`):
  - `pandas`
  - `numpy`
  - `scipy`
  - `matplotlib`
  - `glob`
  - `csv`
  - `os`

  ## Usage

1. Clone this repository to your local computer.

- https://github.com/mesve/SpecieScan

2. Navigate to the project directory.

3. Install the required Python packages if not already installed.
- The algorithm is in Jupyter notebook - see section ## Usage with Jupyter Notebook below on installation and use.

4. Make 3 folders on your computer and copy all preprocessed samples (a batch) in each of them.
    Naming: 
        folder 1: /Batch_name/_all_d 
        folder 2: /Batch_name/_corr6 # correlation scores for the first 6 most likely species
        folder 3: /Batch_name/_contaminants

5. Place your preprocessed CSV files of unknown samples in the specified directory.

6. Run the Python script to process and analyse the data. 

- Modify the script as needed for your specific data and analysis parameters based on the comments.

6.1 First the code will process one file (first chosen as df), then it will batch process all files from the chosen folder.

6.2 Once all all_d are batch processed, the next code snippet deletes the non-all_d files from the folder.

6.3 all_d files are merged

6.4 Same unknown file is read in as df (copy from above)

6.5 run transform_correlation to get the first n number of more likely species and the correlation scores

6.6 Batch process corr6 files (path needs to be populated with the preprocessed files, not the all_d files!)

6.7 Merge corr6 files 

6.8 Merge the merged all_d files and merged corr6 files.

7. Review the generated output files, including species identification and correlation scores.

8. Scan for contamination in all samples

8.1 Read in one preprocessed file and the contamination reference file included in the Zenodo documentation.

8.2 Batch processing: define path with the preprocessed files (folder 3, ending in _contaminants) and run the code. This will save all files with their contaminants.

8.3 Run code to delete the non-contaminants files in the contaminants folder.

8.4 Run code to concatenate all contaminants files.

## Usage with Jupyter Notebook

1. Install Jupyter Notebook:

- If you don't have Jupyter Notebook installed, you can install it using `pip`:

  ```
  pip install notebook
  ```

2. Start Jupyter Notebook:

- Open a terminal or command prompt and navigate to the project directory (`SpecieScan`).

- Run the following command:

  ```
  jupyter notebook
  ```

3. In your web browser, open the provided Jupyter Notebook URL (usually `http://localhost:8888`) to access the notebook interface.

4. Open the provided Jupyter Notebook file (`SpecieScan.ipynb`) to interactively run and modify the code. Follow the instructions within the notebook.


## Configuration

- Modify the script and adjust parameters as needed to suit your specific data and analysis requirements.

## Output

- Processed data files will be saved in the specified directories.
- Identified species and correlation scores will be displayed or saved in output files.

## Contributing

Contributions are welcome and greatly appreciated! If you have suggestions, improvements, or bug fixes, please open an issue or submit a pull request.

## License

This project is licensed under the MIT License. See LICENSE.txt for more information.
Cite paper: Végh & Douka 2023. SpecieScan: Automated Species Identification Algorithm for Bone Fragments from MALDI-ToF MS Spectra.
            DOI:
& cite code: DOI: 10.5281/zenodo.8055426


## Contact

For questions or inquiries, please contact [Dr Emese Végh](mailto:mesve3@gmail.com) -- Twitter: @emese_vegh
