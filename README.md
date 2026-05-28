# Melange Filter

### A filtering tool for [Melange](https://github.com/sandragodinhosilva/melange)

Melange Filter allows for the selection and filtering of specific annotations from previously calculated Melange results
It does not modify your Melange results

if you want to calculate PCOAs from Melange results check [PCOA_maker](https://github.com/JFrAlmeida/PCOA_maker.git)
It can be run in the conda environment generated with this tool with small additions (check Instalation)

### Provides:
- Single table output of selected features, with DNA and amino acid sequences for each feature
- Addition of external dbcan results and merging into the Melange results
- DNA or amino acid Fasta files of the selected features
- Stacked bar plots of user defined groups of annotations (e.g. endochitinases vs exochitinases vs LPMOs etc.)


### Files you should interact with
There are 2 files you should interact with to get full use of this:
- Config.py (Mandatory)
- tree order.csv (Optional, In MF statistics, explained in Config.py)

### Config.py
General controls over the entire tool

### Tree_order.csv
Optional metadata file.

Allows setting the order of the stacked bar plot, and setting the color for each label

Requires two columns columns:
    Tree_order: Controls sample order in the figure.

    Color: Controls label colors.


## Instalation

Navigate in your terminal window to the folder you want to do you analysis in.
Move the "Annotation results" and "Annotation" folders containing the Melange
results to that folder

git clone https://github.com/JFrAlmeida/Melange_filter.git

conda env create -f environment.yml
conda activate Melange_filter


### Bonus, combine with PCOA_maker

-> after activating the conda environment Melange_filter

Additionally install the following
```
conda install -c conda-forge scipy=1.10.1 scikit-bio=0.6.3
```
git clone https://github.com/JFrAlmeida/PCOA_maker.git

and run PCOA_maker (after checking instructions there) with
```
python PCOA_Maker.py
```

conda create -n melangextra
conda activate melangextra
conda install -c conda-forge python=3.13.9 
conda install -c conda-forge pandas




pip install biopython==1.86

---USAGE---

```
python Melange_filter.py
```


