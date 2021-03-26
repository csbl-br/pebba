# Pebba


TODO: description

## Install
To install pebba, first clone this repository and install its dependencies with `pip install -r requirements.txt`.
Then, install pebba in editable mode (unless you intend to deploy it in production, **do not** deploy it in editable mode): `pip install -e .`


## Using pebba
Once installed, pebba can be used as a standalone module:

`python -m pebba <deg_file> <gmt_file>`


For more options use the --help flag: `python -m pebba --help`

You can also use the test data to run a pebba analysis: 

`python -m pebba tests/data/GSE49757_Septic_vs_Healthy.txt tests/data/Reactome_2016_15and100Genes.gmt `
