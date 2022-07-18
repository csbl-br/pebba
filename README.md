<p align="center">
  <img  src="./PEBBA_banner.png">
</p>

---------------------------------------

Over-representation analysis (ORA) is a critical technique to determine if a set of differentially expressed genes (DEGs) is enriched with genes from specific gene sets or pathways. However, the cut-off used to define the number of DEGs that are utilised significantly impacts ORA results. To overcome the arbitrary choice of a cut-off and identify cut-off-independent enriched pathways, we developed PEBBA. This user-friendly tool ranks genes based on their statistical and biological significance and then systematically performs ORA for different cut-offs. There is no need to shortlist genes or waste time fine-tuning parameters. By simplifying ORA, PEBBA can be employed to lighten users’ burdens concerning parameter choice and decrease false positives. By visually exploring the parameter space, users can draw more precise conclusions about their dataset.

## Install
To install pebba, first clone this repository and install its dependencies with `pip install -r requirements.txt`.
Then, install pebba in editable mode (unless you intend to deploy it in production, **do not** deploy it in editable mode): `pip install -e .`

(soon pebba will be uploaded to pypi, thus simplifying this process)

## Using pebba
Once installed, pebba can be used as a standalone module:

`python -m pebba <deg_file> <gmt_file>`


For more options use the --help flag: `python -m pebba --help`

You can also use the test data to run a pebba analysis: 

`python -m pebba tests/data/GSE49757_Septic_vs_Healthy.txt tests/data/Reactome_2016_15and100Genes.gmt `

PEBBA is also available as an online tool at [pebba.sysbio.tools](https://pebba.sysbio.tools/)
