<p align="center">
  <img  src="https://raw.githubusercontent.com/csbl-br/pebba/master/PEBBA_banner.png">
</p>

---------------------------------------

Over-representation analysis (ORA) is a critical technique to determine if a set of differentially expressed genes (DEGs) is enriched with genes from specific gene sets or pathways. However, the cut-off used to define the number of DEGs that are utilised significantly impacts ORA results. To overcome the arbitrary choice of a cut-off and identify cut-off-independent enriched pathways, we developed PEBBA. This user-friendly tool ranks genes based on their statistical and biological significance and then systematically performs ORA for different cut-offs. There is no need to shortlist genes or waste time fine-tuning parameters. By simplifying ORA, PEBBA can be employed to lighten users’ burdens concerning parameter choice and decrease false positives. By visually exploring the parameter space, users can draw more precise conclusions about their dataset.

## Install
PEBBA may be installed through pip: `pip install pebba`.
Alternatively, PEBBA is also available as an online tool at [pebba.sysbio.tools](https://pebba.sysbio.tools/)


## Using pebba
Once installed, pebba can be used as a standalone module:

`python -m pebba <deg_file> <gmt_file>`


For more options use the --help flag: `python -m pebba --help`
