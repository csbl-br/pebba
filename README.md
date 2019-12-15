# Python-PEBBA


Fischer exact test + FDR correction BH

A ideia é:
Dado um numero de genes diferencialmente expressos N (top_genes no meu caso) e um numero total M de genes no estudo, qual a probabilidade de que dada uma via de tamanho n eu tenha X genes diferencialmente expressos?

Eu vou percorrendo via por via e realizando o teste para ver se o numero de genes diferencialmente expressos nela é suficiente para dizer que ela esta sendo diferencialmente expressa.




implementacao do teste na biblioteca q eu to usando (pypathway):
https://github.com/iseekwonderful/PyPathway/blob/master/pypathway/analysis/ora/__init__.py#L58


teste hipergeometrico (hipergeom sf):

https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html
https://github.com/scipy/scipy/blob/v0.14.0/scipy/stats/_discrete_distns.py#L347
https://github.com/scipy/scipy/blob/master/scipy/stats/_discrete_distns.py#L478

fdr_bh:
https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
https://www.statsmodels.org/stable/_modules/statsmodels/stats/multitest.html#multipletests
