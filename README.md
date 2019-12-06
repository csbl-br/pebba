# Python-PEBBA


Fischer exact test + FDR correction BH



Teste exato de fischer=?=teste hipergeometrico

#2X2
#vc tem duas amostras(doente e saudavel) e compara se a proporcao de genes super expressos em cada uma é diferente de alguma #maneira relevante.




Lista de genes associados a uma via metabolica  e alguns genes dessa via q estao super-expressos. Usando essas informacoes se faz um teste para ver se pode se dizer q a via como um todo esta alterada

implementacao do teste na biblioteca q eu to usando (pypathway):
https://github.com/iseekwonderful/PyPathway/blob/master/pypathway/analysis/ora/__init__.py#L58


teste hipergeometrico (hipergeom sf):
https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.hypergeom.html
https://github.com/scipy/scipy/blob/v0.14.0/scipy/stats/_discrete_distns.py#L347

fdr_bh:
https://www.statsmodels.org/stable/generated/statsmodels.stats.multitest.multipletests.html
https://www.statsmodels.org/stable/_modules/statsmodels/stats/multitest.html#multipletests
