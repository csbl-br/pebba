import time
import pandas as pd
from pebba.main import pebba

# Since python doesn't allow imports over parent directory, change sys path or location of this file
def timing_pebba(n):
    deg_in = "data/GSE49757_Septic_vs_Healthy.txt"
    gmt_in = "data/Reactome_2016_15and100Genes.gmt"
    run_times = []
    for i in range(n):
        print("Iteration number " + str(i + 1) + " out of " + str(n))
        time_start = time.time()
        pebba(deg_in, gmt_in, force=True)
        run_times.append(time.time() - time_start)
    print(run_times)
    print("acabei")
    df = pd.DataFrame.from_dict(
        {str(i + 1): run_times[i] for i in range(len(run_times))},
        orient="index",
        columns=["run_time: "],
    )
    df.to_csv("timing_performance.csv")
    return


timing_pebba(10)
# for n = 10 we got a mean of 273.2
# near 4,5 minutes
