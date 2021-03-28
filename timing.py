import os
import sys
import time
import pandas as pd
from datetime import datetime
from pebba.main import pebba

print("starting")


def timing_pebba(n=10):
    run_times = run_n_times(n)
    save_results(run_times, n)

    print("Results Saved")
    return


def run_n_times(n):
    deg_in = "tests/data/GSE49757_Septic_vs_Healthy.txt"
    gmt_in = "tests/data/Reactome_2016_15and100Genes.gmt"
    run_times = []
    for i in range(n):
        print("Iteration number " + str(i + 1) + " out of " + str(n))
        time_start = time.time()
        pebba(deg_in, gmt_in, force=True)
        run_times.append(time.time() - time_start)
    return run_times


def save_results(run_times, n):
    df_new_run = get_df_for_new_run(run_times)
    try:
        df = pd.read_csv("timing_performances_n" + str(n) + ".csv", index_col=0)
        df.index = df.index.astype(int)
        df = pd.merge(df, df_new_run, left_index=True, right_index=True)
    except FileNotFoundError:
        df = df_new_run
    finally:
        df.to_csv("timing_performances_n" + str(n) + ".csv")


def get_df_for_new_run(run_times):
    time_stamp = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

    df_new_run = pd.DataFrame.from_dict(
        {str(i + 1): run_times[i] for i in range(len(run_times))},
        orient="index",
        columns=[time_stamp + " run_time: "],
    )
    df_new_run.index.rename("Run Number i=", inplace=True)
    df_new_run.index = df_new_run.index.astype(int)
    return df_new_run


if __name__ == "__main__":
    timing_pebba()

# In the future, add a clean up step to delete the results folder
