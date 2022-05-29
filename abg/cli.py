import click
import math
from pathlib import Path

import pandas as pd
from tabulate import tabulate

from abg.compute import ABGComputer, ABGResults
from abg import logger


def get_formatted_results(r: ABGResults):
    return [
        round(r.a * 180.0 / math.pi, 1),
        round(r.b * 180.0 / math.pi, 1),
        round(r.g * 180.0 / math.pi, 1),
        round(r.rmsd1, 2),
        round(r.rmsd2, 2),
        round(r.x, 2),
        round(r.y, 2),
        round(r.z, 2),
    ]


@click.command(
    help="a program to calculate the euler angles of alpha beta gamma of "
    "rna junctions. Can supply either a pdb or a csv file with a column 'pdb' "
    "to supply multiple pdbs"
)
@click.argument("input")
@click.option(
    "-h1", "--helix-1-res", help="residue nums in helix 1 seperated by ,", required=True
)
@click.option(
    "-h2", "--helix-2-res", help="residue nums in helix 1 seperated by ,", required=True
)
def main(input, helix_1_res, helix_2_res):
    logger.setup_applevel_logger()
    log = logger.get_logger("main")
    target_res_1, target_res_2 = [], []
    try:
        target_res_1 = [int(x) for x in helix_1_res.split(",")]
        target_res_2 = [int(x) for x in helix_2_res.split(",")]
    except:
        log.error("residues list was not supplied correctly!")
        exit()
    if len(target_res_1) < 2 or len(target_res_2) < 2:
        log.error("must supply at least 1 basepair per helix as a list of two residues")
        exit()
    log.info(f"residues in helix 1: {target_res_1}")
    log.info(f"residues in helix 2: {target_res_2}")
    extension = Path(input).suffix
    if extension == ".pdb":
        log.info(f"{input} is determined to be in PDB format")
        abg_computer = ABGComputer()
        r = abg_computer.compute(input, target_res_1, target_res_2)
        data = [get_formatted_results(r)]
        df = pd.DataFrame(data, columns="a,b,g,rmsd1,rmsd2,x,y,z".split(","))
        print(tabulate(df, headers="keys", tablefmt="simple", showindex=False))
    elif extension == ".csv":
        log.info(f"{input} is determined to be in csv format")
        log.info("will compute abg for all pdbs in 'pdb' column")
        df = pd.read_csv(input)
        if "pdb" not in df.columns:
            log.error("must have column 'pdb' in supplied csv")
            exit()
        log.info(f"{len(df)} rows detected!")
        df = df.join(
            pd.DataFrame(
                {
                    "a": [[] for _ in range(len(df))],
                    "b": [[] for _ in range(len(df))],
                    "g": [[] for _ in range(len(df))],
                    "rmsd1": [[] for _ in range(len(df))],
                    "rmsd2": [[] for _ in range(len(df))],
                    "x": [[] for _ in range(len(df))],
                    "y": [[] for _ in range(len(df))],
                    "z": [[] for _ in range(len(df))],
                },
                index=df.index,
            )
        )
        for i, row in df.iterrows():
            abg_computer = ABGComputer()
            r = abg_computer.compute(row['pdb'], target_res_1, target_res_2)
            data = get_formatted_results(r)
            df.at[i, ["a", "b", "g", "rmsd1", "rmsd2", "x", "y", "z"]] = data
        log.info("outputing data to abg_output.csv")
        df.to_csv("abg_output.csv", index=False)

    else:
        exit()


if __name__ == "__main__":
    main()
