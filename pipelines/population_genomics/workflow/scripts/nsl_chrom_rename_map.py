# Maps each chromosome name to its 1-based position in the chromosome list.
with open(snakemake.input.chromosomes) as f, open(snakemake.output.map, "w") as out, \
        open(snakemake.log[0], "w") as logf:
    n = 0
    for i, name in enumerate((l.strip() for l in f if l.strip()), start=1):
        out.write(f"{name}\t{i}\n")
        n = i
    logf.write(f"Wrote {n} chromosome renaming entries to {snakemake.output.map}\n")