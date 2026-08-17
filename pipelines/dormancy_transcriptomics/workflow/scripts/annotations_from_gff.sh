#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <input.gff> <output.tsv>"
    exit 1
fi

input="$1"
output="$2"

# Extract geneID/feature/product from GFF attribute column, URL-decoding the
# percent-encoded characters GFF uses for reserved symbols (commas, semicolons, etc.)
# in the product description.
awk -F'\t' 'BEGIN{OFS="\t"; print "geneID","feature","product"}
  $3 ~ /^(snoRNA|snRNA|rRNA|pseudogene|tRNA|transcript|lnc_RNA|mRNA)$/ {
    match($9, /Parent=([^;]+)/, g)
    match($9, /product=([^;]+)/, p)
    gsub(/%2C/, ",",  p[1])
    gsub(/%3B/, ";",  p[1])
    gsub(/%25/, "%",  p[1])
    gsub(/%20/, " ",  p[1])
    gsub(/%3D/, "=",  p[1])
    gsub(/%26/, "\&", p[1])
    gsub(/%2B/, "+",  p[1])
    gsub(/, transcript variant X[0-9]+$/, "", p[1])
    gsub(/, variant [0-9]+$/, "", p[1])
    print g[1], $3, p[1]
  }' "$input" > "$output"
