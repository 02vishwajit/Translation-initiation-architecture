from pathlib import Path
import csv
from collections import Counter

path = Path(
    "14_expanded_study/orthology/six_genome_matrix/"
    "strict_six_way_orthologs_annotated_v2.tsv"
)
outdir = path.parent

strains = [
    "ZM4", "NOVO", "RSP", "ZM_NCIMB11163",
    "TSA_JW_SL_YS485", "CT_ATCC27405"
]

def classify(products):
    text = " ".join(products).lower()

    categories = [
        ("Central carbon metabolism",
         ["glycolysis", "glyceraldehyde-3-phosphate dehydrogenase",
          "pyruvate", "phosphofructokinase", "enolase",
          "phosphoglycerate", "aldolase", "triosephosphate",
          "glucose", "gluconate", "6-phosphogluconate"]),
        ("Pentose phosphate / sugar metabolism",
         ["pentose phosphate", "transketolase", "transaldolase",
          "ribose", "xylose", "arabinose", "sugar transporter"]),
        ("Pyruvate / acetyl-CoA metabolism",
         ["acetyl-coa", "pyruvate dehydrogenase",
          "pyruvate decarboxylase", "acetaldehyde dehydrogenase"]),
        ("Translation / ribosome",
         ["ribosomal protein", "translation initiation factor",
          "translation elongation factor", "aminoacyl-trna",
          "trna synthetase", "ribosome"]),
        ("Energy conversion / ATP synthesis",
         ["atp synthase", "atpase", "electron transfer",
          "ferredoxin", "cytochrome", "oxidoreductase"]),
        ("DNA replication / repair",
         ["dna helicase", "dna polymerase", "uvr", "recombinase",
          "holliday junction", "dna repair", "gyrase"]),
        ("Protein quality control",
         ["protease", "endopeptidase", "chaperone", "heat shock"]),
        ("Other / review annotation",
         [])
    ]

    for category, keywords in categories[:-1]:
        if any(keyword in text for keyword in keywords):
            return category
    return "Other / review annotation"

with path.open(newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

results = []
for row in rows:
    products = [row.get(f"{s}_product", "") for s in strains]
    result = {
        "ortholog_group": row["ortholog_group"],
        "functional_category": classify(products),
        "ZM4_product": row.get("ZM4_product", ""),
        "NOVO_product": row.get("NOVO_product", ""),
        "RSP_product": row.get("RSP_product", ""),
        "ZM_NCIMB11163_product": row.get("ZM_NCIMB11163_product", ""),
        "TSA_JW_SL_YS485_product": row.get("TSA_JW_SL_YS485_product", ""),
        "CT_ATCC27405_product": row.get("CT_ATCC27405_product", ""),
    }
    results.append(result)

out = outdir / "six_way_functional_classification.tsv"
with out.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(results[0]), delimiter="\t")
    writer.writeheader()
    writer.writerows(results)

counts = Counter(r["functional_category"] for r in results)

print(f"Saved: {out}")
print(f"Groups classified: {len(results)}")
print("\nPreliminary category counts:")
for category, count in counts.most_common():
    print(f"{category}\t{count}")

print("\nPotential carbon-metabolism candidates:")
for r in results:
    if r["functional_category"] in (
        "Central carbon metabolism",
        "Pentose phosphate / sugar metabolism",
        "Pyruvate / acetyl-CoA metabolism"
    ):
        print(
            f'{r["ortholog_group"]}\t'
            f'{r["functional_category"]}\t{r["ZM4_product"]}'
        )
