import requests

SERVER = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json"}

def get_mane_select(symbol, species="homo_sapiens"):
    # Step 1: find gene
    url = f"{SERVER}/lookup/symbol/{species}/{symbol}?expand=1"
    r = requests.get(url, headers=HEADERS)

    if not r.ok:
        return None

    data = r.json()

    # Step 2: scan transcripts
    for tx in data.get("Transcript", []):
        if "MANE_Select" in tx.get("tag", []):
            return {
                "gene": symbol,
                "ensembl_transcript": tx["id"],
                "mane_refseq": tx.get("RefSeq_match"),
                "biotype": tx.get("biotype")
            }

    return None

genes = ["AKT1", "AKT2", "TP53", "EGFR"]

results = []
for g in genes:
    print(g)
    mane = get_mane_select(g)
    results.append(mane)

for r in results:
    print(r)