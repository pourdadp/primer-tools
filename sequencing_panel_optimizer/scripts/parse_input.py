
import yaml
import json

with open("config/config.yaml", encoding='utf-8') as f:
    config = yaml.safe_load(f)

data = {
    "sequence": config["dna_sequence"],
    "primers": config["primers"],
    "max_mispairing": config["max_mispairing"],
    "max_product_length": config["max_product_length"],
    "min_product_length": config.get("min_product_length", 0),
    "max_tm_diff": config.get("max_tm_diff", 5.0)
}

with open("results/input_data.json", "w", encoding='utf-8') as out:
    json.dump(data, out, indent=2)
