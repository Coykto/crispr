import re

from itertools import product, permutations, groupby



genomes = {
    "genomes": {
        "hg38":
            {
                "genome": "hg38",
                "species": "human",
                "chromosomes": "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y,M",
            },
        "m39":
            {
                "genome": "m39",
                "species": "mouse",
                "chromosomes": "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,X,Y,M",
            },
        "hg19":
            {
                "genome": "hg19",
                "species": "human",
                "chromosomes": "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,X,Y,M",
            },
    }
}

config = [
    {
        "field": "organism",
        "options": [
            {
                "name": "human",
                "genomes": [
                    {
                        "name": "hg38",
                        "chromosomes": [
                            "1", "2", "3", "4", "5", "6", "7", "8", "9", "10","11", "12", "13",
                            "14", "15", "16", "17", "18", "19", "20","21", "22", "X", "Y", "M"
                        ]
                    },
                    {
                        "name": "hg19",
                        "chromosomes": [
                            "1", "2", "3", "4", "5", "6", "7", "8", "9", "10","11", "12", "13",
                            "14", "15", "16", "17", "18", "19", "20","21", "22", "X", "Y", "M"
                        ]
                    }
                ],
            },
            {
                "name": "mouse",
                "genomes": [
                    {
                        "name": "m39",
                        "chromosomes": [
                            "1", "2", "3", "4", "5", "6", "7", "8", "9", "10","11", "12", "13",
                            "14", "15", "16", "17", "18", "19", "X", "Y", "M"
                        ]
                    }
                ],
            },
        ]
    }
]

if __name__ == '__main__':
    genome_options = []
    genome_list = sorted([v for k, v in genomes["genomes"].items()], key=lambda x: x["species"])
    for species, genome_items in groupby(genome_list, key=lambda x: x["species"]):
        genome_options.append(
            {
                "name": species,
                "genomes": [
                    {
                        "name": genome["genome"],
                        "chromosomes": genome["chromosomes"].split(",")
                    } for genome in genome_items
                ]
            }
        )
    result = {
        "field": "organism",
        "options": genome_options
    }


    awd = 23