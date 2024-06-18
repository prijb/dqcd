d = {
    "schema_version": 2,
    "corrections": [
        {
            "name": "NUM_Data_DEN_MC_sublead_pt_dimuon_eta_2022_syst",
            "description": "NUM_Data_DEN_MC_subleading_pt_dimuon_eta_2022_syst",
            "version": 1,
            "inputs": [
                {
                    "name": "sublead_pT",
                    "type": "real",
                    "description": "Subleading muon pT"
                },
                {
                    "name": "dimuon_eta",
                    "type": "real",
                    "description": "Dimuon eta"
                },
                {
                    "name": "ValType",
                    "type": "string",
                    "description": "sf or syst (currently 'sf' is nominal, and 'systup' and 'systdown' are up/down variations with total stat+syst uncertainties"
                }
            ],
            "output": {
                "name": "weight",
                "type": "real",
                "description": "Output scale factor (nominal) or uncertainty"
            },
            "data": {
                "nodetype": "binning",
                "input": "sublead_pT",
                "edges": [
                    0.,
                    3.0,
                    5.0,
                    10.0,
                    15.0,
                    30.0,
                    50.0,
                    float("inf")
                ],
                "content": [],
                "flow": "error"
            }
        }
    ]
}

"""
d_eta = {
    "nodetype": "binning",
    "input": "dimuon_eta",
    "edges": [
        -2.5,
        -2.0,
        -1.5,
        -1.0,
        -0.5,
        0.0,
        0.5,
        1.0,
        1.5,
        2.0,
        2.5,
    ],
    "content": [],
    "flow": "error"
}
"""
d_eta = {
    "nodetype": "binning",
    "input": "dimuon_eta",
    "edges": [
        -float("inf"),
        -2.5,
        -2.0,
        -1.5,
        -1.0,
        -0.5,
        0.0,
        0.5,
        1.0,
        1.5,
        2.0,
        2.5,
        float("inf")
    ],
    "content": [],
    "flow": "error"
}

d_syst = {
    "nodetype": "category",
    "input": "ValType",
    "content": [
        {
            "key": "sf",
            "value": 0
        },
        {
            "key": "systup",
            "value": 0
        },
        {
            "key": "systdown",
            "value": 0
        }
    ]
}

from copy import deepcopy as copy
from analysis_tools.utils import import_root
ROOT = import_root()
tf = ROOT.TFile.Open("/vols/cms/pb4918/scouting_sf_snt/kinematic_reweight_2022.root")
histo = tf.Get("kinematic_reweight_2022")

"""
for ib in range(len(d["corrections"][0]["data"]["edges"]) - 1):
    new_d_eta = copy(d_eta)
    for ibeta in range(len(new_d_eta["edges"]) - 1):
        new_d_syst = copy(d_syst)
        print(ib + 1, ibeta + 1, histo.GetBinContent(ib + 1, ibeta + 1))
        if ib == 6:
            print("Filling pT bin 50-inf with 1 +- 0")
            content = 1.
            error = 0.
            new_d_syst["content"][0]["value"] = content
            new_d_syst["content"][1]["value"] = content + error
            new_d_syst["content"][2]["value"] = content - error
        else:
            content = histo.GetBinContent(ib + 1, ibeta + 1)
            error = histo.GetBinError(ib + 1, ibeta + 1)
            new_d_syst["content"][0]["value"] = content
            new_d_syst["content"][1]["value"] = content + error
            new_d_syst["content"][2]["value"] = content - error
        new_d_eta["content"].append(new_d_syst)
    d["corrections"][0]["data"]["content"].append(new_d_eta)
"""

for ib in range(len(d["corrections"][0]["data"]["edges"]) - 1):
    new_d_eta = copy(d_eta)
    ibeta = 0
    for ibeta_iter in range(len(new_d_eta["edges"]) - 1):
        new_d_syst = copy(d_syst)
        if ib == 6:
            print("Filling pT bin 50-inf with 1 +- 0 for ({ib}, {ibeta})".format(ib=ib, ibeta=ibeta_iter))
            content = 1.
            error = 0.
            new_d_syst["content"][0]["value"] = content
            new_d_syst["content"][1]["value"] = content + error
            new_d_syst["content"][2]["value"] = content - error
        else:
            if (ibeta_iter==0) | (ibeta_iter==11):
                print("Filling eta bin with inf edge with 1 +- 0 for ({ib}, {ibeta})".format(ib=ib, ibeta=ibeta_iter))
                content = 1.
                error = 0.
                new_d_syst["content"][0]["value"] = content
                new_d_syst["content"][1]["value"] = content + error
                new_d_syst["content"][2]["value"] = content - error
            else:
                #Only case where a non-1 value is needed
                print(ib + 1, ibeta + 1, histo.GetBinContent(ib + 1, ibeta + 1))
                content = histo.GetBinContent(ib + 1, ibeta + 1)
                error = histo.GetBinError(ib + 1, ibeta + 1)
                new_d_syst["content"][0]["value"] = content
                new_d_syst["content"][1]["value"] = content + error
                new_d_syst["content"][2]["value"] = content - error
                ibeta += 1

        new_d_eta["content"].append(new_d_syst)
    d["corrections"][0]["data"]["content"].append(new_d_eta)

import json
with open("data/NUM_Data_DEN_MC_subleading_pt_dimuon_eta_2022_syst.json", "w+") as f:
    json.dump(d, f, indent=4)
    