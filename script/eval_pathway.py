import os
from rdkit import Chem

block_list = []
for line in open("./dataset/building_block/building_block.csv"):
    smi = line.strip().replace("[CoA]", "[33S]")
    smi = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
    block_list.append(smi.replace("[33S]", "[CoA]"))  

gt_route = {}
n = 0
for line in open("./dataset/test_data/biosynthetic_set.txt"):
    data = line.strip().split("\t")
    molid = data[0]
    route = []
    
    blocks = []
    branches = 0
    mols = []
    rxns = data[1].split("|")

    for rxn in rxns:
        if "." in rxn.split(">>")[1]:
            branches += 1
            mols.extend(rxn.split(">>")[1].split("."))
        else:
            mols.append(rxn.split(">>")[1])
        mols.append(rxn.split(">>")[0])
        #rxn = rxn.replace("[CoA]", "*").replace("[SAH]", "")
        src = rxn.split(">>")[0]
        tar = rxn.split(">>")[1]
        
        for item in tar.split("."):
            if (item in block_list) or (item.upper().count("C") < 4):
                blocks.append(item)
        #tar = Chem.MolToSmiles(Chem.MolFromSmiles(tar))
        route.append(tar)
    if n in gt_route.keys():
        gt_route[n].append((route, blocks))
    else:
        gt_route[n] = [(route, blocks)]
    n += 1
with open("result/result.txt", "w") as f:
    f.write("mol_id\tsuccess\tpathway_hit\trank\tblock_hit\tintersect\tlength\tsolution\n")
for mol_id in gt_route.keys():
    #print(gt_route[mol_id][0][0])
    #print(gt_route[mol_id][0][1])
    match = False
    block_match = False
    max_intersect = 0
    max_pred = 0
    solution = 0
    succeed = False
    match_list = []
    max_id = 1
    k = 0
    rank = k
    if os.path.exists("result/mol_"+str(mol_id)+".txt"):
        succeed = True
        m = 0
        for line1 in open("result/mol_"+str(mol_id)+".txt"):
            m += 1
            if m >5:
                break
            k += 1
            if ">" not in line1:
                continue
            solution += 1
            pred_route = []
            data1 = line1.strip().split("\t")
            rxns = data1[0].split("|")
            pred_blocks = []
            for rxn in rxns:
                rxn = rxn.replace("[CoA]", "[33S]").replace("[SAH]", "[32S]")
                #rxn = rxn.replace("*", "[33S]").replace("R", "*")  ### for old version
                src = rxn.split(">")[0]
                tar = rxn.split(">")[2]
                for item in tar.split("."):
                    if (item in block_list) or (item.upper().count("C") < 4):
                        pred_blocks.append(item.replace("[33S]", "[CoA]").replace("[32S]", "[SAH]"))
                        #pred_blocks.append(item.replace("[33S]", "[CoA]"))   ### for old version
                #src = Chem.MolToSmiles(Chem.MolFromSmiles(src))
                tar = Chem.MolToSmiles(Chem.MolFromSmiles(tar))
                pred_route.append(tar.replace("[33S]", "[CoA]").replace("[32S]", "[SAH]"))
                #pred_route.append(tar.replace("[33S]", "[CoA]"))   ### for old version
            i = 1
            for pair in gt_route[mol_id]:
                #print(set(pred_blocks))
                #print(set(pair[1]))
                #print(set(pair[1]) == set(pred_blocks))
                if set(pair[1]) == set(pred_blocks):
                    block_match = True
                intersect = len(set(pair[0]) & set(pred_route))
                if set(pair[0]) == set(pred_route):
                    match = True
                    match_list.append(i)
                    rank = k
                if intersect > max_intersect:
                    max_intersect = intersect
                    max_id = i
                if len(rxns) > max_pred:
                    max_pred = len(rxns)
                i += 1
    with open("result/result.txt", "a") as f:
        f.write("mol_"+str(mol_id) + "\t" + str(succeed) + "\t" + str(match) + "\t" + str(rank) + "\t" + str(block_match) + "\t" + str(max_intersect) + "\t" + str(max_pred) + "\t" + str(solution) + "\n")
