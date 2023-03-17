import logging

from retro_star.model import ValueMLP
from rxnfp.transformer_fingerprints import (RXNBERTFingerprintGenerator, get_default_model_and_tokenizer, generate_fingerprints)
from retro_star.utils import setup_logger
from retro_star.common import prepare_starting_molecules, \
    prepare_molstar_planner, smiles_to_fp, onmt_trans
from sklearn.metrics.pairwise import euclidean_distances
from retro_star.utils.load_bio_rxn import prepare_rxn_fps
import numpy as np
import torch


class RSPlanner(object):
    def __init__(self, gpu, expansion_topk, iterations, use_value_fn, do_rerank, buliding_block_path, fp_dim,
                 one_step_model_path, value_fn_model_path, rxn_fps_path, beam_size, route_topk, viz, viz_dir):

        setup_logger()
        device = torch.device('cuda:%d' % gpu if gpu >= 0 and torch.cuda.is_available() else 'cpu')
        print(device)
        starting_mols = prepare_starting_molecules(buliding_block_path)
        print("number of starting mols: ", len(starting_mols))
        beam_size = beam_size if beam_size > expansion_topk else expansion_topk
        one_step_handler = lambda x: onmt_trans(
            x,
            topk=expansion_topk,
            model_path=one_step_model_path,
            beam_size=beam_size,
            device=gpu
        )

        self.top_k = route_topk

        if do_rerank:
            model, tokenizer = get_default_model_and_tokenizer()
            rxnfp_generator = RXNBERTFingerprintGenerator(model, tokenizer)
            rxn_templates = prepare_rxn_fps(rxn_fps_path)
            rxn_fps = [pair[1] for pair in rxn_templates]
            rxn_ids = [pair[0] for pair in rxn_templates]
            
            def rerank_fn(subs, prod, scores):
                result = {}
                pred_rxns = [reactant+">>"+prod for reactant in subs]
                pred_rxns_fp = rxnfp_generator.convert_batch(pred_rxns)
                dis = euclidean_distances(pred_rxns_fp, rxn_fps)
                dis_norm = dis.T / dis.T.max(axis=0)
                id_list = np.argmin(dis_norm, axis=0)
                #scores = (np.log(dis_norm.min(0)) + np.array(scores))/2
                scores = (np.exp(dis_norm.min(0)) + scores)/2
                #scores = 1-dis_norm.min(0)
                ordered_tuples = sorted(enumerate(scores), key=lambda x:x[1])
                #length = min(len(subs), expansion_topk)
                ordered_reactants = [subs[each[0]] for each in ordered_tuples]
                ordered_scores = [scores[each[0]] for each in ordered_tuples]
                result['reactants'] = ordered_reactants
                result['scores'] = ordered_scores
                result['templates'] = [rxn_ids[idx] for idx in id_list]
                return result
        else:
            def rerank_fn(subs, prod, scores):
                result = {}
                length = min(len(subs), expansion_topk)
                result['reactants'] = subs[:length]
                result['scores'] = scores[:length]
                result['templates'] = ["" for i in range(0,length)]
                return result
        if use_value_fn:
            model = ValueMLP(
                n_layers=1,
                fp_dim=fp_dim,
                latent_dim=128,
                dropout_rate=0.1,
                device=device
            ).to(device)
            logging.info('Loading value nn from %s' % value_fn_model_path)
            model.load_state_dict(torch.load(value_fn_model_path, map_location=device))
            model.eval()

            def value_fn(mol):
                fp = smiles_to_fp(mol, fp_dim=fp_dim).reshape(1, -1)
                fp = torch.FloatTensor(fp).to(device)
                v = model(fp).item()
                return v
            
        else:
            value_fn = lambda x: 0.
        self.plan_handle = prepare_molstar_planner(
            expansion_handler=one_step_handler,
            value_fn=value_fn,
            rerank_fn=rerank_fn,
            starting_mols=starting_mols,
            iterations=iterations,
            viz=viz,
            viz_dir=viz_dir,
            route_topk=route_topk
        )
    def plan(self, target_mol):
        succ, msg = self.plan_handle(target_mol)

        #if succ:
        ori_list = msg[3]
        routes_list = []
        for i in ori_list:
            routes_list.append(i.serialize_with_score())
        if succ:
            return succ, routes_list[:self.top_k]
        else:
            return succ, routes_list[1:self.top_k +1]
        #return succ, routes_list
            
        #else:
        #    logging.info('Synthesis path for %s not found. Please try increasing '
        #                 'the number of iterations.' % target_mol)
        #    return None


if __name__ == '__main__':
    planner = RSPlanner(
        gpu=0,
        use_value_fn=True,
        iterations=100,
        expansion_topk=50
    )

    result = planner.plan('CCCC[C@@H](C(=O)N1CCC[C@H]1C(=O)O)[C@@H](F)C(=O)OC')
    print(result)

    result = planner.plan('CCOC(=O)c1nc(N2CC[C@H](NC(=O)c3nc(C(F)(F)F)c(CC)[nH]3)[C@H](OC)C2)sc1C')
    print(result)

    result = planner.plan('CC(C)c1ccc(-n2nc(O)c3c(=O)c4ccc(Cl)cc4[nH]c3c2=O)cc1')
    print(result)

