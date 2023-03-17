import os

from tqdm import tqdm
from rdkit import Chem

from config import Config
from retro_star.api import RSPlanner



def run(conf):
    #os.environ['CUDA_VISIBLE_DEVICES'] = '0, 1, 2, 3'

    # canonicalization
    mol = Chem.MolToSmiles(Chem.MolFromSmarts(conf.target_mol))

    planner = RSPlanner(
        gpu=conf.gpu,
        use_value_fn=conf.use_value_fn,
        value_fn_model_path=conf.value_fn_model_path,
        fp_dim=conf.fp_dim,
        iterations=conf.iterations,
        expansion_topk=conf.expansion_topk,
        route_topk=conf.route_topk,
        buliding_block_path=conf.buliding_block_path,
        one_step_model_path=conf.one_step_model_path,
        beam_size=conf.beam_size,
        viz=conf.viz,
        viz_dir=conf.viz_dir,
        do_rerank = conf.do_rerank, 
        rxn_fps_path=conf.rxn_fps_path,
    )

    succ, result = planner.plan(mol)
    result_routes_list = []

    if result is None:
        return None

    for i, route in enumerate(result):
        route_dict = {
            'route_id': i,
            'route': route[0],
            'route_score': route[1]
        }
        result_routes_list.append(route_dict)
    return result_routes_list


if __name__ == '__main__':
    conf = Config('config/bionavi_conf.yaml')


    viz_dir = f"viz/"
    if not os.path.exists(viz_dir):
        os.makedirs(viz_dir)
    conf.viz_dir = viz_dir

    try:
        result = run(conf)
    except Exception as e:
        result = None
        print(e)

    if result is not None:
        for route in result:
            print(f"route id: {route['route_id']}\n"
                    f"route score:  {route['route_score']}\n"
                    f"route: {route['route']}\n")
