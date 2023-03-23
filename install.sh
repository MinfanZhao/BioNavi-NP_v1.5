conda install -y rdkit=2020.03.3 -c rdkit 
conda install pytorch==1.10.0 torchvision==0.11.0 torchaudio==0.10.0 cudatoolkit=11.3 -c pytorch -c conda-forge
pip install networkx graphviz pyaml pynvml tqdm torchtext==0.6.0 configargparse rxnfp
pip install -e onmt/