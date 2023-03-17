conda install -y rdkit=2020.03.3 -c rdkit 
conda install -y pytorch torchvision torchaudio  cudatoolkit=11.1 -c pytorch-lts -c nvidia
pip install networkx graphviz pyaml pynvml tqdm torchtext==0.6.0 configargparse rxnfp
pip install -e onmt/