DDA-BERT

   
## Environment Setup

Create a new conda environment first:

```
conda create --name DDABert python=3.10
```

This will create an anaconda environment

Activate this environment by running:

```
conda activate DDABert
```

then install dependencies:

```
pip install -r ./requirements.txt
```

## Model Overview

### Step 1: Model Checkpoint Dir

/DDA_BERT/checkpoints/DDA_BERT.pt

### Step 2: Train DDA-BERT

deepspeed --bind_cores_to_rank   /DDA_BERT/train.py --deepspeed --deepspeed_config /DDA_BERT/ds_config.json --node_num 1 --gpu_num 8 --config /DDA_BERT/yaml/model.yaml


### Step 3: Eval DDA-BERT

cd /DDA_BERT; python eval.py --config /DDA_BERT/yaml/eval_model.yaml