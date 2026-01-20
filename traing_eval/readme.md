# DDA-BERT

## Environment Setup

1. Create a new conda environment:

```
conda create --name DDABert python=3.10
```

2. Activate this environment by running:

```
conda activate DDABert
```

3. Install dependencies:

```
pip install -r ./requirements.txt
```

## Training and Evaluation Pipeline

Model path: /DDA_BERT/software/resource/model/mp_rank_00_model_states.pt

### Step 1: Train DDA-BERT

deepspeed --bind_cores_to_rank   /DDA_BERT/traing_eval/train.py --deepspeed --deepspeed_config /DDA_BERT/traing_eval/ds_config.json --node_num 1 --gpu_num 8 --config /DDA_BERT/traing_eval/yaml/model.yaml


### Step 2: Evaluate DDA-BERT to determine the best model

cd /DDA_BERT/traing_eval; python eval.py --config /DDA_BERT/traing_eval/yaml/eval_model.yaml