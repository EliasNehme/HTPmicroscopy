## Folder containing pre-trained model weights 
Here there should be 2 pre-trained model weights, one for the focus finder and one for the conditional 3D segmentor. The folder layout should look like this:
```
Models
    |___Focus_finder_TP_model_best_ckpt.pth
    |___Cellsnap_3D_segmentor_TP_model_best_ckpt.pth
    |___README.md
```
These weights are inside the compressed folder `Tetrapod_Pretrained_models_for_CellSnap.zip` which can be downloaded from Zenodo:
```
https://zenodo.org/records/10928122/files/Tetrapod_Pretrained_models_for_CellSnap.zip?download=1
```