# CMML_2

The FastCAR_Base.R script was sourced directly from the FastCAR repository (https://github.com/Nawijn-Group-Bioinformatics/FastCAR), since it is not distributed as an installable package and must be run standalone to access its functions.

# CMML_2 Ambient RNA Correction Project

This repository contains code and configurations for benchmarking FastCAR, DecontX and CellBender on the 10x PBMC4k dataset.

---

## Requirements

### R (FastCAR & DecontX)
- R 4.4.2  
- Packages: `Seurat`, `celda`, `scater`, `qlcMatrix`, `Matrix`, `pheatmap`, `ggplot2`, etc.

### CellBender Environment
- **OS**: Ubuntu 22.04  
- **Python**: 3.12  
- **PyTorch**: 2.3.1  
- **CUDA**: 12.1  
- **cuDNN**: 8  
- **NVCC** (CUDA compiler)  
- **VNC** (for remote desktop)

---

## CellBender Patch

Before running CellBender, apply the following changes to `cellbender/remove_background/checkpoint.py` as suggested in [GitHub Issue #296](https://github.com/broadinstitute/CellBender/issues/296):

```diff
-    torch.save(model_obj, filebase + '_model.torch')              # line 115
+    torch.save(model_obj.state_dict(), filebase + '_model.torch') # line 115
-    torch.save(scheduler, filebase + '_model.torch')              # line 116
+    scheduler.save(filebase + '_optim.torch')                     # line 116
