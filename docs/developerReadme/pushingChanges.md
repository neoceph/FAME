# How to push changes as a developer.

## Generating/Updating environment
### Windows ###
```powershell
conda env export | Select-String -NotMatch "^prefix: " | Out-File -Encoding utf8 environment.yaml
```

### Linux ###
```bash
conda env export | grep -v "^prefix: " > environment.yml
```
